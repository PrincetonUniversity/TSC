MODULE tsc_fmcfm_transport_mod
!
!   Module tsc_fmcfm_transport_mod provides an interface between the TSC code
! and the FMCFM Common Transport Interface.
!
!   The first part of this module consists of declarations of data
! structures needed by the FMCFM interface to transport models.
!
!   Routines included in this module are:
! tsc_transport : computes transport coefficients for the tsc code
! tsc_fmcfm_transport_allocate : allocates data structures
!                                   needed by the FMCFM interface
! tsc_fmcfm_transport   : computes transport coefficients using FMCFM
! tsc_trcdef_neoclass7  : neoclassical transport corresponding to option 7
!                           from file trcdef.f90
! tsc_trcdef_sawtooth   : sawtooth model as it was implemented in trcdef.f90
!
!   Record of revisions:
! Date        Programmer  Description of change
! Dec 2009    Bateman     Corrected dchii_dtip(j) and similar expressions
!                         Corrected fm_sp(j,jn)%*%gradTemp for Newton's method
!                         Corrected fm_eq(j,jn)%gradRho
! Jan 2009    Bateman     First draft of module
!
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7----:---

USE type_mod       ! Data needed tor transport calculations within FMCFM
USE anomalous_type_mod
USE flow_shear     ! To compute flow shear rate and other non-local arrays
USE glf23_mod

! Modules not yet used:

! USE allmodels_mod  ! All of the available transport models
! USE flow_shear     ! Flow shear rate computation
! USE mmm95_mod

IMPLICIT NONE

PRIVATE  ! The variables in this data structure are private so that
         ! they will be used only within this module

!-----------------------------------------------------------------------------

! INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

! In all of the following allocatable arrays, 
! the first index refers to the radial zone boundary,
! the second index refers to the number of times that the 
! transport model has to be called during each Newton's iteration 

! Equilibrium geometry in FMCFM
! Data structure MagGeom is in file type_mod.f90

TYPE(MagGeom), ALLOCATABLE      :: fm_eq(:,:)

! Species variables in FMCFM
! Data structure AllSpecies is in file type_mod.f90

TYPE(AllSpecies), ALLOCATABLE   :: fm_sp(:,:)

! Transport fluxes in FMCFM
! Data structure TransSpecies is in file type_mod.f90

TYPE(TransSpecies), ALLOCATABLE :: fm_flux(:,:)

! Transport coefficients in FMCFM
! Data structure TransSpeciesCoeff is in file type_mod.f90

TYPE(TransSpeciesCoeff), ALLOCATABLE :: fm_coef(:,:)

! Transport details in FMCFM
! Data structure AnomTransDetails is in file type_mod.f90

TYPE(AnomTransDetails), ALLOCATABLE :: fm_ad(:,:)

! Variables passed to transport models
! Data structure AnomSurfVars is in file type_mod.f90

TYPE(AnomSurfVars), ALLOCATABLE :: fm_var(:,:)

! Moved here for easy control  oct-15-2009
REAL*8 :: &
   zchimin  & ! minimimum value of diffusivity [m^2/sec]
 , zchimax  & ! maximimum value of diffusivity [m^2/sec] (chiitima, chietema, diffnema)
 , zrelfactor &  ! relaxation factor to deal with stiff equations
 , zrelfactor2   ! relaxation factor used in time smoothing of chiitima, chietema, diffnema

INTEGER :: &
   ipmax      ! outer radial index over which GLF23 is computed

INTEGER, SAVE :: &
   initialized = 0 &  !
 , istep = 1          ! local time step counter

! Parameters for the GLF23 model

TYPE( Glf23Flags ) :: glfflags

      REAL*8, ALLOCATABLE, DIMENSION(:) :: formf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: formi

PUBLIC tsc_transport &
 , tsc_fmcfm_transport_allocate &
 , tsc_fmcfm_transport

!-----------------------------------------------------------------------------

CONTAINS
!-----------------------------------------------------------------------------
SUBROUTINE tsc_transport ( k_model, kerror )

!  Sbrtn tsc_transport computes transport coefficients for the tsc code
!    using FMCFM as well as models that are not in FMCFM
!  This routine was patterned after .../tsc/tr_tsc/source/trcdef.f90
!
!  The output arrays are:
!
!  chiitima(j) ion thermal diffusivity chi_i at zone boundary j [cm^2/sec]
!  chietema(j) electron thermal diffusivity chi_e at zone boundary j [cm^2/sec]
!  diffnema(j) particle diffusivity D_i at zone boundary j [cm^2/sec]
!
!---------
!
!  Input variables from common blocks in clinam are

USE clinam            ! TSC common blocks in file module_clinam.f90
USE newplot           ! Allocatable arrays in file module_newplot.f90
USE saprop            ! Allocatable arrays in file module_saprop.f90

IMPLICIT NONE
INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

INTEGER, INTENT(IN)    :: k_model !Choice of transport model

! k_model = 0 for MMM95 transport model
! k_model = 26 for GLF23 transport model together with neoclassical
!              corresponding to model 6 in .../tsc/tr_tsc/source/trcdef.f90

INTEGER, INTENT(INOUT) :: &
   kerror   ! On input, kerror is the number of Newton's iterations requested
            ! On output, kerror is an error condition
            ! kerror = 0 for a normal return

!cj REAL*8 :: &
!cj    zrelfactor2 &     ! used in time smoothing of chiitima, chietema, diffnema
!cj  , zchimax           ! maximum allowed value of chiitima, chietema, diffnema

!cj INTEGER, SAVE :: &
!cj   istep = 1 &       ! time step index
!cj    i_newton_max = 5  ! number of times the transport model is called
                     ! for each iteration of Newton's method

INTEGER :: iprint = 0       ! controls amount of printout

INTEGER :: &
  j          ! Do loop index

REAL*8 :: &
  sum, sum2, sum3,sum3i,sum4,sum5,sum6,sum7,sum8,sum9 &
 ,cswitchre,enercn,ptot,ptoti,wdot,dvol,term,term2,term3,term4
 
!-----------------------------------------------------------------
! Set constants
itrmode = itrmod   !cj added oct-15-2009 to fix chie plot

ipmax = INT(pwidthc * REAL(npsit))

zchimin = 1.0e-8
zchimax = 50.E6_R8/ max(apl,1._R8)

zrelfactor = acoef(3010)
zrelfactor2 = acoef(3009)

!-----------------------------------------------------------------
! Initialize arrays as needed

IF ( kcycle <= 0 ) THEN

  chiitima = 1.0E-8_R8
  chietema = 1.0E-8_R8
  chiicopi = 0.0_R8
  chiecopi = 0.0_R8
  chiinca  = 0.0_R8
  chienca  = 0.0_R8

ELSE

! save previous diffusiviies

  chiitimao = chiitima
  chietemao = chietema
  diffnemao = diffnema

  chiitimao(1) = chiitimao(2)
  chietemao(1) = chietemao(2)
  diffnemao(1) = diffnemao(2)

!.....copied from tecdef to fix popcon plots, cj oct-23-09 
!.....calculate global integrals and form factors

      IF(.not.ALLOCATED(formf)) ALLOCATE( formf(ppsi), STAT=kerror)
      IF(.not.ALLOCATED(formi)) ALLOCATE( formi(ppsi), STAT=kerror)

      sum = 0._R8
      sum2 = 0._R8
      sum3 = 0._R8
      sum3i= 0._R8
      sum4 = 0._R8
      sum5 = 0._R8
      sum6 = 0._R8
      sum7 = 0._R8
      sum8 = 0._R8
      sum9 = 0._R8
      cswitchre = 0._R8
      formf(1) = 0._R8
      formi(1) = 0._R8
      do 300 j=2,npsit
!
!
      dvol = vary(j)-vary(j-1)
!
      term = (.5_R8*(as(j-1)+as(j)) +vlooph(j) )/tpi                     &
     &     *(gxmja2(j)-gxmja2(j-1))
      term2 = dvol*(savee(j) + savea(j)*ialpha                           &
!    1 - sradion(j)*usdp/usdt - savebre(j) - savecyc(j) - saveimp(j)
!.....above line commented out 10/31/96....scj
     & + savefw(j) + savebm(j)  + savelh(j))
      term3 = dvol*(savei(j) + savia(j)*ialpha + savibm(j) + savilh(j)   &
     & +savifw(j) )
      term4 =                                                            &
     &  dvol*(equila(j)*((1._R8+avez(j))*ade(j)-avez(j)*adp(j)))/vpg(j)
      if(term .le. 0) term=0._R8
!     if(term2.le. 0) term2=0.
!
      sum = sum + term
      sum2 = sum2 + 1.5_R8*dvol*adp(j)/vpg(j)
      sum3 = sum3 + term2  + term + term3
      sum3i = sum3i + term3 + term4
      sum4 = sum4 + te(j)*dvol
      sum5 = sum5 + ti(j)*dvol
      sum6 = sum6 + dvol
      sum7 = sum7 + dvol*(savea(j) + savia(j))
      sum8 = sum8 + dvol*(savei(j)+savee(j)+savifw(j)                    &
     & + savefw(j) + savebm(j) + savibm(j) + savelh(j) + savilh(j))
      sum9 = sum9 + dvol*(                                               &
     & - sradion(j)*usdp/usdt - savebre(j) - savecyc(j) - saveimp(j) )
!
  304 continue
      formf(j) = formf(j-1) + term + term2 + term3
      if(formf(j).lt.0.0_R8) formf(j) = 0._R8
      formi(j) = formi(j-1)  + term3 + term4
      if(formi(j).lt.0.0_R8) formi(j) = 0._R8
  300 continue
      if(sum4.gt.0.0_R8.and.sum5.gt.0.0_R8) go to 305
      ineg=43
      return
  305 pohmic = sum*udsp/udst
      if(pohmic.lt.0.0_R8) pohmic = 0.0_R8
      enercn = sum2*udsp
      teform = te(2)*sum6/sum4
      tiform = ti(2)*sum6/sum5
      palpha = sum7*ialpha*udsp/udst
      palphap = sum7*udsp/udst
      prad = sum9*udsp/udst
      paux = sum8*udsp/udst
!
      if(sum3.le.1.E-12_R8) sum3 = 1.E-12_R8
      if(sum3i.le.1.E-12_R8) sum3i = 1.E-12_R8
      ptot = sum3*udsp/udst
      ptoti = sum3i*udsp/udst
!
!.....NOTE:  use Coppi/Tang during ohmic phase for itrmod .ge.8
      itrmode = itrmod
      if((paux .eq. 0.0_R8 .and. palpha .lt. 10.E6_R8)                   &
     &                .and. itrmod.ge.8)       itrmode = 2
      if(itrmod.eq.12) itrmode = 2
      if(itrmod.eq.13) itrmode = 2
      if(itrmod.eq.14) itrmode = 2
!
!...==> special diagnostic
      if(iplt .ge. nskipl) then
      if(ialpha.gt.0) write(nterm,6661) paux,palpha,itrmod,itrmode
 6661 format("paux,palpha,itrmod,itrmode=",1p2e12.4,0p2i6)
      endif
!
!.....special fix for disruption modeling
      if(times .gt. acoef(95)) qsaw = acoef(96)
!
      wdot = wdotmw*1.E6_R8*usdp/usdt
      tauems = (sum2/(sum3-wdot))*udst*1000._R8

ENDIF  ! This was 309 continue  in trcdef.f90

!.....increase chi inside max tempsurface

!.....If ilhcd is not zero, a power dependent contribution to the
!.....resistivity, resulting from the hot plasma population, is
!.....calculated.  sighot is the ratio of sigma_hot obtained from
!.....Eq. 11 of Fisch, Phys. Fluids 28, 245 (1985), to sigma_spitzer
!.....obtained from the Hirshman et al. paper refered to below.
!.....sigma_spitzer is the Hirshman result for conductivity when trapping
!.....is absent.

!-----------------------------------------------------------------
!  Use FMCFM interface to compute anomalous transport
!
!..Skip update of transport coefficients if kcycle <= 0

IF ( kcycle > 0 ) THEN

!..Update GLF23 transport coefficients only if iskipsf == 1
  IF ( 1 == iskipsf ) THEN

!..When using GLF23 (itrmod == 26) and acoef(4958) > 0.0 
!  i_newton_max = 5 because  transport model is
!  called 5 times for each iteration of Newton's method

    IF ( iprint > 0 ) THEN
      WRITE(*,*) ' call tsc_fmcmf_transport with k_model = ',k_model
    ENDIF

    CALL tsc_fmcfm_transport ( k_model, kerror )

  ENDIF

!-----------------------------------------------------------------
! Time smoothing and limits on chiitima(j), chietema(j)
! and diffnema(j)

  DO j = 2, npsit-1
    chiitima(j) = ( 1.0_R8 - zrelfactor2 ) * chiitima(j) &
      + 0.5_R8 * zrelfactor2 * ( chiitimao(j+1) + chiitimao(j-1) )
    chietema(j) = ( 1.0_R8 - zrelfactor2 ) * chietema(j) &
      + 0.5_R8 * zrelfactor2 * ( chietemao(j+1) + chietemao(j-1) )
    diffnema(j) = ( 1.0_R8 - zrelfactor2 ) * diffnema(j) &
      + 0.5_R8 * zrelfactor2 * ( diffnemao(j+1) + diffnemao(j-1) )

    chiitima(j) = MAX ( 1.0E-8_R8, MIN ( zchimax, chiitima(j) ) )
    chietema(j) = MAX ( 1.0E-8_R8, MIN ( zchimax, chietema(j) ) )
    diffnema(j) = MAX ( 1.0E-8_R8, MIN ( zchimax, diffnema(j) ) )
  ENDDO

!-----------------------------------------------------------------
! Compute neoclassical transport
!   corresponding to option 7 in trcdef.f90

  CALL tsc_trcdef_neoclass7

ENDIF  ! end of updating transport coefficients when kcycle <= 0

!-----------------------------------------------------------------
! Call sawtooth model

CALL tsc_trcdef_sawtooth

!-----------------------------------------------------------------
! Diagnostic printout

IF ( iprint > 0 .AND. istep < 5 ) THEN

  WRITE (*,*) ' End of sbrtn tsc_transport, istep = ',istep

    WRITE(*,*)
    WRITE(*,*) '   j   rminora' &
      , '  ti(j)   te(j)'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, rminora(j) &
        , ti(j), te(j)
    ENDDO

    WRITE(*,*)
    WRITE(*,*) '   j   rminora' &
      , '  rminora(j)  chiitima(j)   chietema(j)   diffnema(j)'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, rminora(j) &
        , chiitima(j), chietema(j), diffnema(j)
    ENDDO

    WRITE(*,*)
    WRITE(*,*) '   j   rminora' &
      , '  dchii_dtip  dchii_dtep  dchie_dtep  dchie_dtip ' &
      , '  dte_dpsi    dti_dpsi'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, rminora(j) &
        , dchii_dtip(j), dchii_dtep(j), dchie_dtep(j), dchie_dtip(j) & 
        , dte_dpsi(j), dti_dpsi(j)
    ENDDO

ENDIF

RETURN
END SUBROUTINE tsc_transport
!-----------------------------------------------------------------------------
SUBROUTINE tsc_fmcfm_transport_allocate ( &
   k_radial_dim, k_species_dim, k_newton_max, kerror )

! The data structures needed by the FMCFM interface are allocated

IMPLICIT NONE
INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

INTEGER, INTENT(IN)    :: k_radial_dim  ! dimension of radial arrays

INTEGER, INTENT(IN)    :: k_species_dim ! dimension of spcies arrays

!cj INTEGER, INTENT(INOUT) :: k_newton_max  ! number of times that the
INTEGER, INTENT(IN) :: k_newton_max  ! number of times that the
  ! transport model is called during each Newton iteration
  ! Note: At the present time, either k_newton_max = 1 or k_newton_max = 5

INTEGER, INTENT(INOUT) :: kerror ! = 0 for normal return, error flag

kerror = 0
!cj k_newton_max = MAX ( 1, k_newton_max )
!cj IF ( k_newton_max > 1 ) k_newton_max = 5

! Allocate radial and species arrays

IF ( .NOT. ALLOCATED( fm_eq ) ) THEN
  ALLOCATE ( fm_eq( 1:k_radial_dim, 1:k_newton_max ) &
  , STAT = kerror )
ENDIF

IF ( .NOT. ALLOCATED( fm_sp ) ) THEN
  ALLOCATE ( fm_sp( 1:k_radial_dim, 1:k_newton_max ) &
  , STAT = kerror )
ENDIF

IF ( .NOT. ALLOCATED( fm_flux ) ) THEN
  ALLOCATE ( fm_flux( 1:k_radial_dim, 1:k_newton_max ) &
  , STAT = kerror )
ENDIF

IF ( .NOT. ALLOCATED( fm_coef ) ) THEN
  ALLOCATE ( fm_coef( 1:k_radial_dim, 1:k_newton_max ) &
  , STAT = kerror )
ENDIF

IF ( .NOT. ALLOCATED( fm_ad ) ) THEN
  ALLOCATE ( fm_ad( 1:k_radial_dim, 1:k_newton_max ) &
  , STAT = kerror )
ENDIF

IF ( .NOT. ALLOCATED( fm_var ) ) THEN
  ALLOCATE ( fm_var( 1:k_radial_dim, 1:k_newton_max ) &
  , STAT = kerror )
ENDIF

RETURN
END SUBROUTINE tsc_fmcfm_transport_allocate
!-----------------------------------------------------------------------------
!cj SUBROUTINE tsc_fmcfm_transport ( k_model, k_newton_max, kerror )
SUBROUTINE tsc_fmcfm_transport ( k_model, kerror )

!  Sbrtn tsc_fmcfm_transport computes transport coefficients using FMCFM
!  This routine was patterned after .../tsc/tr_tsc/source/trcdef.f90
!
!  The output arrays are:
!
!  chiitima(j) ion thermal diffusivity chi_i at zone boundary j [cm^2/sec]
!  chietema(j) electron thermal diffusivity chi_e at zone boundary j [cm^2/sec]
!  diffnema(j) particle diffusivity D_i at zone boundary j [cm^2/sec]
!
!  dchii_dtip(j) ! derivative of ion thermal diffusivity
!                  wrt ion temperature gradient, used for Newton's method
!  dchie_dtip(j) ! derivative of electron thermal diffusivity
!                  wrt ion temperature gradient, used for Newton's method
!  dchii_dtep(j) ! derivative of ion thermal diffusivity
!                  wrt electron temperature gradient, used for Newton's method
!  dchie_dtep(j) ! derivative of electron thermal diffusivity
!                  wrt electron temperature gradient, used for Newton's method
!  dti_dpsi(j)   ! ion temperature gradient
!  dte_dpsi(j)   ! electron temperature gradient
!
!---------
!
!  Input variables from common blocks in clinam are
!
!  kcycle   

USE clinam            ! TSC common blocks in file module_clinam.f90
USE saprop            ! Allocatable arrays in file module_saprop.f90

IMPLICIT NONE

#ifdef HAVE_MPI
!...  JCdebug april-03-2009 
include 'mpif.h'
INTEGER :: myPE, totPEs, ierr
INTEGER(kind=ispec) :: mpi_gr_world
integer stat(MPI_STATUS_SIZE)

#ifdef I_HAVE_MPI
!cj   !... JCdebug april-09-2009
INTEGER :: my_ipmax, my_start_index, my_end_index, counts
INTEGER, ALLOCATABLE :: array_index_start(:), array_index_end(:), array_index_counts(:)

INTEGER, ALLOCATABLE :: mp_ranks(:)
INTEGER :: mp_gr_gnode, mp_comm_gnode
#ifdef I_HAVE_MPI_DEBUG
REAL*8, ALLOCATABLE :: cjtmp(:,:)
#endif
#endif
#endif

INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

INTEGER, INTENT(IN)    :: k_model !Choice of transport model

! k_model = 0 for MMM95 transport model
! k_model = 26 for GLF23 transport model together with neoclassical
!              corresponding to model 6 in .../tsc/tr_tsc/source/trcdef.f90

!cj INTEGER, INTENT(INOUT) :: k_newton_max
INTEGER :: k_newton_max

! k_newton_max is the number of times that the transport model is
!  called for each iteration of Newton's method

INTEGER, INTENT(INOUT) :: &
   kerror   ! On input, kerror is the number of Newton's iterations requested
            ! On output, kerror is an error condition
            ! kerror = 0 for a normal return

REAL*8 :: &
   ztfluxb &      ! toroidal flux within plasma
 , ztflux  &      ! toroidal flux within flux surface boundary
 , zarho_exp &    ! toroidal flux at last closed flux surface (LCFS)
 , zfacs          ! 2 * normalized radius

REAL*8, PARAMETER :: &
   zepsilon = 1.D-34  ! small number

REAL*8, SAVE :: &
   zdelta_t  ! differential temperature used in Newton's computation

REAL*8, SAVE, DIMENSION(5) :: &
   zdte(5)  & ! increment added to electron temperature gradient
 , zdti(5)    ! increment added to ion temperature gradient

REAL*8 :: &
   zchiitim & ! ion thermal diffusivity [m^2/sec]
 , zchietem & ! electron thermal diffusivity [m^2/sec]
 , zdiffnem !cj & ! ion particle diffusivity [m^2/sec]
!cj  , zchimin  & ! minimimum value of diffusivity [m^2/sec]
!cj  , zchimax  & ! maximimum value of diffusivity [m^2/sec]
!cj  , zrelfactor &  ! relaxation factor to deal with stiff equations
!cj  , zrelfactor2   ! relaxation factor

REAL*8, DIMENSION(npsit) :: &
   zrho     & ! normalized minor radius variable
              ! sqrt ( \psi_tor / \psi_tor_edge )
 , zdrhodr    ! d rho / d rminor


!..Diffusivities for diagnostic output

REAL*8, DIMENSION(npsit,5) :: &
   zchii   & ! ion thermal diffusivities [m^2/sec]
 , zchie   & ! electron thermal diffusivities [m^2/sec]
 , zdiff     ! ion particle diffusivity [m^2/sec]

REAL*8 :: ztemp

INTEGER :: &
   ierror  &          ! local error flag
 , iprint = 0 &       ! controls amount of printout
 , iskip  = 100 &     ! controls frequency of printout
 , i_species_dim = 5  ! Dimension of number of species

! INTEGER, SAVE :: &
!    initialized = 0 &  !
!  , istep = 1          ! local time step counter

!cj INTEGER :: &
!cj    ipmax   &   ! outer radial index over which GLF23 is computed
!cj  , ipsihm      ! index of inner edge of pedestal

INTEGER :: &
   j       &   ! Radial do-loop index
 , jn          ! Newton's do-loop index

!-----------------------------------------------------------------
! Set constants

IF ( iprint > 0 ) THEN
  WRITE(*,*) ' sbrtn tsc_fmcfm_transport with k_newton_max = '&
    , k_newton_max
ENDIF

!cj zchimin = 1.0e-8
!cj zchimax = 50.E6_R8/ max(apl,1._R8)
!cj zrelfactor = acoef(3010)

!cj ipsihm = INT(pwidthc * REAL(npsit))
!cj ipmax  = ipsihm

!-----------------------------------------------------------------
! Set up Neton's iteration variables

k_newton_max = 1
!cj k_newton_max = MAX ( 1, k_newton_max )
!cj IF ( k_newton_max > 1 ) k_newton_max = 5

zdelta_t = -1.0_R8
IF ( acoef(4958) > 0.0 ) zdelta_t = acoef(4958)

dchie_dtip = 0.0_R8
dchii_dtip = 0.0_R8
dchie_dtep = 0.0_R8
dchii_dtep = 0.0_R8
dti_dpsi   = 0.0_R8
dte_dpsi   = 0.0_R8

zchii      = 0.0_R8
zchie      = 0.0_R8
zdiff      = 0.0_R8

IF ( zdelta_t > 0.0_R8 .and. acoef(4957) > 0.0_R8 ) THEN

  k_newton_max = 5

  zdte(1) = 0.0_R8
  zdte(2) = - zdelta_t
  zdte(3) =   zdelta_t
  zdte(4) = 0.0_R8
  zdte(5) = 0.0_R8

  zdti(1) = 0.0_R8
  zdti(2) = 0.0_R8
  zdti(3) = 0.0_R8
  zdti(4) = - zdelta_t
  zdti(5) =   zdelta_t

!cj ELSE
!cj 
!cj   k_newton_max = 1

ENDIF

IF ( iprint > 0 ) THEN
  WRITE(*,*) ' sbrtn tsc_fmcfm_transport with k_newton_max = '&
    , k_newton_max
ENDIF

write(*,*) "JCdebug tsc_fmcfm_transport: iskipsf=", iskipsf, &
           "newton_terms=", k_newton_max, &
           "npsi=", npsi, &
           "npsit=", npsit, &
           "ppsi=", ppsi, &
           "ipmax=", ipmax

!-----------------------------------------------------------------
! Allocate the arrays in the data structure in this module

CALL tsc_fmcfm_transport_allocate ( &
  npsi, i_species_dim, k_newton_max, kerror )

!-----------------------------------------------------------------
! Initialize switches used by GLF23 model

  ierror = 0
  IF ( initialized < 1) THEN

    CALL FmInitFlagsGlf23( glfflags, ierror )

  ENDIF

!..Set glfflags%basic%tranFlags(j) corresponding to itport_pt(j), j=1,5

glfflags%basic%tranFlags(1) = int(acoef(4981))
glfflags%basic%tranFlags(2) = int(acoef(4982))
glfflags%basic%tranFlags(3) = int(acoef(4983))
glfflags%basic%tranFlags(4) = int(acoef(4984))
glfflags%basic%tranFlags(5) = int(acoef(4985))

glfflags%basic%BtFlag  = 1  ! corresponds to ibt_flag

glfflags%basic%rotationFlag  =  int(acoef(4988))

glfflags%basic%calpha   = 0.0_R8  ! No alpha stabilization

glfflags%extra%xky0     = 0.2_R8

!..Bateman, 15 March 2009, try turning alpha_e = 0.0_R8

glfflags%extra%alpha_e  = 1.35_R8
glfflags%basic%cExB     = 1.35_R8

!-----------------------------------------------------------------
! Diagnostic output

IF ( iprint > 8 .AND. initialized < 1 ) THEN
WRITE (*,*) &
' Data from glfflags after FmInitFlagsGlf23 in tsc_fmcfm_transport_mod.f90'
WRITE (*,*) ' glfflags%basic%eigen  = ', glfflags%basic%eigen
WRITE (*,*) ' glfflags%basic%nroot  = ', glfflags%basic%nroot
WRITE (*,*) ' glfflags%basic%iglf   = ', glfflags%basic%iglf
WRITE (*,*) ' glfflags%basic%tranFlags  = ', glfflags%basic%tranFlags
WRITE (*,*) ' glfflags%basic%rotationFlag  = ', glfflags%basic%rotationFlag
WRITE (*,*) ' glfflags%basic%BtFlag  = ', glfflags%basic%BtFlag
WRITE (*,*) ' glfflags%basic%UnitNum   = ', glfflags%basic%UnitNum 
WRITE (*,*) ' glfflags%basic%ifElec  = ', glfflags%basic%ifElec
WRITE (*,*) ' glfflags%basic%ifETG  = ', glfflags%basic%ifETG
WRITE (*,*) ' glfflags%basic%cExB  = ', glfflags%basic%cExB
WRITE (*,*) ' glfflags%basic%calpha  = ', glfflags%basic%calpha
WRITE (*,*) ' glfflags%extra%ikymax  = ', glfflags%extra%ikymax
WRITE (*,*) ' glfflags%extra%xkymin  = ', glfflags%extra%xkymin
WRITE (*,*) ' glfflags%extra%xkymax  = ', glfflags%extra%xkymax
WRITE (*,*) ' glfflags%extra%xky0  = ', glfflags%extra%xky0
WRITE (*,*) ' glfflags%extra%cnorm  = ', glfflags%extra%cnorm
WRITE (*,*) ' glfflags%extra%adamp  = ', glfflags%extra%adamp
WRITE (*,*) ' glfflags%extra%park  = ', glfflags%extra%park
WRITE (*,*) ' glfflags%extra%alphaP = ', glfflags%extra%alphaP
WRITE (*,*) ' glfflags%extra%cbetae  = ', glfflags%extra%cbetae
WRITE (*,*) ' glfflags%extra%cnu  = ', glfflags%extra%cnu
WRITE (*,*) ' glfflags%extra%iflagin  = ', glfflags%extra%iflagin
WRITE (*,*) ' glfflags%extra%xparam  = ', glfflags%extra%xparam
WRITE (*,*) ' glfflags%extra%alpha_e  = ', glfflags%extra%alpha_e
WRITE (*,*)
WRITE (*,*) glfflags

ENDIF  ! end of diagnostic printout

IF ( iprint > 7 .AND. initialized < 1 ) THEN
  WRITE (*,*) &
  '  jmm, rlti_gf, rlne_gf, rlni_gf, rlnimp_gf, chiitim, rmin_gf'
ENDIF

!-----------------------------------------------------------------
! Initialize derived data types

ierror = 0
IF ( initialized < 1 ) THEN

  DO j = 1, npsi             ! Radial index

    DO jn = 1, k_newton_max  ! Newton's index

      CALL FmInitMagGeom ( fm_eq(j, jn), ierror )
        IF ( ierror > 0 ) &
          WRITE(*,*) ' After CALL FmInitMagGeom, ierror = ', ierror
      CALL FmInitAllSpecies( fm_sp(j, jn), 1, ierror )
        IF ( ierror > 0 ) &
          WRITE(*,*) ' After CALL FmInitAllSpecies, ierror = ', ierror
      CALL FmInitTransSpecies( fm_flux(j, jn), 1, ierror )
        IF ( ierror > 0 ) &
          WRITE(*,*) ' After CALL FmInitTransSpecies, ierror = ', ierror
      CALL FmInitAnomSurfVars( fm_var(j, jn), ierror )
        IF ( ierror > 0 ) &
          WRITE(*,*) ' After CALL FmInitAnomSurfVars, ierror = ', ierror

    ENDDO

  ENDDO

ENDIF

!-----------------------------------------------------------------
! Initialize arrays as needed

! Code has been moved to sbrtn tsc_transport

!-----------------------------------------------------------------
! Loop over index j for radial zone boundaries
! to define local arrays

DO j = 1, npsit

  zrho(j) = sqrt( REAL(j-1) / REAL( npsit-1 ) )

ENDDO

ztfluxb = (npsit-1)*dpsi  ! toroidal flux within plasma
zarho_exp = sqrt ( ztfluxb * xplas / ( pi * gzero ) )
     ! \f$arho=\sqrt{\frac{\Phi_{total}}{\pi B_T}}~[m]\f$

!-----------------------------------------------------------------
! Header for diagnostic printout

IF ( iprint > 0 .AND. initialized < 1 ) THEN

  WRITE(*,*) ' zarho_exp=', zarho_exp,' ztfluxb=',ztfluxb &
    , ' zepsilon=', zepsilon

  WRITE (*,*)
  WRITE (*,*) '  npsit = ', npsit
  WRITE (*,*) ' Diangostic output from tsc_fmcfm_transport' &
    ,' at time step istep = ', istep
!  WRITE (*,*) '  j    radius   chi_i    chi_e   diff_h'
!  WRITE (*,*) '  j   ztemp   gradDen  gradDen_imp  chii  rminor'

ENDIF

IF ( iprint > 0 ) THEN

  WRITE(*,*) &
   'j, rminora(j),   ti(j),    te(j)'

  DO j=1,ipmax
    WRITE(*,"(' 1# ',I5,9ES12.4)") &
      j, rminora(j), ti(j), te(j)
  ENDDO

ENDIF

!-----------------------------------------------------------------
! Loop over index j for radial zone boundaries
!  (was do 4000 j=1,npsit in trcdef.f90

DO j = 1, ipmax       ! Radial index

  DO jn = 1, k_newton_max ! Newton's index


!-----------------------------------------------------------------
! define surface averaged resistivity and equilibration

!......neoclassical conductivity from hirshman,hawryluk, nuclear fusion

! The following statement in trcdef.f90 is used to select 
! the combination of transport models that will be computed
!
!      go to(1,2,3,4,5,6,7,2,1,2,1),itrmode
!
!-----------------------------------------------------------------
! Set elements of the data structures for FMCFM

!..fm_eq is an element in data structure MagGeom in file type_mod.f90

  fm_eq(j,jn)%rmin       = rminora(j)  ! minor radius [m]
  fm_eq(j,jn)%rmaj       = rmajora(j)  ! major radius [m]
  fm_eq(j,jn)%kappa      = elonga(j)   ! elongation
  fm_eq(j,jn)%q          = qprof2(j)   ! magnetic q-value

  IF (j.EQ.1) THEN
    fm_eq(j,jn)%DlnQDlnRho &
      = j * ( qprof2(j+1) - qprof2(j) ) / qprof2(j)
  ELSE
    fm_eq(j,jn)%DlnQDlnRho &
      = j * ( qprof2(j+1) - qprof2(j-1) ) / qprof2(j)
  ENDIF
     ! d ln q / d ln rho
     ! This expression is based on trcdef.f90 and sglf.f
     ! in the TSC code
     ! there might be a problem at j = npsit

  zfacs   = 2.0_R8 * sqrt( REAL(j-1) / REAL(npsit-1) )

  fm_eq(j,jn)%arho       = zarho_exp
     ! \f$arho=\sqrt{\frac{\Phi_{total}}{\pi B_T}}~[m]\f$
  fm_eq(j,jn)%rho        = zrho(j)
     ! square root of normalized torodial flux
  fm_eq(j,jn)%gradRho &
    = 0.5_R8 * SQRT( gja2(j) * tpi * qprof2(j) * xplas &
      / ( vp2(j) * pi * gzero * (j-1) * dpsi ) ) / zarho_exp
     ! \f$< | \nabla \rho | > [1/m]\f$
     !     = ( zrho(j) - zrho(j-1) ) &
     !       / ( rminora(j) - rminora(j-1) + zepsilon )
     ! computed as it is computed in callglf2d
     ! but without the factor zarho_exp
  fm_eq(j,jn)%gradRhoSq  = fm_eq(j,jn)%gradRho * fm_eq(j,jn)%gradRho
     ! Warning: probably incorrect way it is defined in trcdef.f90
     ! Warning: At present this has the wrong units
     ! \f$< | \nabla \rho|^2 > [1/m^2]\f$
  fm_eq(j,jn)%rmajor     = xmag
     ! geometrical major radius at  magnetic axis \f[m]\f$
  fm_eq(j,jn)%Bt         = gzero / xplas
     ! device vaccum field, Bt = constant(a)/R(a)  \f$[T]\f$
  IF (j.EQ.1) THEN
    fm_eq(j,jn)%DrDrho    &
      =   ( rminora(j+1) - rminora(j) )        &
       / ( zrho(j+1)    - zrho(j) + zepsilon)
  ELSE
    fm_eq(j,jn)%DrDrho    &
      =   ( rminora(j) - rminora(j-1) )        &
       / ( zrho(j)    - zrho(j-1) + zepsilon)
  ENDIF

! fm_sp is an element in data structure AllSpecies in file type_mod.f90

!..Ion species in fm_sp(j)%ion(1)%...

  fm_sp(j,jn)%ion(1)%charge      = 1.0_R8
    ! charge state of isotope, (-1) for electrons
  fm_sp(j,jn)%ion(1)%nprotons    = 1.0_R8
    ! number of protons in nucleus, (0) for electrons
  fm_sp(j,jn)%ion(1)%amu         = amgas
    ! atomic mass number unit, (1/1827)~0 for electrons
  fm_sp(j,jn)%ion(1)%temperature &
    = 0.5_R8 * ( ti(j+1) * ( 1.0_R8 + zdti(jn) ) &
               + ti(j) * ( 1.0_R8 - zdti(jn) ) )
    ! ion temperature at zone boundary in eV
    ! Need to implement the following bound used in trcdef.f90
    ! if(ti_m(i) .le. 0.1_R8*te_m(i)) ti_m(i) = 0.1_R8*te_m(i)
  fm_sp(j,jn)%ion(1)%gradTemp    &
    = zfacs * ( npsit - 1 ) * ( ti(j+1) - ti(j) &
        + ( ti(j+1) + ti(j) ) * zdti(jn) )
    ! Ion temperature gradient w.r.t. \rho
  IF ( 1 == jn &
    .AND. fm_sp(j,jn)%ion(1)%gradTemp >= -zepsilon ) THEN
    fm_sp(j,jn)%ion(1)%gradTemp    &
      = 0.1_R8 * zfacs * ( npsit - 1 ) * ( te(j+1) - te(j) ) &
        * ( ti(j+1) + ti(j) ) / ( te(j+1) + te(j) )
  ENDIF
  fm_sp(j,jn)%ion(1)%density     &
    =  0.5_R8 * ( ane(j+1)/avez(j+1) + ane(j)/avez(j) )
    ! ion density [m^{-3}]
  fm_sp(j,jn)%ion(1)%gradDen     &
    = zfacs * ( npsit - 1 ) * (ane(j+1)/avez(j+1)-ane(j)/avez(j))
    ! Ion density gradient w.r.t. \rho
  fm_sp(j,jn)%ion(1)%vphi        = 0.0_R8
    ! toroidal velocity  \f$[m/s]\f$
  fm_sp(j,jn)%ion(1)%vpara       = 0.0_R8
    ! parallel velocity \f$[m/s]\f$
  fm_sp(j,jn)%ion(1)%vperp       = 0.0_R8
    ! perpendicular velocity \f$[m/s]\f$

!..Hydrogenic species in fm_sp(j)%AvgHydrogeni%...

  fm_sp(j,jn)%AvgHydrogenic%charge      = 1.0_R8
    ! charge state of isotope, (-1) for electrons
  fm_sp(j,jn)%AvgHydrogenic%nprotons    = 1.0_R8
    ! number of protons in nucleus, (0) for electrons
  fm_sp(j,jn)%AvgHydrogenic%amu         = amgas
    ! atomic mass number unit, (1/1827)~0 for electrons
  fm_sp(j,jn)%AvgHydrogenic%temperature &
    = 0.5_R8 * ( ti(j+1) * ( 1.0_R8 + zdti(jn) ) &
               + ti(j) * ( 1.0_R8 - zdti(jn) ) )
    ! ion temperature at zone boundary in eV
    ! Need to implement the following bound used in trcdef.f90
    ! if(ti_m(i) .le. 0.1_R8*te_m(i)) ti_m(i) = 0.1_R8*te_m(i)
  fm_sp(j,jn)%AvgHydrogenic%gradTemp    &
    = zfacs * ( npsit - 1 ) * ( ti(j+1) - ti(j) &
        + ( ti(j+1) + ti(j) ) * zdti(jn) )
    ! Ion temperature gradient w.r.t. \rho
  IF ( 1 == jn &
    .AND. fm_sp(j,jn)%AvgHydrogenic%gradTemp >= -zepsilon ) THEN
    fm_sp(j,jn)%AvgHydrogenic%gradTemp    &
      = 0.1_R8 * zfacs * ( npsit - 1 ) * ( te(j+1) - te(j) ) &
        * ( ti(j+1) + ti(j) ) / ( te(j+1) + te(j) )
  ENDIF
  fm_sp(j,jn)%AvgHydrogenic%density     &
    =  0.5_R8 * ( ane(j+1)/avez(j+1) + ane(j)/avez(j) )
    ! ion density [m^{-3}]
  fm_sp(j,jn)%AvgHydrogenic%gradDen     &
    = zfacs * ( npsit - 1 ) * (ane(j+1)/avez(j+1)-ane(j)/avez(j))
    ! Ion density gradient w.r.t. \rho
  fm_sp(j,jn)%AvgHydrogenic%vphi        = 0.0_R8
    ! toroidal velocity  \f$[m/s]\f$
  fm_sp(j,jn)%AvgHydrogenic%vpara       = 0.0_R8
    ! parallel velocity \f$[m/s]\f$
  fm_sp(j,jn)%AvgHydrogenic%vperp       = 0.0_R8
    ! perpendicular velocity \f$[m/s]\f$

  ztemp = - zfacs * npsit * (ane(j+1)/avez(j+1)-ane(j)/avez(j)) &
          / ( 0.5_R8 * ( ane(j+1)/avez(j+1) + ane(j)/avez(j) ) )

! Average impurity species in fm_sp(j)%AvgImpurity%...

  fm_sp(j,jn)%AvgImpurity%charge      = zimp
  fm_sp(j,jn)%AvgImpurity%speciesType = 'avgimpurity'
  fm_sp(j,jn)%AvgImpurity%amu         = 2*zimp
  fm_sp(j,jn)%AvgImpurity%density     = 1.0_R8
    ! WARNING: Must replace with computed value
  fm_sp(j,jn)%AvgImpurity%gradDen     = -1.0_R8
    ! Not sure (?)
  fm_sp(j,jn)%AvgImpurity%Temperature = 0.5_R8 * ( ti(j+1) + ti(j) )
  fm_sp(j,jn)%AvgImpurity%gradTemp    = 0.0_R8

  fm_sp(j,jn)%Zeff                    = 0.5_R8*(zeffa(j+1)+zeffa(j))

!..Electrons in fm_sp(j)%electron%...

  fm_sp(j,jn)%electron%density        = 0.5_R8*(ane(j+1)+ane(j))
  fm_sp(j,jn)%electron%gradDen        &
    = zfacs * ( npsit - 1 ) * (ane(j+1)-ane(j))
  fm_sp(j,jn)%electron%temperature &
    = 0.5_R8 * ( te(j+1) * ( 1.0_R8 + zdte(jn) ) &
               + te(j) * ( 1.0_R8 - zdte(jn) ) )
    ! electron temperature at zone boundary in eV
  fm_sp(j,jn)%electron%gradTemp       &
    = zfacs * ( npsit - 1 ) * (te(j+1) - te(j) &
        + ( te(j+1) + te(j) ) * zdte(jn) )

! Other factors in fm_var(j)%... in data structure AnomSurfVars
! Have to provide computed values

  fm_var(j,jn)%velAng    = 0.0_R8
  fm_var(j,jn)%wExb      = 0.0_R8
  fm_var(j,jn)%wExbDia   = 0.0_R8
  fm_var(j,jn)%wPara     = 0.0_R8
  fm_var(j,jn)%alphaMhd  = 0.0_R8

! Finalize species to calculate averages for impurities and fast ions

  ierror = 0
!  CALL FinalizeAllSpecies( fm_sp(j,jn), ierror )

!  WRITE (*,*) ' After CALL FinalizeAllSpecies, ierror = ', ierror
!-----------------------------------------------------------------
! 

  ENDDO  ! end of do-loop over Newton's index

ENDDO    ! end of do-loop over radial index
! (was 4000 continue in trcdef.f90)

!-----------------------------------------------------------------
! Compute flow shear rate and other non-local arrays

DO jn = 1, k_newton_max ! Newton's index

  ierror = 0
  CALL calcFlowShearGlf23( glfflags & ! Parameters for the GLF23 model
 , fm_eq(1:ipmax, jn) &  ! Equilibrium: MagGeom in file type_mod.f90
 , fm_sp(1:ipmax, jn) &  ! Species: AllSpecies in file type_mod.f90
 , fm_var(1:ipmax, jn) & ! AnomSurfVars in file type_mod.f90
 , ipmax & ! number of flux surfaces
 , ierror )

  IF ( ierror > 0 ) THEN
     WRITE (*,*) ' ierror = ', ierror &
       , '  after CALL calcFlowShearGlf23 in tsc_fmcmf_transport'
     STOP
  ENDIF

ENDDO

!-----------------------------------------------------------------
! Compute transport fluxes
! fm_flux is an element in data structure TransSpecies in file type_mod.f90

#ifdef HAVE_MPI
!cj   !...  JCdebug april-03-2009
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, totPEs, ierr )
      CALL MPI_COMM_GROUP(mpi_comm_world,mpi_gr_world,ierr)
      IF ( .NOT. ALLOCATED( mp_ranks ) ) THEN
        ALLOCATE ( mp_ranks( 1:totPEs ) , STAT = kerror )
      ENDIF
     mp_ranks(1) = myPE
     CALL MPI_GROUP_INCL(mpi_gr_world,1,mp_ranks(1:1),mp_gr_gnode,ierr)
     CALL MPI_COMM_CREATE(mpi_comm_world,mp_gr_gnode,mp_comm_gnode,ierr)
     glfFlags%mpiComm = mp_comm_gnode

!      glfFlags%radialComm = mpi_comm_world
!      glfFlags%wavenumComm = mpi_gr_world
!      glfFlags%mpiComm = mpi_comm_world

#ifdef I_HAVE_MPI
      my_ipmax = ipmax / totPEs
      my_start_index =  myPE   *my_ipmax + 1
      my_end_index   = (myPE+1)*my_ipmax
      if(myPE==0     )   my_start_index=2        ! the first processor
      if(myPE==totPEs-1) my_end_index  =ipmax    ! the last  processor
      IF ( .NOT. ALLOCATED( array_index_start ) ) THEN
        ALLOCATE ( array_index_start( 1:totPEs ) , STAT = kerror )
      ENDIF
      IF ( .NOT. ALLOCATED( array_index_end ) ) THEN
        ALLOCATE ( array_index_end( 1:totPEs ) , STAT = kerror )
      ENDIF
      IF ( .NOT. ALLOCATED( array_index_counts ) ) THEN
        ALLOCATE ( array_index_counts( 1:totPEs ) , STAT = kerror )
      ENDIF
      array_index_start(:)=0
      array_index_end(:)=0
      array_index_counts(:)=0
      array_index_start(myPE+1)=my_start_index
      array_index_end(myPE+1)=my_end_index
      call MPI_Allgather(array_index_start(myPE+1), 1, MPI_INTEGER, &
                         array_index_start, 1, MPI_INTEGER, mp_comm_gnode, ierr);
!                        array_index_start, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr);
      call MPI_Allgather(array_index_end(myPE+1), 1, MPI_INTEGER, &
                         array_index_end, 1, MPI_INTEGER, mp_comm_gnode, ierr);
!                        array_index_end, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr);
      DO j = 1, totPEs      ! Radial index
         !... Broadcast 1 int from process myPE to every process in the group.  
!        call MPI_Bcast(array_index_start(j),1,MPI_INTEGER, j-1,MPI_COMM_WORLD,ierr)
!        call MPI_Bcast(array_index_end(j),1,MPI_INTEGER, j-1,MPI_COMM_WORLD,ierr)

         array_index_counts(j) = array_index_end(j) - array_index_start(j) + 1
         IF ( initialized < 1) &
         write(*,*) myPE+1, "JCdebug tsc_fmcfm_transport array_index is ", &
                    array_index_start(j), array_index_end(j), array_index_counts(j)
      ENDDO    ! end of do-loop over total pe
!cj      call MPI_Bcast(fm_ad(array_index_start(j):array_index_end(j),jn)%tmomentum%phi , &
!cj                     counts,MPI_DOUBLE_PRECISION,j-1,MPI_COMM_WORLD,ierr)
!cj      call MPI_Bcast(fm_ad(array_index_start(j):array_index_end(j),jn)%tmomentum%par , &
!cj                     counts,MPI_DOUBLE_PRECISION,j-1,MPI_COMM_WORLD,ierr)
!cj      call MPI_Bcast(fm_ad(array_index_start(j):array_index_end(j),jn)%tmomentum%perp, &
!cj                     counts,MPI_DOUBLE_PRECISION,j-1,MPI_COMM_WORLD,ierr)
!cjint MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
!cj                  void *recvbuf, int recvcount, MPI_Datatype recvtype, 
!cj                  MPI_Comm comm ) 
!cjint MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
!cj                   void *recvbuf, int *recvcounts, int *displs, 
!cj                   MPI_Datatype recvtype, MPI_Comm comm )

#ifdef I_HAVE_MPI_DEBUG
      IF ( .NOT. ALLOCATED( cjtmp ) ) THEN
        ALLOCATE ( cjtmp( 1:npsit, 1:k_newton_max ) , STAT = kerror )
      ENDIF
      cjtmp(:,:)=0.
#endif 
#endif
#endif

#ifdef I_HAVE_MPI
DO j = my_start_index, my_end_index      ! Radial index
#else
DO j = 2, ipmax      ! Radial index
#endif 
  DO jn = 1, k_newton_max ! Newton's index

    ierror = 0

#ifdef I_HAVE_MPI_DEBUG
        cjtmp(j,jn)=myPE+1
#else
    CALL InitTransSpeciesCoeff( &
        fm_coef(j,jn) & ! Transport coefficients
      , 1             & ! 1 ion species here
      , ierror )

    IF ( ierror > 0 ) THEN
       WRITE (*,*) ' ierror = ', ierror &
         , '  after CALL InitTransSpeciesCoeff in tsc_fmcmf_transport'
!   ELSE
!      WRITE (*,*) 'JCdebug tsc_fmcfm_transport: after InitTransSpeciesCoeff'
    ENDIF

    CALL CalcAnomTranGlf23( &
        fm_eq(j,jn)   & ! Equilibrium geometry
      , fm_var(j,jn)  & ! Variables passed to transport models
      , fm_sp(j,jn)   & ! Species variables
      , fm_coef(j,jn) & ! Transport coefficients
      , fm_ad(j,jn)   & ! Transport details
      , glfflags      & ! Parameters for the GLF23 model
!#ifdef HAVE_MPI
!      , glfFlags%radialComm    & ! cj added for mpi
!      , glfFlags%wavenumComm   & ! cj added for mpi
!#endif
      , ierror )

    IF ( ierror > 0 ) THEN
       WRITE (*,*) ' ierror = ', ierror &
         , '  after CALL CalcAnomTranGlf23 in tsc_fmcmf_transport'
!   ELSE
!      WRITE (*,*) 'JCdebug tsc_fmcfm_transport: after CalcAnomTranGlf23'
    ENDIF 
#endif

  ENDDO  ! end of do-loop over Newton's index 
ENDDO    ! end of do-loop over radial index

#ifdef I_HAVE_MPI_DEBUG
  DO j = 1, ipmax
  DO jn = 1, k_newton_max
     write(*,*) myPE+1, &
             "JCdebug tsc_fmcfm_transport before: cjtmp(", j, ",", jn, ")=", cjtmp(j,jn)
  ENDDO
  ENDDO

  DO jn = 1, k_newton_max ! Newton's index
  !cj...  DO j = 1, totPEs      ! total pe
  !cj...  counts = array_index_end(j) - array_index_start(j) + 1
  !cj...  call MPI_Bcast(cjtmp(array_index_start(j):array_index_end(j),jn), &
  !cj...                 array_index_counts(j),MPI_DOUBLE_PRECISION,j-1,MPI_COMM_WORLD,ierr)
  !cj...  ENDDO ! end of do-loop over total pe

     call MPI_Allgatherv(cjtmp(array_index_start(myPE+1),jn), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                         cjtmp, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                         mp_comm_gnode,ierr) 
!                        MPI_COMM_WORLD,ierr) 
  ENDDO  ! end of do-loop over Newton's index

  DO j = 1, ipmax
  DO jn = 1, k_newton_max
     write(*,*) myPE+1, &
               "JCdebug tsc_fmcfm_transport after: cjtmp(", j, ",", jn, ")=", cjtmp(j,jn)
  ENDDO
  ENDDO
  DEALLOCATE(cjtmp)

  write(*,*) myPE+1, "JCdebug tsc_fmcfm_transport after: MPI_Allgatherv"
  CALL MPI_FINALIZE(ierr)
  stop
#endif

go to 555
   jn=1
#ifdef HAVE_MPI
   DO j = my_start_index, my_end_index      ! Radial index
!cjif(myPE==totPEs-1) then
!cj   DO j = my_end_index, my_end_index      ! Radial index
#else
   DO j = 2, ipmax       ! Radial index
!cj   DO j = ipmax, ipmax       ! Radial index
#endif
      WRITE(*,*) j, k_newton_max, jn                 &
        , fm_eq(j,jn)%gradRho                        &
        , fm_eq(j,jn)%rmin                           &
        , fm_coef(j,jn)%ion(1)%thermal%diffusivity   &
        , fm_coef(j,jn)%electron%thermal%diffusivity &
        , fm_coef(j,jn)%ion(1)%particle%diffusivity  
   ENDDO    ! end of do-loop over radial index
 2222 format("JCdebug tsc_fmcfm_transport: j=", I4, x, 5p4e12.4)
#ifdef HAVE_MPI
!cjendif
#endif
555 continue
!cj   stop

!-----------------------------------------------------------------
! Diagnostic output

IF ( iprint > 6 .AND. initialized < 1 ) THEN

  DO j = 2, npsit-1       ! Radial index

    DO jn = 1, k_newton_max ! Newton's index

      zfacs   = 2.0_R8 * sqrt( REAL(j-1) / REAL(npsit-1) )

      ztemp = - zfacs * npsit * (ane(j+1)/avez(j+1)-ane(j)/avez(j)) &
          / ( 0.5_R8 * ( ane(j+1)/avez(j+1) + ane(j)/avez(j) ) )

      WRITE(*,"(' 3# ',I5,9ES12.4)") &
        j &
      , ztemp &
      , fm_sp(j,jn)%AvgHydrogenic%gradDen &
      , fm_sp(j,jn)%AvgHydrogenic%density &
      , zfacs &
      , fm_eq(j,jn)%rmin 

    ENDDO  ! end of do-loop over Newton's index

  ENDDO    ! end of do-loop over radial index

ENDIF


!-----------------------------------------------------------------
! Compute diffusivities and their derivatives that are passed to
! the rest of the TSC code

IF ( iprint > 0 .AND. initialized < 1 ) THEN

  WRITE(*,*)
  WRITE(*,*) ' k_newton_max = ',  k_newton_max
  WRITE(*,*) ' npsit = ', npsit
  WRITE(*,*) ' dpsi  = ', dpsi
  WRITE(*,*) ' ipmax = ', ipmax
  WRITE(*,*) ' zdelta_t = ', zdelta_t
  WRITE(*,*) ' zdti(j) = ', (zdti(j), j=1,5)
  WRITE(*,*) ' zdte(j) = ', (zdte(j), j=1,5)
  WRITE(*,*) ' 2#  j  zchii(j,2), zchii(j,3), zfacs ' &
    , ' gradTemp(j,2)  gradTemp(j,3) dchii_dtip(j)'

ENDIF

#ifdef I_HAVE_MPI
DO j = my_start_index, my_end_index      ! Radial index
#else
DO j = 2, ipmax       ! Radial index
#endif 
  DO jn = 1, k_newton_max ! Newton's index

!..Transfer diffusivities from FMCFM data structures
!  to local arrays for processing

!    zdrhodr(j) = ( zrho(j) - zrho(j-1) ) &
!      / ( rminora(j) - rminora(j-1) + zepsilon )

    zdrhodr(j) = 0.5_R8 * SQRT( gja2(j) * tpi * qprof2(j) * xplas &
         / ( vp2(j) * pi * gzero * (j-1) * dpsi ) )

    zchii(j,jn) = fm_coef(j,jn)%ion(1)%thermal%diffusivity
    zchie(j,jn) = fm_coef(j,jn)%electron%thermal%diffusivity
    zdiff(j,jn) = fm_coef(j,jn)%ion(1)%particle%diffusivity

!..Put limits on diffusivities

    zchii(j,jn) = MAX ( zchimin, MIN ( zchimax, zchii(j,jn) ) )
    zchie(j,jn) = MAX ( zchimin, MIN ( zchimax, zchie(j,jn) ) )
    zdiff(j,jn) = MAX ( zchimin, MIN ( zchimax, zdiff(j,jn) ) )

  ENDDO ! end of loop over Newton's index

!..apply relaxation factor to deal with stiff equations

  chiitima(j) = ( 1.0_R8 - zrelfactor ) * chiitima(j) &
                + zrelfactor * zchii(j,1)
  chietema(j) = ( 1.0_R8 - zrelfactor ) * chietema(j) &
                + zrelfactor * zchie(j,1)
  diffnema(j) = ( 1.0_R8 - zrelfactor ) * diffnema(j) &
                + zrelfactor * zdiff(j,1)

!..Compute derivatives needed for Newton's method

  IF ( 5 == k_newton_max ) THEN

  zfacs   = 2.0_R8 * SQRT( REAL(j-1) / REAL(npsit-1) )

  dchie_dtep(j) = dpsi * ( zchie(j,2) - zchie(j,3) ) &
                  * ( zfacs * ( npsit - 1 ) &
    / ( ( fm_sp(j,2)%electron%gradTemp / fm_sp(j,2)%electron%temperature &
        - fm_sp(j,3)%electron%gradTemp / fm_sp(j,3)%electron%temperature ) &
          * 0.5_R8 * ( te(j+1) + te(j) ) ) )

  dchii_dtep(j) = dpsi * ( zchii(j,2) - zchii(j,3) ) &
                  * ( zfacs * ( npsit - 1 ) &
    / ( ( fm_sp(j,2)%electron%gradTemp / fm_sp(j,2)%electron%temperature &
        - fm_sp(j,3)%electron%gradTemp / fm_sp(j,3)%electron%temperature ) &
          * 0.5_R8 * ( te(j+1) + te(j) ) ) )

  dchie_dtip(j) = dpsi * ( zchie(j,4) - zchie(j,5) ) &
                  * ( zfacs * ( npsit - 1 ) &
    / ( ( fm_sp(j,4)%ion(1)%gradTemp / fm_sp(j,4)%ion(1)%temperature &
        - fm_sp(j,5)%ion(1)%gradTemp / fm_sp(j,5)%ion(1)%temperature ) &
          * 0.5_R8 * ( ti(j+1) + ti(j) ) ) )

  dchii_dtip(j) = dpsi * ( zchii(j,4) - zchii(j,5) ) &
    * ( zfacs * ( npsit - 1 ) &
    / ( ( fm_sp(j,4)%ion(1)%gradTemp / fm_sp(j,4)%ion(1)%temperature &
        - fm_sp(j,5)%ion(1)%gradTemp / fm_sp(j,5)%ion(1)%temperature ) &
          * 0.5_R8 * ( ti(j+1) + ti(j) ) ) )

  dti_dpsi(j) = ( ti(j+1) - ti(j) ) / dpsi

  dte_dpsi(j) = ( te(j+1) - te(j) ) / dpsi

    IF ( iprint > 1 ) THEN

WRITE(*,*) &
' j  dti_dpsi(j)   ti(j+1)   ti(j)   dpsi'
WRITE(*,"(I5,9ES12.4)") &
  j, dti_dpsi(j), ti(j+1), ti(j), dpsi

      WRITE(*,"(' 2# ',I5,9ES12.4)") &
        j, zchii(j,2), zchii(j,3), zfacs &
        , fm_sp(j,4)%ion(1)%gradTemp &
        , fm_sp(j,5)%ion(1)%gradTemp &
        , dchii_dtip(j)

    ENDIF

  ENDIF

ENDDO    ! end of do-loop over radial index

#ifdef I_HAVE_MPI
  IF ( iprint > 0 ) &
  write(*,*) myPE+1, "JCdebug tsc_fmcfm_transport start mpi_allgatherv" 
  call MPI_Allgatherv(chiitima(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      chiitima, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  call MPI_Allgatherv(chietema(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      chietema, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  call MPI_Allgatherv(diffnema(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      diffnema, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  call MPI_Allgatherv(dchii_dtip(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      dchii_dtip, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  call MPI_Allgatherv(dchie_dtip(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      dchie_dtip, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  call MPI_Allgatherv(dchii_dtep(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      dchii_dtep, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  call MPI_Allgatherv(dchie_dtep(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      dchie_dtep, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  call MPI_Allgatherv(dti_dpsi(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      dti_dpsi, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  call MPI_Allgatherv(dte_dpsi(array_index_start(myPE+1)), array_index_counts(myPE+1), MPI_DOUBLE_PRECISION, &
                      dte_dpsi, array_index_counts, array_index_start-1, MPI_DOUBLE_PRECISION, &
                      mp_comm_gnode,ierr)
!                     MPI_COMM_WORLD,ierr)
  DEALLOCATE(array_index_start, array_index_end, array_index_counts, mp_ranks)
  IF ( iprint > 0 ) &
  write(*,*) myPE+1, "JCdebug tsc_fmcfm_transport end mpi_allgatherv" 
#endif


!-----------------------------------------------------------------
! Diffusivities outside of the radial range from j = 2, ipmax

chiitima(1) = 1.E-8_R8
chietema(1) = 1.E-8_R8
diffnema(1) = 1.E-8_R8

IF ( ipmax < npsit ) THEN

  DO j = ipmax+1, npsit

    chiitima(j) = 1.E-8_R8
    chietema(j) = 1.E-8_R8
    diffnema(j) = 1.E-8_R8

  ENDDO

ENDIF

!-----------------------------------------------------------------
! pedestal modification to GLF23

DO j = 2, npsit  ! loop over radial index

  ztflux = REAL(j-1) * dpsi  ! toroidal flux within zone boundary
!  ipsihm = INT( pwidthc * REAL(npsit) ) ! was set earlier

  IF ( 0 /= chiped .AND. 0 /= pwidthc &
       .AND. ztflux >= ipmax * dpsi ) THEN

    chiitima(j) = chiped * (ztflux/(FLOAT(ipmax-1)*dpsi))**4
    chietema(j) = chiped * (ztflux/(float(ipmax-1)*dpsi))**4

  ENDIF

ENDDO

!-----------------------------------------------------------------
! Apply smoothing in time, controlled by acoef(3009)
! and bounds on diffusivities

!cj oct-26-2009 fix segmentation bug DO j = 1, npsit  ! loop over radial index
DO j = 2, npsit-1  ! loop over radial index

!cj   zrelfactor2 = acoef(3009)

      chiitima(j) = (1._R8-zrelfactor2)*chiitima(j)             &  
     &      + .5_R8*zrelfactor2*(chiitimao(j+1)+chiitimao(j-1))
      chietema(j) = (1._R8-zrelfactor2)*chietema(j)             &  
     &      + .5_R8*zrelfactor2*(chietemao(j+1)+chietemao(j-1))
      diffnema(j) = (1._R8-zrelfactor2)*diffnema(j)             &  
     &      + .5_R8*zrelfactor2*(diffnemao(j+1)+diffnemao(j-1))

  chiitima(j) = MAX ( 1.0e-9_R8, MIN ( zchimax, chiitima(j) ) )
  chietema(j) = MAX ( 1.0e-9_R8, MIN ( zchimax, chietema(j) ) )
  diffnema(j) = MAX ( 1.0e-9_R8, MIN ( zchimax, diffnema(j) ) )

ENDDO

!-----------------------------------------------------------------
! Diagnostic printout

IF ( iprint > 0 .AND. initialized < 1 ) THEN

  jn = 1

  WRITE(*,*) ' zrelfactor = ', zrelfactor
  WRITE(*,*) ' zrelfactor2 = ', zrelfactor2

  IF ( k_newton_max > 1 ) THEN
    WRITE(*,*)
    WRITE(*,*) ' In tsc_fmcfm_transport_mod, istep = ', istep
    WRITE(*,*)
     WRITE(*,*) '   j   zchii(j,1)    zchii(j,2)   zchii(j,3)' &
      , '   zchii(j,4)   zchii(j,5)'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, (zchii(j,jn),jn=1,5)
    ENDDO

  ENDIF

  IF ( k_newton_max > 1 ) THEN
    WRITE(*,*)
    WRITE(*,*) ' In tsc_fmcfm_transport_mod, istep = ', istep
    WRITE(*,*)
    WRITE(*,*) '   j   dchii_dtip  dchii_dtep  dchie_dtep  dchie_dtip ' &
      , '  dte_dpsi    dti_dpsi'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, dchii_dtip(j), dchii_dtep(j), dchie_dtep(j), dchie_dtip(j) & 
        , dte_dpsi(j), dti_dpsi(j)
    ENDDO

  ENDIF

  WRITE(*,*) ' 1#  j      chi_i     chi_e     diff_i      rminor'
  DO j=1,npsit
    WRITE(*,"(' 1# ',I5,9ES12.4)") &
      j, zchii(j,1), zchie(j,1), zdiff(j,1) &
      , fm_eq(j,1)%rmin
  ENDDO

  WRITE(*,*) ' 2#  j  chiitima(j)  chietema(j)  diffnema(j)  rminor'
  DO j=1,npsit
    WRITE(*,"(' 2# ',I5,9ES12.4)") &
      j, chiitima(j), chietema(j), diffnema(j) &
      , fm_eq(j,1)%rmin
  ENDDO

ENDIF

!-----------------------------------------------------------------
! Increment local increment time step counter

initialized = 1    ! no further initializations
istep = istep + 1  ! increment time step counter
#ifdef I_HAVE_MPI_DEBUG
  IF ( iprint > 0 ) &
  write(*,*) myPE+1, "JCdebug tsc_fmcfm_transport end" 
#endif 

go to 556
#ifdef HAVE_MPI
!cj   call MPI_Barrier (MPI_COMM_WORLD,ierr) 
!cj   DO j = my_start_index, my_end_index      ! Radial index
if(myPE==totPEs-1) then
   DO j = my_end_index, my_end_index      ! Radial index
#else
!cj   DO j = 2, ipmax       ! Radial index
   DO j = ipmax, ipmax       ! Radial index
#endif
      WRITE(*,*) j      &
        , chiitima(j)   &
        , chietema(j)   &
        , diffnema(j)   &
        , dchii_dtip(j) &
        , dchie_dtip(j) &
        , dchii_dtep(j) &
        , dchie_dtep(j) &
        , dti_dpsi(j)   &
        , dte_dpsi(j)
   ENDDO    ! end of do-loop over radial index
#ifdef HAVE_MPI
endif
#endif
!  ion      thermal diffusivity chi_i at zone boundary j [cm^2/sec]
!  electron thermal diffusivity chi_e at zone boundary j [cm^2/sec]
!  particle diffusivity D_i at zone boundary j [cm^2/sec]
!  derivative of ion      thermal diffusivity
!            wrt ion temperature gradient, used for Newton's method
!  derivative of electron thermal diffusivity
!            wrt ion temperature gradient, used for Newton's method
!  derivative of ion      thermal diffusivity
!            wrt electron temperature gradient, used for Newton's method
!  derivative of electron thermal diffusivity
!            wrt electron temperature gradient, used for Newton's method
!  ion      temperature gradient
!  electron temperature gradient
556 continue

!cjgo to 557
!cjDO j = ipmax/2, ipmax/2+1     ! Radial index
!cj  DO jn = 1, 1 !!k_newton_max ! Newton's index
!cj  WRITE(45,FMT=*) "kcycle  j  jn = ", kcycle, j, jn
!cj  WRITE(46,FMT=*) "kcycle  j  jn = ", kcycle, j, jn
!cj  WRITE(47,FMT=*) "kcycle  j  jn = ", kcycle, j, jn
!cj  WRITE(48,FMT=*) "kcycle  j  jn = ", kcycle, j, jn
!cj  WRITE(49,FMT=*) "kcycle  j  jn = ", kcycle, j, jn
!cj  WRITE(50,FMT=*) "kcycle  j  jn = ", kcycle, j, jn
!cj  WRITE(51,FMT=*) "kcycle  j  jn = ", kcycle, j, jn
!cj  call FmDumpFlagsGlf23(glfflags,45,ierror)
!cj  call FmDumpMagGeom(fm_eq(j,jn),46,ierror)
!cj  call FmDumpTransSpeciesCoeff(fm_coef(j,jn),47,ierror)
!cj  call FmDumpTransSpeciesFlux(fm_flux(j,jn),48,ierror)
!cj  call FmDumpAllSpecies(fm_sp(j,jn),49,ierror)
!cj  call FmDumpAnomSurfVars(fm_var(j,jn),50,ierror)
!cj  call FmDumpAnomTransDetails(fm_ad(j,jn),51,ierror)
!cj  ENDDO ! end of loop over Newton's index
!cjENDDO    ! end of do-loop over radial index
!cj557 continue

RETURN
END SUBROUTINE tsc_fmcfm_transport
!-----------------------------------------------------------------
SUBROUTINE tsc_trcdef_neoclass7
!
! Implement the neoclassical transport model as it is implemented
! in option 7 from file
! .../tsc/tr_tsc/source/trcdef.f90
!-----------------------------------------------------------------

USE clinam            ! TSC common blocks in file module_clinam.f90
USE saprop            ! Allocatable arrays in file module_saprop.f90
USE newplot           ! Allocatable arrays in file module_newplot.f90

IMPLICIT NONE
INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

!  Local variables from trcdef.f90

REAL*8 :: zfp2, zbtmax, zfcyc, zfcyc2, zfcyci, zfcyci2, zfpi2 &
 , zy2x, zy2rat

REAL*8 :: ztemid, ztemidd, ztimid, ztimidd, ztimide, zanemid &
 , zetpera, zetperac, zallam, zalstar, zalinc, zgps, ztflux &
 , zcr, zsi, zalame, zdelt, zrfac, zw2lim, zdiff, zwratf, zpnorm  &
 , zpsmid, zaimh, zemid, zpmid, zaccn, zface, zfaci  &
 , zchienc, zchiinc, zchiesecj, zchiisecj !cj &
!cj  , zchimin

REAL*8 :: &
   zcoef1, zcoef2, zcoef3, zcoef4, zcoef5, zcoef6, zcoef3nc &
 , zgval, zgpval, zgppval, zbsq, zal0, zrmid

INTEGER :: &
   j       !cj &  ! Do loop index
!cj  , ipsihm

!-----------------------------------------------------------------
! 
!cj       zchimin = 1.0e-8
!cj       ipsihm = int(pwidthc * REAL(npsit))

      IF (ilhcd .GT. 0) THEN
        zfp2 = 806.2_R8* r0 / freqlh**2
        zbtmax = ABS(gs(imag,jmag))/xmag
        zfcyc = 28.0_R8 * zbtmax / freqlh
        zfcyc2 = zfcyc * zfcyc
        zfcyci = zfcyc * zion / (aion * 1836._R8)
        zfcyci2 = zfcyci * zfcyci
        zfpi2 = zfp2 * zion**2 / (aion * 1836._R8)
        zy2x = 1.0_R8/ (zfcyci * zfcyc)
        zy2rat = zy2x / (1.0_R8- zy2x)
        IF ( zy2x .GE. 1.0_R8.OR. zfpi2 .LT. zy2rat ) THEN
          zaccn = SQRT(zfp2)/zfcyc + SQRT(1._R8+ zfp2/zfcyc2 - zfpi2)
        ELSE
          zaccn = SQRT( 1._R8/(1._R8- zy2x) )
        ENDIF
      ENDIF

DO j=2,npsit

!.....calculate neoclassical transport coefficients
!     REF: Hirshman & Jardin PF 22 (1979) p 741

!.....define surface averaged resistivity and equilibration

      ztemid = .5_R8*(te(j+1)+te(j))
      zanemid = .5_R8*(ane(j+1)+ane(j))
      ztimid = .5_R8*(ti(j+1)+ti(j))
      IF ( ztimid.le.0.1_R8*ztemid) ztimid = 0.1_R8*ztemid
!      zanimid = .5_R8*(anhy(j+1)+anhy(j))*1.E-6_R8

      zallam = 24._R8-log(sqrt(zanemid*1.E-6_R8)/ztemid)
      zetperac= 1.03E-4_R8*zallam*ztemid**(-1.5_R8)*usdr
!
!.....add neoclassical corrections
!
!......neoclassical conductivity from hirshman,hawryluk, nuclear fusion
!......                               17 (page 611) 1977
!
      zcr = (0.56_R8/zeff)*(3.0_R8-zeff)/(3.0_R8+zeff)
      zsi = 0.58_R8+ 0.20_R8*zeff
      zalame = (3.40_R8/zeff)*(1.13_R8+zeff)/(2.67_R8+zeff)
      zdelt = sqrt(vary(j)/(tpi*pi*xmag))/xmag
!
      anue(j) = 0._R8
      if(ftrap(j) .gt. 0.0_R8)                                           &  
     &  anue(j) = (xmag*qprof2(j)*zanemid*zallam)                          &  
     &          / (ftrap(j)*zdelt*ztemid**2*10.2E16_R8)
      zrfac = zalame*(1.0_R8-ftrap(j)/(1.0_R8+zsi*anue(j)))                 &  
     &     *       (1.0_R8-zcr*ftrap(j)/(1.0_R8+zsi*anue(j)) )

      zetpera = zetperac/zrfac * acoef(810)
      sighot(j+1)=0._R8

!.....If ilhcd is not zero, a power dependent contribution to the
!.....resistivity, resulting from the hot plasma population, is
!.....calculated.  sighot is the ratio of sigma_hot obtained from
!.....Eq. 11 of Fisch, Phys. Fluids 28, 245 (1985), to sigma_spitzer
!.....obtained from the Hirshman et al. paper refered to below.
!.....sigma_spitzer is the Hirshman result for conductivity when trapping
!.....is absent.
!
      IF (ilhcd .GT. 0 .and. ifk.eq.0) THEN
!.......calculate w2lim based on lh accessibility condition.
!
!.......Normalize w2lim by v-thermal
        zw2lim = 7.15E2_R8/(zaccn*SQRT(ztemid) )
!.......Since raytracing results show that zw2lim become no larger than 10
!.......when zw2lim is greater than 10 we replace it by 10.  Larger values
!.......of zw2lim are likely to correspond to a runaway region.
        IF (zw2lim .GT. 10.0_R8) zw2lim = 10.0_R8
!.....Check that zw2lim is no less than w1lim + 1.0e-3
        zdiff = ABS(w1lim)+ 1.0E-3_R8
        IF (zw2lim .LT. zdiff) zw2lim = zdiff
        IF (w1lim .LT. 0.0_R8) zw2lim = -zw2lim
        zwratf = (zw2lim**4 - w1lim**4) / (4.0_R8*log(zw2lim/w1lim))
        zpnorm = 260.2_R8* savelh(j+1)*(udsp/udst)*SQRT(ztemid)/zanemid**2
        sighot(j+1)=( 8.0_R8/ (7.52_R8*zalame *(zeff+3.0_R8)*zallam*       &
     &  9.11E-28_R8) )*                                                  &  
     &              zpnorm * zwratf
      ENDIF

      zetpera = min( etav , zetpera / (1.0_R8+ sighot(j+1) ))

      if(ztemid.lt.30.0_R8.and. kcycle.gt.0)  then
        etpara(j) = 0.10_R8*(zetpera - etpara(j)) + zetpera
      else
        etpara(j) = zetpera
      endif
      if(whalos.gt.0 .and. etpara(j).gt.etah) etpara(j)=etah
!
      equila(j) = zanemid*3.1E-11_R8*zetperac*zeffa2(j)*udsr/usdt

      zrmid = .5_R8*(adn(j+1)/vp(j+1)+adn(j)/vp(j))
      zemid = .5_R8*(ade(j+1)/vpg(j+1)+ade(j)/vpg(j))
      zpmid = .5_R8*(adp(j+1)/vpg(j+1)+adp(j)/vpg(j))
      ztemidd = zemid / zrmid
      ztimidd = (zpmid-zemid) / zrmid
      if(ztimidd.le.0.0_R8) ztimidd = 0.1_R8*ztemidd

      zface = zemid * acoef(897)
      zfaci = (zpmid-zemid) * acoef(897)

! xsv2(j) is in module_clinam.f90
! zgval, zgpval and zgppval are computed within sbrtn geval
!   which is in .../tsc/tsc_m/source/geval.f90

   IF ( j > 1 ) THEN
      call geval (xsv2(j),2,zgval,zgpval,zgppval,imag,jmag)
      zbsq = bsqar(j)*1.E-8_R8

      zal0 = zrmid * zetpera * gxmja2(j) / xmja2(j) / zbsq
      zalstar = ftrap(j) * zgval**2 * zrmid * zetpera / zbsq

      zalinc  = 61._R8*(ztemidd/ztimidd)**1.5_R8*(zal0            &  
     &    * ( 1._R8 + 2._R8 * qprof2(j)**2 )                   &  
     &       + 0.46_R8*zalstar/(1._R8-0.54_R8*ftrap(j)) )

      cs2(j) = 0._R8
      cs3(j) = 0._R8
      cs4(j) = 0._R8

      zcoef3nc = -(tpi*qprof2(j))**2 * ztimidd * zalinc

!.....note:  gps is the surface averaged gradient of the toroidal flux squared

      zgps = gja2(j) / vp2(j) *( tpi*qprof2(j))
      d2s(j) = zgps

      zchiinc = -zcoef3nc / (udst*zgps)
!
!....NOTE:  Set electron thermal conductivity to ion value.

      zchienc = -zcoef3nc*acoef(73)/(udst*zgps)/sqrt(2._R8*1835._R8)

      zchiisecj  = zchiinc
      chiinca(j) = zchiinc
      zchiesecj  = zchienc
      chienca(j) = zchienc
!
   ENDIF

!     add contribution from GLF23
!
!      if(itrmode.eq.7) go to 4000

      IF ( itrmod .eq. 14 ) THEN

        ztflux = float(j-1)*dpsi
        IF ( ztflux .lt. ipmax*dpsi ) THEN
          zchiesecj = sqrt(chietema(j)**2 + chienca(j)**2)
          zchiisecj = sqrt(chiitima(j)**2 + chiinca(j)**2)
        ELSE
          zchiesecj = chietema(j)
          zchiisecj = chiitima(j)
        ENDIF

      ELSE

!.....enter here for itrmode=8,9,10,11

        zchiesecj = sqrt(chietema(j)**2 + chienca(j)**2 )
!     &                 +chiecopi(j)**2)
        zchiisecj = sqrt(chiitima(j)**2 + chiinca(j)**2 )
!     &                 +chiicopi(j)**2)

      ENDIF

!     set bounds

      IF( chiitima(j)+dchii_dtip(j)*dti_dpsi(j) .lt. chiinca(j) ) THEN
        chiitima(j)   = zchimin
        dchii_dtep(j) = 0.0
        dchii_dtip(j) = 0.0
      ENDIF

      IF( chietema(j)+dchie_dtep(j)*dte_dpsi(j) .lt. chienca(j) ) THEN
        chietema(j)   = zchimin
        dchie_dtep(j) = 0.0
        dchie_dtip(j) = 0.0
      ENDIF

      zcoef3 = -(zchiisecj+chiitima(j)/zchiisecj*dchii_dtip(j)*dti_dpsi(j))
      zcoef4 = -chiitima(j)/zchiisecj*dchii_dtep(j)*dti_dpsi(j)
      zcoef5 = chiitima(j)/zchiisecj*(                                    &
     &                    dchii_dtip(j)*dti_dpsi(j)*dti_dpsi(j)+        &
     &                    dchii_dtep(j)*dti_dpsi(j)*dte_dpsi(j))        
      zcoef1 = -chietema(j)/zchiesecj*dchie_dtip(j)*dte_dpsi(j)
      zcoef2 = -(zchiesecj+chietema(j)/zchiesecj*dchie_dtep(j)*dte_dpsi(j))
      zcoef6 = chietema(j)/zchiesecj*(                                    &
     &                    dchie_dtep(j)*dte_dpsi(j)*dte_dpsi(j)+        &
     &                    dchie_dtip(j)*dte_dpsi(j)*dti_dpsi(j))        
 
!
      zcoef3 = zcoef3*udst*zgps
      zcoef4 = zcoef4*udst*zgps
      zcoef5 = zcoef5*udst*zgps*(ztemidd/ztemid)
!
      zcoef1 = zcoef1*udst*zgps
      zcoef2 = zcoef2*udst*zgps
      zcoef6 = zcoef6*udst*zgps*(ztemidd/ztemid)

      dsi0(j) = zcoef5*zrmid
      dsi1(j) = -(zcoef3*ztimidd+zcoef4*ztemidd)
      dsi2(j) = zcoef3
      dsi3(j) = -zcoef3 + zcoef4
      dse0(j) = zcoef6 * zrmid
      dse1(j) = -zcoef2*ztemidd - zcoef1*ztimidd
      dse2(j) = zcoef1
      dse3(j) = zcoef2 - zcoef1

      cs1(j) = sqrt(diffnema(j)**2 + diffary(j)**2)*zgps*udst / zrmid

!.....add in convection term only for diffary part  (added 03/28/03)

      dse1(j) = dse1(j) - zface*diffary(j)*zgps*udst / zrmid
      dsi1(j) = dsi1(j) - zfaci*diffary(j)*zgps*udst / zrmid

go to 557
  IF ( 1 == iskipsf ) THEN
      WRITE(*,*) j      &
        , dsi0(j)   &
        , dsi1(j)   &
        , dsi2(j)   &
        , dsi3(j)   &
        , dse0(j)   &
        , dse1(j)   &
        , dse2(j)   &
        , dse3(j)
  ENDIF
557 continue

ENDDO

RETURN
END SUBROUTINE tsc_trcdef_neoclass7
!-----------------------------------------------------------------
SUBROUTINE tsc_trcdef_sawtooth
!
! implement sawtooth model to modify:
!                   etpara
!                   dsej , dsij, csj     for j=0,4
!
! Code copied from .../tsc/tr_tsc/source/trcdef.f90

!-----------------------------------------------------------------

USE clinam            ! TSC common blocks in file module_clinam.f90
USE saprop            ! Allocatable arrays in file module_saprop.f90

USE newplot
USE runaway

IMPLICIT NONE
INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

REAL*8 :: zgps, zemid, zpmid, zsum2

INTEGER :: &
  j          ! Do loop index

!-----------------------------------------------------------------
! 
      call sawtooth
!
!-----------------------------------------------------------------
!.....define energy confinement times in physical units

DO j=1,npsit
  zgps = gja2(j) / vp2(j) * (tpi * qprof2(j) )
  IF ( zgps .eq. 0 ) CYCLE
  chiesec(j) = -(dse3(j)+dse2(j))/(udst*zgps)
  chiisec(j) = -dsi2(j)/(udst*zgps)
ENDDO

      chiesec(1)=chiesec(2)
      chiisec(1)=chiisec(2)
      chiecopi(1) = chiecopi(2)
      chiicopi(1) = chiicopi(2)
      chienca(1)  = chienca(2)
      chiinca(1)  = chiinca(2)

IF ( kcycle .lt. 0 ) RETURN

!-----------------------------------------------------------------
!     note:  <J.B> =  (zgval**2/(tpi**2*vp2(j)))
!    2                *(adi(j)*gxmja(j)*xmja(j)
!    3                - adi(j-1)*gxmja(j-1)*xmja(j-1))*rdpsi

DO  j=1,npsit

!.....define loop voltage due to hyperresistivity term

      vlooph(j+1) = -(avhyp(j+1)-avhyp(j))/dpsi

!.....subtract off bootstrap and current-drive terms

      as0(j) = -( tpi**2 / ( 0.5_R8 * (adi(j)+adi(j+1)) *               &
     & xmja2(j))) * etpara(j) * etafac(j) *                             &  
     & tpi * ( ajavbs(j) + ajavcd(j) + ajavfw(j) + ajavlh2(j) + ajavec(j))

!.....now define cs and ds arrays for convenience - - surface centered

      as(j) =                                                            &  
     &    as0(j)                                                         &  
     &  + as1(j)*(adn(j+1)/vp(j+1) - adn(j)/vp(j))*rdpsi                 &  
     &  + as2(j)*(adp(j+1)/vpg(j+1) - adp(j)/vpg(j))*rdpsi               &  
     &  + as3(j)*(ade(j+1)/vpg(j+1) - ade(j)/vpg(j))*rdpsi               &  
     &  + rdpsi*etpara(j)*etafac(j)*tpi**2/(xmja2(j)*                    &
     &    0.5_R8*(adi(j)+adi(j+1)) )**2                                  &  
     &  *(adi(j+1)*gxmja(j+1)*xmja(j+1)                                  &  
     &                            -adi(j)*gxmja(j)*xmja(j))
!
      cs(j) =                                                            &  
     &    cs0(j)                                                         &  
     &  + cs1(j)*(adn(j+1)/vp(j+1) - adn(j)/vp(j))*rdpsi                 &  
     &  + cs2(j)*(adp(j+1)/vpg(j+1) - adp(j)/vpg(j))*rdpsi               &  
     &  + cs3(j)*(ade(j+1)/vpg(j+1) - ade(j)/vpg(j))*rdpsi               &  
     &  + cs4(j)*(adi(j+1)*gxmja(j+1)*xmja(j+1)                          &  
     &            - adi(j)*gxmja(j)*xmja(j))*rdpsi
!
      dsi(j) =                                                           &  
     &     dsi0(j)                                                       &  
     &   + dsi1(j)*(adn(j+1)/vp(j+1)-adn(j)/vp(j))*rdpsi                 &  
     &   + dsi2(j)*(adp(j+1)/vpg(j+1)-adp(j)/vpg(j))*rdpsi               &  
     &   + dsi3(j)*(ade(j+1)/vpg(j+1)-ade(j)/vpg(j))*rdpsi               &  
     &   + dsi4(j)*(adi(j+1)*gxmja(j+1)*xmja(j+1)                        &  
     &              - adi(j)*gxmja(j)*xmja(j))*rdpsi
!
      dse(j) =                                                           &  
     &     dse0(j)                                                       &  
     &   + dse1(j)*(adn(j+1)/vp(j+1)-adn(j)/vp(j))*rdpsi                 &  
     &   + dse2(j)*(adp(j+1)/vpg(j+1)-adp(j)/vpg(j))*rdpsi               &  
     &   + dse3(j)*(ade(j+1)/vpg(j+1)-ade(j)/vpg(j))*rdpsi               &  
     &   + dse4(j)*(adi(j+1)*gxmja(j+1)*xmja(j+1)                        &  
     &              - adi(j)*gxmja(j)*xmja(j))*rdpsi
!
      zemid = .5_R8*(ade(j+1)/vpg(j+1)+ade(j)/vpg(j))
      zpmid = .5_R8*(adp(j+1)/vpg(j+1)+adp(j)/vpg(j))
      IF ( dse(j) .ne. 0 ) ratioe(j) = cs(j)*zemid/dse(j)
      IF (dsi(j) .ne.0 ) ratioi(j) = cs(j)*(zpmid-zemid)/dsi(j)

ENDDO


!.....define cs and ds at j=1 for each transport model

      cs(1) = 0._R8
      dsi(1) = 0._R8
      dse(1) = 0._R8

!.....adjust value of as(npsit) to be continuous

      as(npsit) = as(npsit-1)

      DO j=npsit+1,npsi
        as(j) = as(npsit)
      ENDDO

!.....define as array at j=1

      as(1) = (2._R8*as(2) - as(3))

!.....second def of energy confinement time

      enerst2 = 0.0_R8

      DO j=2,npsit
        enerst2 = enerst2 + &
          udsp * 1.5_R8 * ( vary(j)-vary(j-1) ) * adp(j) / vpg(j)
      ENDDO

RETURN
END SUBROUTINE tsc_trcdef_sawtooth
!-----------------------------------------------------------------------------
END MODULE tsc_fmcfm_transport_mod
