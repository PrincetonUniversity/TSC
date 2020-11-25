!## ############################################################
!## -*-f90-*-
!## File:  glf23_mod.F
!##
!## Purpose:  Initialize: glf23Flags
!##                          glf23 internal variables
!## 
!## $Id: glf23_mod.F 24 2007-02-16 00:17:02Z srinath $
!## 20091213 Bateman changed gfac to
!##   eqMg%gradRho / ( eqMg%drdRho * ( eqMg%gradRhoSq + epsilon ) )
!## 20071020 pankin
!##   - added Fortran directives to allow MPI option
!## 20070310 pankin
!##   - added computation of collisionality, gyro-Bohm factor,
!##     velocities, and rates for ExB flow shear stabilization
!## 20070309 pankin
!##   - descriptions of the GLF23 flags and parameters is added
!## 20070307 pankin
!##   - method to compute diffusivities and fluxes 
!##     glf23CalFlux is added
!## 20070302 pankin
!##   - conversion of FluxSurfQuantity for glf2d call  
!## 20070301 pankin
!##   - glf23 riutines are moved to the module and declared 
!##     as private
!##   - new private variables glf23_flags_extra is added
!##   - subroutines glf23SetFlags is added
!##   - additional type for GLF Flags is created and methods for
!##     its initialization are added
!## ############################################################
c> This module provides the specific data structure for GLF23 flags
c! and the methods to access GLF23 calculations of anomouls transport fluxes.
      MODULE glf23_mod
      USE glf23_data_mod !< trying to capture use statments in doxygen
        USE anomalous_type_mod
        IMPLICIT NONE

!  configure stuff
#ifdef HAVE_CONFIG_H
#include "config.h"
#undef HAVE_CONFIG_H
#endif        
!=======================================================================
!  Module Structure
!               1)Initialize internal variables 
!               2)Initialize glf23 flags
!               3)Gather info and call glf2d                
!=======================================================================
        TYPE :: Glf23FlagsBasic            
          INTEGER :: eigen           !< eigenvalue solver: (0) for cgg, 1 for tomsqz, 2 for zgeev
          INTEGER :: nroot           !< no. roots,(8) for default, 12 for impurity dynamics
          INTEGER :: iglf            !< 0 for original GLF23, (1) for retuned version
          INTEGER, DIMENSION(5) :: tranFlags    !< 1:5 transport flags
          INTEGER :: rotationFlag    !< rotation flag
          INTEGER :: BtFlag          !< switch for effective toroidal field use in rhosda
          INTEGER :: UnitNum         !< unit number for outputs
          INTEGER :: ifElec          !< Electric field switch
          INTEGER :: ifETG           !< ETG transport switch
          REAL(kind=r8)    :: cExB            !< effect of ExB stabilization [0..1]
          REAL(kind=r8)    :: calpha          !< effect of alpha stabilization [0..1]
        END TYPE Glf23FlagsBasic
!
        TYPE :: Glf23FlagsExtra       
          INTEGER          :: ikymax          !< maximum ky index
          REAL(kind=r8)    :: xkymin          !< minimum ky value
          REAL(kind=r8)    :: xkymax          !< maximum ky value
          REAL(kind=r8)    :: xky0            !< related to ky spectrum
          REAL(kind=r8)    :: cnorm           !< ITG normalization
          REAL(kind=r8)    :: adamp           !< radial mode damping exponent
          REAL(kind=r8)    :: park            !< parallel ion motion fit
          REAL(kind=r8)    :: alphaP         !< parallel velocity shear fit
          REAL(kind=r8)    :: alpha_e         
          REAL(kind=r8)    :: cbetae          !< set cbetae=1 for fully electromagnetic
          REAL(kind=r8)    :: cnu             !< factor for the collisionality
          INTEGER, DIMENSION(5) :: iflagin    !< Transport flags internal to the GLF module
          REAL(kind=r8), DIMENSION(30)   :: xparam      !< Normalizers and calibration coefficients internal to the module 
        END TYPE Glf23FlagsExtra
!

        TYPE :: Glf23Flags !#SETGET
            TYPE(Glf23FlagsBasic)       :: basic
            TYPE(Glf23FlagsExtra)       :: extra
            INTEGER                     :: mpiComm
!            INTEGER                     :: mpiCommGroup
            TYPE(CompFluxes)            :: fluxesAvail
        END TYPE Glf23Flags

        TYPE(Glf23Flags), PRIVATE       :: glf23_flags
        INTEGER, PRIVATE :: ifInitialized = 0
        !INTEGER, PRIVATE ::  mpiCommGroup=0, mpiComm=0
        INTEGER, PRIVATE ::  mpiComm=0
        REAL(kind=r8), PRIVATE, PARAMETER ::  epsilon = 1.d-34 
        !
        PRIVATE :: CalcDiffusivity

!
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
! Original GLF23 subroutines and functions
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
        PRIVATE :: 
     &   dscal, zaxpy, zcopy, zdscal, zscal, zswap, ZGEMV, ZTRMV,
     &   ZTRSV, ZGERC, ZGEMM, ZTRMM, glf2d, cgg_glf, cbabk2, cbal,
     &   cdiv, comqr, comqr2, corth, csroot,
     &   R8TOMSQZ, R8CQZHES, R8CQZVAL, R8CQZVEC,
     &   ZGEEV, ZGEBAK, ZGEBAL, ZGEHD2, ZGEHRD, ZHSEQR, ZLACGV, ZLACPY,
     &   ZLAHQR, ZLAHRD, ZLARF, ZLARFB, ZLARFG, ZLARFT, ZLARFX, ZLASCL,
     &   ZLASET, ZLASSQ, ZLATRS, ZTREVC, ZUNG2R, ZUNGHR, ZUNGQR, DLABAD,
     &   DLADIV, DLAMC1, DLAMC2, DLAMC4, DLAMC5, XERBLA, idamax, dzasum,
     &   DZNRM2, izamax, zdotc, zdotu, pythag, ZLADIV, ZLANGE,
     &   ZLANHS, DLAMCH, DLAMC3, IEEECK, ILAENV, dcabs1
        PRIVATE :: DLAPY2, DLAPY3, DLAPY3GF
        PRIVATE :: LSAME
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c

      CONTAINS

!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!       SUBROUTINE FmInitFlagsGlf23 
!>          Set GLF flags to the default values
         SUBROUTINE FmInitFlagsGlf23(this,error) !#WRAP:INIT 
           TYPE(Glf23Flags), INTENT(OUT) :: this
           INTEGER                       :: error
!
           CALL InitCompFluxes(glf23_flags%fluxesAvail) 

           glf23_flags%basic%eigen                 = 0
           glf23_flags%basic%nroot                 = 8 
           glf23_flags%basic%iglf                  = 1
           glf23_flags%basic%tranFlags(1)          = 0
           glf23_flags%basic%tranFlags(2)          = 1
           glf23_flags%basic%tranFlags(3)          = 1
           glf23_flags%basic%tranFlags(4)          = 1
           glf23_flags%basic%tranFlags(5)          = 0
           glf23_flags%basic%rotationFlag          = 0
           glf23_flags%basic%BtFlag                = 0
           glf23_flags%basic%cExB                  = 0
           glf23_flags%basic%calpha                = -1
           glf23_flags%basic%unitnum               = 6
           glf23_flags%basic%ifElec                = 1
           glf23_flags%basic%ifETG                 = 1
!
!
           glf23_flags%extra%iflagin(1) = 0
           glf23_flags%extra%iflagin(2) = 1
           glf23_flags%extra%iflagin(3) = 1
           glf23_flags%extra%iflagin(4) = 0
           glf23_flags%extra%iflagin(5) = 3
!
           glf23_flags%extra%xparam     = 0.
           glf23_flags%extra%xparam(3)  = .7
           glf23_flags%extra%xparam(7)  = 1.
           glf23_flags%extra%xparam(9)  = 1.
!ap
           glf23_flags%extra%xparam(10) = 12.
           glf23_flags%extra%xparam(13) = 0.2
           glf23_flags%extra%xparam(14) = 1.
           glf23_flags%extra%xparam(15) =-0.1
           glf23_flags%extra%xparam(17) = 0.1
           glf23_flags%extra%xparam(23) = 1.
!
           glf23_flags%extra%ikymax     = 10
           glf23_flags%extra%xkymin     = .02
           glf23_flags%extra%xkymax     = .5
           glf23_flags%extra%xky0       = .3
!
           glf23_flags%extra%cnorm      = 100.
           glf23_flags%extra%adamp      = .50
           glf23_flags%extra%park       = 0.7
           glf23_flags%extra%alphaP    = 0.50
!
           glf23_flags%extra%cbetae     = 1.e-6
           glf23_flags%extra%cnu        = 1.

           glf23_flags%fluxesAvail%bulkion%particle=1
           glf23_flags%fluxesAvail%bulkion%thermal=1
           glf23_flags%fluxesAvail%electron%thermal=1
           glf23_flags%fluxesAvail%impurity%particle=1

!
           this = glf23_flags
                                           
         END SUBROUTINE FmInitFlagsGlf23
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!       SUBROUTINE Glf23SetFlags
!>          Update GLF flags initialized in InitGLfFlags with the
!!          values by a user
         SUBROUTINE SetFlagsGlf23(this)
!
           TYPE(Glf23Flags), INTENT(IN) :: this
!                    
           REAL(kind=r8) zpi
!
           zpi = atan2 (0.D0, -1.D0)
!
           ifInitialized = 1
! 
           glf23_flags = this 
!
           eigen_gf     = this%basic%eigen
           nroot_gf     = this%basic%nroot
!
           IF (glf23_flags%basic%iglf.eq.1) THEN
             glf23_flags%extra%iflagin(5) = 5
             glf23_flags%extra%xparam(10) = 12.
             glf23_flags%extra%xparam(13) = 0.15_r8
             glf23_flags%extra%xparam(16) = 0.15_r8
             glf23_flags%extra%xparam(17) = 0.25_r8
             glf23_flags%extra%xparam(19) = 1.
             glf23_flags%extra%cnorm      = 50.
             glf23_flags%extra%adamp      = 0.7_r8
             glf23_flags%extra%park       = 0.8_r8
             glf23_flags%extra%alphaP    = 0.35_r8
           ENDIF
!
           IF (this%basic%ifETG.eq.0) glf23_flags%extra%xparam(10)=0.
!
           IF (this%basic%ifElec.eq.-1) THEN
             xnu_gf=1000.
             glf23_flags%extra%xparam(15)=-1.
           ENDIF
!
!                    
           iflagin_gf(1) = glf23_flags%extra%iflagin(1)
           iflagin_gf(2) = glf23_flags%extra%iflagin(2)
           iflagin_gf(3) = glf23_flags%extra%iflagin(3)
           iflagin_gf(4) = glf23_flags%extra%iflagin(4)
           iflagin_gf(5) = glf23_flags%extra%iflagin(5)
!
           xparam_gf     = glf23_flags%extra%xparam
!
           ikymax_gf     = glf23_flags%extra%ikymax
           xkymin_gf     = glf23_flags%extra%xkymin
           xkymax_gf     = glf23_flags%extra%xkymax
           xky0_gf       = glf23_flags%extra%xky0
!     
           cnorm_gf      = glf23_flags%extra%cnorm
           adamp_gf      = glf23_flags%extra%adamp
           park_gf       = glf23_flags%extra%park
           alpha_p_gf    = glf23_flags%extra%alphaP
!----

!           mpiCommGroup = glf23_flags%mpiCommGroup
           mpiComm  = glf23_flags%mpiComm
           !WRITE(*,*) "mpiComm = ", mpiComm
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!        Initialize other internal variables        
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
! 
           rms_theta_gf  = zpi/3.D0
           ghat_gf       = 1.
           gchat_gf      = 1.
           xwell_gf      = 0.
!
           alpha_star_gf = 0.
           alpha_mode_gf = 0.
           xkdamp_gf     = 0.
 
           gamma_e_gf    = 0.
!
         END SUBROUTINE SetFlagsGlf23

!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!       Call GLF2d to set basic flags and to compute diffusivities
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!#ifdef HAVE_MPI
!       SUBROUTINE CalcAnomTranGlf23(eqMg,surfVar,myspecies,fluxCoeff
!     &               ,ad,flags, mp_gcomm, mp_ggroup,ierr)
!#else
          SUBROUTINE CalcAnomTranGlf23(eqMg,surfVar,myspecies,fluxCoeff
     &              ,ad,flags,ierr)
!#endif
           !
           TYPE(MagGeom), INTENT(IN)             :: eqMg
           TYPE(MagGeom)                         :: eqMGInternal
           TYPE(AnomSurfVars), INTENT(INOUT)     :: surfVar
           TYPE(AllSpecies)                      :: myspecies
           TYPE(TransSpeciesCoeff)               :: fluxCoeff
           TYPE(Glf23Flags), INTENT(IN)          :: flags
           TYPE(AnomTransDetails)                :: ad 
           INTEGER(kind=ispec), INTENT(OUT)      :: ierr
!#ifdef HAVE_MPI
!           INTEGER(kind=ispec), INTENT(IN)      :: mp_gcomm
!           INTEGER(kind=ispec), INTENT(IN)      :: mp_ggroup
!#endif
           !
           ierr = 0
           !
           CALL SetFlagsGlf23(flags)
           !  Convert from FMCFM convention to GLF23
           eqMgInternal = eqMg
           eqMgInternal%rho  =   eqMg%rho !*eqMg%arho
           eqMgInternal%gradRho = eqMg%gradRho * eqMg%arho
           eqMgInternal%drdRho  = eqMg%drdRho /  eqMg%arho
           eqMgInternal%gradRhoSq   = eqMg%gradRhoSq * (eqMg%arho**2)
 
          CALL CalcDiffusivity(eqMGInternal,surfVar,myspecies,fluxCoeff,
     &                              ad,ierr)
!           PRINT *,  'ion-th = ', fluxCoeff%ion(1)%thermal%diffusivity
!           PRINT *,  'ele-th = ', fluxCoeff%electron%thermal%diffusivity
!           PRINT *,  'ion-pt = ', fluxCoeff%ion(1)%particle%diffusivity
!           PRINT *,  'imp-pt = ', 
!      &               fluxCoeff%AvgImpurity%particle%diffusivity
          IF (ierr>0) THEN
            ierr = 100 + ierr
            RETURN
          ENDIF        
!#endif
         END SUBROUTINE CalcAnomTranGlf23
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!       Calc conductive flux terms. 
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!>GLF23 units
         SUBROUTINE CalcFluxGlf23(eqMg, myspecies, myFlux, sV,
     &                              glfFlags, ierr)
           TYPE(MagGeom), INTENT(IN)             :: eqMg
           TYPE(AnomSurfVars), INTENT(INOUT)     :: sV
           TYPE(AllSpecies)                      :: myspecies
           TYPE(TransSpecies)                    :: myFlux
           TYPE(AnomTransDetails)                :: ad
           TYPE(Glf23Flags), INTENT(IN)          :: glfFlags
           INTEGER(kind=ispec), INTENT(OUT)      :: ierr
           INTEGER(kind=ispec)                   ::icount,nspecies
           TYPE(TransSpeciesCoeff)                :: fluxCoeff
!#ifdef HAVE_MPI           
!           INTEGER(kind=ispec), INTENT(IN)          :: mp_gcomm
!           INTEGER(kind=ispec), INTENT(IN)       :: mp_ggroup
!#endif

!           ALLOCATE(fluxCoeff%ion(SIZE(myFlux%ion)))
!#ifdef HAVE_MPI           
!           CALL CalcAnomTranGlf23(eqMG,sV,myspecies,fluxCoeff,ad,
!     &           glfFlags,mp_gcomm, mp_ggroup, ierr)
!#else     
          ierr = 0
          !
          CALL InitTransSpeciesCoeff(fluxCoeff,SIZE(myFlux%ion),ierr)
          IF (ierr > 0) THEN
            ierr = ierr+100
            RETURN
          ENDIF
          CALL CalcAnomTranGlf23(eqMg,sV,myspecies,fluxCoeff,ad,
     &                 glfFlags,ierr)
          IF (ierr > 0) THEN
            ierr = ierr+200
            RETURN
          ENDIF

!#endif     
           !
          nspecies=SIZE(myspecies%ion)

          myFlux%ion(1)%thermal=
     &            -fluxCoeff%ion(1)%thermal%diffusivity
     &            *myspecies%ion(1)%gradTemp

          myFlux%electron%thermal=
     &            -fluxCoeff%electron%thermal%diffusivity
     &            *myspecies%electron%gradTemp

          myFlux%ion(1)%particle= 
     &            -fluxCoeff%ion(1)%particle%diffusivity
     &            *myspecies%ion(1)%gradDen

          myFlux%AvgImpurity%particle=
     &            -fluxCoeff%AvgImpurity%particle%diffusivity
     &            *myspecies%AvgImpurity%gradDen

          myFlux%AvgImpurity%thermal= 0.

          CALL DelTransSpeciesCoeff(fluxCoeff,ierr)
          IF (ierr > 0) THEN
            ierr = ierr+300
            RETURN
          ENDIF

         END SUBROUTINE CalcFluxGlf23
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!       This is with metric terms included with unit conversion
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!> Caclulates *Only* the diffusive part of the flux.
         SUBROUTINE FmCalcFluxGlf23(eqMg, myspecies, myFlux, sV,
     &                              glfFlags, ierr) !#WRAP:CALCFLUX
          TYPE(MagGeom), INTENT(IN)             :: eqMg
          TYPE(AllSpecies)                      :: myspecies
          TYPE(TransSpecies)                    :: myFlux
          TYPE(AnomSurfVars), INTENT(INOUT)     :: sV
          TYPE(Glf23Flags), INTENT(IN)          :: glfFlags
          INTEGER(kind=ispec), INTENT(OUT)      :: ierr
          TYPE(AnomTransDetails)                :: ad
          INTEGER(kind=ispec)                   :: icount,nspecies
          TYPE(TransSpeciesCoeff)               :: fluxCoeff
          !
          ierr = 0
          !
          CALL InitTransSpeciesCoeff(fluxCoeff,SIZE(myFlux%ion),ierr)
          IF (ierr > 0) THEN
            ierr = ierr+100
            RETURN
          ENDIF
          CALL CalcAnomTranGlf23(eqMg, sV, myspecies,
     &                       fluxCoeff, ad, glfFlags, ierr)
          IF (ierr > 0) THEN
            ierr = ierr+200
            RETURN
          ENDIF
          !
          nspecies=SIZE(myspecies%ion)
          ! 
          myFlux%ion(1)%particle=(-1.)*(eqMg%gradRhoSq
     &            *fluxCoeff%ion(1)%particle%diffusivity 
     &            *myspecies%ion(1)%gradDen 
     &            +eqMg%gradRho*fluxCoeff%ion(1)%particle%convectiveVel
     &            *myspecies%ion(1)%Density)
      

          myFlux%ion(1)%thermal=(-1.)*(eqMg%gradRhoSq 
     &           *fluxCoeff%ion(1)%thermal%diffusivity 
     &           *myspecies%ion(1)%gradTemp*ev2Joule
     &           *myspecies%ion(1)%Density 
     &           +eqMg%gradRho*fluxCoeff%ion(1)%thermal%convectiveVel 
     &           *myspecies%ion(1)%Density 
     &           *myspecies%ion(1)%Temperature*ev2Joule)
     &           + (5./2.)*myspecies%ion(1)%temperature*ev2Joule
     &           *myFlux%ion(1)%particle

          myFlux%electron%thermal=(-1.)*(eqMg%gradRhoSq
     &           *fluxCoeff%electron%thermal%diffusivity 
     &           *myspecies%electron%gradTemp*ev2Joule
     &           *myspecies%electron%Density
     &           +eqMg%gradRho*fluxCoeff%electron%thermal%convectiveVel
     &           *myspecies%electron%Density
     &           *myspecies%electron%Temperature*ev2Joule)
     &           +(5./2.)*myspecies%electron%temperature*ev2Joule
     &           *myFlux%electron%particle

          myFlux%AvgImpurity%particle=(-1.)*(eqMg%gradRhoSq
     &       *fluxCoeff%AvgImpurity%particle%diffusivity 
     &       *myspecies%AvgImpurity%gradDen      
     &       +eqMg%gradRho*fluxCoeff%AvgImpurity%particle%convectiveVel
     &       *myspecies%AvgImpurity%Density)

          CALL DelTransSpeciesCoeff(fluxCoeff,ierr)
          IF (ierr > 0) THEN
            ierr = ierr+300
            RETURN
          ENDIF
          CALL FmCalcAvgTransSpecies(myFlux,mySpecies,ierr)
          IF (ierr > 0) THEN
            ierr = ierr+400
            RETURN
          ENDIF
!
        END SUBROUTINE FmCalcFluxGlf23


!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!       Call GLF2d to compute diffusivities
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!#ifdef HAVE_MPI
!         SUBROUTINE CalcDiffusivity(eqMG,sV,myspecies,fluxCoeff,ad, 
!     &                              mp_gcomm, mp_ggroup, ierr)
!#else
         SUBROUTINE CalcDiffusivity(eqMg,sV,myspecies,fluxCoeff,ad,ierr)
!#endif
           !
           TYPE(MagGeom), INTENT(IN)             :: eqMg
           TYPE(AnomSurfVars), INTENT(INOUT)     :: sV
           TYPE(AllSpecies), INTENT(IN)          :: myspecies
           TYPE(AnomTransDetails), INTENT(OUT)   :: ad ! anomalous diffusivity 
           TYPE(TransSpeciesCoeff), INTENT(INOUT):: fluxCoeff
           INTEGER(kind=ispec), INTENT(OUT)      :: ierr
           INTEGER(kind=ispec)                   :: nspecies, icount 
           INTEGER(kind=ispec)                   :: bad_count
!#ifdef HAVE_MPI
!           INTEGER(kind=ispec), INTENT(IN) :: mp_gcomm
!           INTEGER(kind=ispec), INTENT(IN) :: mp_ggroup
!#endif
!
!...local variables

           REAL(kind=r8) :: vnewk3x
           REAL(kind=r8) :: nem, nim, alpha_m, betae_m
           REAL(kind=r8) :: drhodrrrho
           REAL(kind=r8) :: cgyrobohm_m
           REAL(kind=r8) :: gfac
           INTEGER       :: im, jm
           REAL(kind=r8), PARAMETER :: znunit = 1E19_r8 !   Because GLF23 expects inputs of this order                                                  
           REAL(kind=r8), PARAMETER :: tunit = 1E-3_r8 !   Because GLF23 expects inputs of this order                                                  
! Note:  Temperature input unit is [keV], but FMCFM is [eV] unit.
!
! Compute variables that dependes on FluxQuantities
           
           REAL(kind=r8) :: zdensz,zavezimp_numerator

           nspecies=SIZE(myspecies%ion)
!
           ierr = 0
           IF ( (eqMg%gradRho*eqMg%gradRhoSq) .eq. 0 ) THEN
              ! gradRho or/and gradRhoSq are no set
              ierr = 1
              RETURN
           ENDIF
!
           rmin_gf       = eqMg%rmin/eqMg%arho
           rmaj_gf       = eqMg%rmaj/eqMg%arho
           elong_gf      = eqMg%kappa
           q_gf          = eqMg%q
!
           nem   = myspecies%electron%density/znunit
           nim   = myspecies%AvgHydrogenic%density/znunit
!            
           dil_gf        = 1_r8 - nim/nem

           zimp_gf=mySpecies%AvgImpurity%charge    
             
           amassimp_gf=mySpecies%AvgImpurity%amu
            
           amassgas_gf  = myspecies%AvgHydrogenic%amu
!
           drhodrrrho = eqMg%rmin / eqMg%drdRho /
     &                  (eqMg%rho+epsilon) / eqMg%arho
!
           shat_gf       = eqMg%DlnQDlnRho*drhodrrrho
!           
           taui_gf = myspecies%AvgHydrogenic%temperature / 
     &               myspecies%electron%temperature
!
           rlti_gf   = -(myspecies%AvgHydrogenic%gradTemp/
     &                myspecies%AvgHydrogenic%temperature)
     &                  / eqMg%drdRho
           rlte_gf   = -(myspecies%electron%gradTemp/
     &                  myspecies%electron%temperature)
     &                  / eqMg%drdRho
           rlne_gf   = -(myspecies%electron%gradDen/
     &                  myspecies%electron%density)
     &                  / eqMg%drdRho
           rlni_gf   = -(myspecies%AvgHydrogenic%gradDen/
     &                   myspecies%AvgHydrogenic%density)
     &                  / eqMg%drdRho
            
            rlnimp_gf= -(mySpecies%AvgImpurity%gradDen/
     &                 mySpecies%AvgImpurity%density)
     &                 / eqMg%drdRho

           apwt_gf   = 1.
           aiwt_gf   = 0.
!... flow shear stabilization
           alpha_e_gf = glf23_flags%basic%cExB
!
!
! ... plasma collisionality
!
           vnewk3x = 0.117_r8*nem
     &       /SQRT(((myspecies%electron%temperature*tunit)**3)
     &      *(myspecies%AvgHydrogenic%temperature*tunit))
     &              *eqMg%arho*SQRT(amassgas_gf/2.)
           xnu_gf  = glf23_flags%extra%cnu*
     &              vnewk3x*SQRT(taui_gf/2.)*myspecies%Zeff
!
! ... plasma beta
!
           betae_m  = 4.e-3_r8*nem
     &      *(myspecies%electron%temperature*tunit)/eqMg%Bt**2
           betae_gf = max(glf23_flags%extra%cbetae*betae_m, 1.e-6_r8)
!
!
! ... alpha stabilization
!
          cgyrobohm_m = 1.e-4*
     &         9.79e5_r8 * SQRT(myspecies%electron%temperature) / 
     &                     ( eqMg%arho*100. )
     &         * ( 1.02D2 * SQRT(myspecies%electron%temperature) 
     &                     / eqMg%Bt / 1.e4 )**2*
     &         SQRT(amassgas_gf)
!
           IF (glf23_flags%basic%calpha.lt.0.) THEN 
             alpha_m = -1./eqMg%drdRho * 
     &             q_gf**2 * eqMg%rmaj * betae_m * 
     &             ( taui_gf*nim/nem*
     &             ( myspecies%AvgHydrogenic%gradDen/
     &                 myspecies%AvgHydrogenic%density
     &              + myspecies%AvgHydrogenic%gradTemp/
     &                 myspecies%AvgHydrogenic%temperature )
     &              + myspecies%electron%gradDen/
     &                 myspecies%electron%density
     &              + myspecies%electron%gradTemp/
     &                 myspecies%electron%temperature )
             alpha_gf  = -glf23_flags%basic%calpha * alpha_m
           ELSE
             alpha_gf  = glf23_flags%basic%calpha * sV%alphaMhd
         ENDIF

! ... ExB stabilization
!
         IF ( (glf23_flags%basic%tranFlags(4).eq.-1).OR.
     &        (glf23_flags%basic%rotationFlag .gt. 0) ) THEN
           gamma_e_gf = sV%wExb
           gamma_p_gf = sV%wPara
           alpha_e_gf = 1
         ENDIF
         !
         IF (sV%wExbDia>0.) gamma_star_gf=sV%wExbDia
!
! Call glf23 module
!
#ifdef HAVE_MPI
          CALL glf2d(glf23_flags%basic%iglf, mpiComm) 
#else
          CALL glf2d(glf23_flags%basic%iglf) 
#endif
!
! GLF23 output
!
!..... Create new structures for glf23 output
!
          CALL InitAnomTransDetails(ad,0_ispec,4_ispec,
     &                               glf23_flags%basic%nroot,20_ispec)
!  
!..... Some geometrical factor
!
          gfac = 
     &  eqMg%gradRho / ( eqMg%drdRho * ( eqMg%gradRhoSq + epsilon ) )
!
!..... Total effective diffusivities
!
!
          fluxCoeff%ion(1)%thermal%diffusivity     = 
     &                          gfac * chii_gf * cgyrobohm_m
          fluxCoeff%electron%thermal%diffusivity   = 
     &                          gfac * chie_gf * cgyrobohm_m
          fluxCoeff%ion(1)%particle%diffusivity    = 
     &                          gfac * diff_gf * cgyrobohm_m

          fluxCoeff%AvgImpurity%particle%diffusivity=
     &                          gfac * diff_im_gf * cgyrobohm_m



          fluxCoeff%AvgImpurity%thermal%diffusivity  = 0. ! is not computed in GLF
!
!..... Toroidal momentum diffusivities
!
          ad%tmomentum%phi    = gfac * eta_phi_gf * cgyrobohm_m
          ad%tmomentum%par    = gfac * eta_par_gf * cgyrobohm_m
          ad%tmomentum%perp   = gfac * eta_per_gf * cgyrobohm_m

          DO im = 1, 4
            ad%freq(im)%gamma = gamma_gf(im)
            ad%freq(im)%omega = freq_gf(im)
            DO jm = 1, 20
              ad%spectrum(im,jm)%gamma = gamma_k_gf(im,jm)
              ad%spectrum(im,jm)%omega = freq_k_gf(im,jm)
            ENDDO
          ENDDO
!
        END SUBROUTINE CalcDiffusivity


          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        SUBROUTINE FmDumpFlagsGlf23(inflags,outunit,error) !#WRAP:DUMP
        IMPLICIT NONE
        TYPE(GLF23Flags) :: inflags
        INTEGER      :: outunit
        INTEGER      :: error
        INTEGER      :: funit
!-----------------------------------------------------------------------
        
          funit=outunit

        WRITE(UNIT=funit,FMT=*) inflags

        END SUBROUTINE FmDumpFlagsGlf23

!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!       Subroutines moved from original NTCC GLF23 module        
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
#include "glf2d.F"
#include "r8tomsqz.F"
#include "zgeev.f"
#include "blas_zgeev.f"


!================================================================================
       END MODULE glf23_mod

