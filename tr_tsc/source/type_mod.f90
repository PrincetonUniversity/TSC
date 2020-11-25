!## ############################################################
!##
!## File:  type_mod.f90
!##
!## Purpose:  Declare the common derived types
!##                          
!##           
!##
!## $Id: type_mod.f90 24 2007-02-16 00:17:02Z srinath $
!## ############################################################
!> This module containes the fundemental data structure for transport
!! calculations. 
!! Data structures in this module include:
!!    1)MagGeom: common magnetic geometry features for all models
!!    2)Species: These are more important characteristics associated
!!        with a species.  Note, anything other than electrons are
!!        ions.
!!    3)Simple Flux construct: We have transport associated
!!      with the three lowest moments of the Boltzman distribution
!!      equaiton (density, momentum, temperature).  Then we consider
!!      the associated particles, which are ions and electron.  This
!!      "allspecies" derived type. FluxFactors is a derived type with
!!      "allspecies"-type members distinquishing the components
!!      of flux associate with a diffusive process on movement due to 
!!      convective velocity.
!!    4) number of ion species: Internal flag  

        MODULE type_mod
        USE constants_mod  !< Specifies kind for data types.
        IMPLICIT NONE


!=======================================================================
!       General Data Structure on Surface
!
!-----------------------------------------------------------------------       

!> The species derived type contains the primary attributes of a spieces.
!! @note \f$\rho =\frac{\Phi_{total}}{\Phi}\f$ is unitless.
!! \f$arho=\sqrt{\frac{\Phi_{total}}{\pi B_T}}\f$ 
!! and \f$\tilde{\rho}=\sqrt{\frac{\Phi}{\pi B_T}}\f$ in GLF23.
!! Thus our \f$\rho=\frac{\tilde{\rho}}{arho}\f$.  
        TYPE :: Species 
            REAL(kind=r8)   :: charge            !< charge state of isotope, (-1) for electrons
            REAL(kind=r8)   :: nprotons          !< number of protons in nucleus, (0) for electrons
            REAL(kind=r8)   :: amu               !< atomic mass number unit, (1/1827)~0 for electrons
            REAL(kind=r8)   :: density           !< \f$ [m^{-3}] \f$
            REAL(kind=r8)   :: gradDen           !< gradient w.r.t. \f$\rho\f $
            REAL(kind=r8)   :: temperature       !< \f$[eV]\f$
            REAL(kind=r8)   :: gradTemp          !< gradient  w.r.t. \f$\rho\f$,
            REAL(kind=r8)   :: vphi              !< toroidal velocity  \f$[m/s]\f$ 
            REAL(kind=r8)   :: vpara             !< parallel velocity \f$[m/s]\f$ 
            REAL(kind=r8)   :: vperp             !< perpendicular velocity \f$[m/s]\f$
            CHARACTER(80)    :: speciesType       !< electron, thermalion, impurity, beamion, rfminority, fusionion, avghydrogenic, avgimpurity, avgfastIon. convention of name is all lowercase
        END TYPE Species
!> A collection of all species within the simulation involved in the transport
!!calculations.
!! @note All averages are density weighted.
        TYPE :: AllSpecies !#SETGET
            TYPE(Species), ALLOCATABLE     :: ion(:)   !< array of ion attributes
            TYPE(Species)              :: electron          !< electron attributes
            TYPE(Species)              :: AvgHydrogenic     !< averaged hyrdogenic ion attributes
            TYPE(Species)              :: AvgImpurity       !< averaged impurity ion attributes
            TYPE(Species)              :: AvgFastIon        !< averaged fast ion attributes
            REAL(kind=r8)              :: Zeff              !< Z effective for all ions
            INTEGER                    :: numIon            !< number of ion species
            INTEGER                    :: numImpurity       !< number of impurity ion species
            INTEGER                    :: numFastIon        !< number of fast ion species
        END TYPE AllSpecies


!> The derived contains magenetic infromation associated with a specific flux
!! surface and device scalars.       
!   DEVELOPERS: Any thing added here needs to be added to Initialization and
!   Dump
        TYPE :: MagGeom !#SETGET
          REAL(kind=r8)  :: rmin                      !< local minor radius in \f$r_{\min}=\frac{R_{\max}(\rho)-R_{\min}(\rho)}{2}~[m]\f$
          REAL(kind=r8)  :: rminor                    !< minor radius of seperatrix in \f$r_{minor}=\frac{R_{\max}(\rho=1)-R_{\min}(\rho=1)}{2}~[m]\f$
          REAL(kind=r8)  :: Rmaj                      !< local major radius in \f$R_{maj}=\frac{R_{\max}(\rho)+R_{\min}(\rho){2}~[m]\f$
          REAL(kind=r8)  :: Rmajor                    !< geometrical major radius at  magnetic axis \f[m]\f$
          REAL(kind=r8)  :: DRmajDrho                 !< shift
          REAL(kind=r8)  :: rho                       !< square root of normalized torodial flux surface i
          REAL(kind=r8)  :: DrDrho                    !< \f$\frac{dr}{d\rho}\f4, where r is the midplane radius 
          REAL(kind=r8)  :: gradRho                   !< \f$< | \nabla \rho | > [1/m]\f$
          REAL(kind=r8)  :: gradRhoSq                 !< \f$< | \nabla \rho|^2 > [1/m^2]\f$
          REAL(kind=r8)  :: arho                      !< \f$arho=\sqrt{\frac{\Phi_{total}}{\pi B_T}}~[m]\f$
          REAL(kind=r8)  :: kappa                     !< elongation, \f$\kappa\f$
          REAL(kind=r8)  :: DkappaDrho                !< \f$\frac{d\kappa}{d\rho}\f$
          REAL(kind=r8)  :: delta                     !< \f$\delta\f$ triangulation
          REAL(kind=r8)  :: deltaMiller               !< \f$\delta_{miller}\f$ triangulation used in GYRO
          REAL(kind=r8)  :: DdeltaDrho                !< \f$\frac{d\delta}{d\rho}\f$
          REAL(kind=r8)  :: DdeltaMillerDrho          !< \f$\frac{d\delta_{miller}}{d\rho}\f$
          REAL(kind=r8)  :: q                         !< safety factor
          REAL(kind=r8)  :: DlnQDlnRho                !< magnetic shear \f$d \ln q / d \ln \rho\f$
          REAL(kind=r8)  :: Bt                        !< device vaccum field, Bt = constant(a)/R(a)  \f$[T]\f$
          REAL(kind=r8)  :: RBtor                     !< such that RBt(rho)=constant (rho) for a given flux surface, then constant for a flux surface                         
          REAL(kind=r8) :: fluxAvgB                  !< \f$<B>~[T]\f$  
          REAL(kind=r8) :: BSq                       !< \f$<B^2>~[T^2]\f$  
          REAL(kind=r8) :: invBsq                    !< \f$<1/B^2>~[1/T^2]\f$
          REAL(kind=r8) :: invBSqGradRhoSq           !< \f$<\frac{| \nabla \rho|^2}{B^2} > ~[\frac{1}{m^2 T^2}]\f$
          REAL(kind=r8) :: coulombLog                !< \f$ \ln \Lambda \f$ Coulomb logarithm
          REAL(kind=r8) :: volumePrime               !< \f$ \frac{dV}{d \rho} ~ [m^{-3}] \f$
          REAL(kind=r8) :: localInvAspect            !< \f$ \frac{r}{R_0} \f$
          REAL(kind=r8) :: gradRhoSqInvRmajSq         !< \f$ <\frac{| \nabla \rho|^2}{R_{maj}^2}> \f$
          REAL(kind=r8) :: invRmajSq                 !< \f$ <\frac{1}{R_{maj}^2}> \f$

!          CHARACTER(67)  :: getandset='rmin rmaj rho kappa q DlnQDlnRho gradRho gradRhoSq arho rmajor Bt'  !*FD User changable variables.
        END TYPE MagGeom
!!
!!-----------------------------------------------------------------------        
!!    OUTPUT
!! sv 7/15/08:  Many changes to the structure.  No use of simple
!!              language interoperability.  Leaning on Babel for this.
!!           The type we want here are the fluxes because they are used in
!!       the solve.  However, diffusivities and frequencies/rates can be
!!       useful outputs as well so we allow them to be initialized if
!!       available.  
!! sv 3/13/07: Changing pointer to allocatable and will use SIZEOF(array)=0
!!               as the logic choice
!! sv 4/03/07:  Need to change to have 2 derive types (diffs and fluxes)
!! sv 5/06/08:  Changed to used allocatable nested derived types.        
!!-----------------------------------------------------------------------
!
!!=================PRIMARY OUTPUT===================================================
!> The moment derived type refers to the zero-th, first and second moment of the
!! distribution function \f$f_s\f$.  
        TYPE :: moment
            REAL(kind=r8)   :: particle !< Units: \f$1/{m^2 s} \f$
            REAL(kind=r8)   :: thermal  !< Units: \f$J/{m^2 s} \f$
            REAL(kind=r8)   :: momentum !< Units: \f$kg/{m^2 s} \f$ 
        END TYPE moment
!
!> The TransSpecies derived type is intented to be titled "flux", exampled as
!! "flux%electron%particle". 
        TYPE :: TransSpecies  !#SETGET
            TYPE(moment),  ALLOCATABLE     :: ion(:)  !< transport flux associated with ions
            TYPE(moment)               :: electron          !< transport flux associated with electrons
            TYPE(moment)               :: AvgHydrogenic     !< transport flux associated with averaged hydrogenic ions
            TYPE(moment)               :: AvgImpurity       !< transport flux associated with averaged impurity ions
            TYPE(moment)               :: AvgFastIon        !< transport flux associated with averaged fast ions
        END TYPE TransSpecies     
!>If a specific model can compute specific fluxes
!! 0: no flux computed
!! 1: only 1 flux will be returned, this may require an average before flux calculation 
!! 2: the model returns quantities equal to number of species         
        TYPE :: Avail
            INTEGER   :: particle 
            INTEGER   :: thermal  
            INTEGER   :: momentum 
        END TYPE Avail

        TYPE :: CompFluxes
                TYPE(Avail)            :: bulkion
                TYPE(Avail)            :: electron
                TYPE(Avail)            :: impurity
        END TYPE  CompFluxes
!
!
!!================Secondary OUTPUT======================================================
!> Type of flux factor.
        TYPE :: FluxFactors
           REAL(kind=r8) :: diffusivity             !< \f$[m^2/s]\f$
           REAL(kind=r8) :: convectiveVel           !< \f$[m/s]\f$
        END TYPE FluxFactors
!> Type of flux moment.
        TYPE :: TransType
            TYPE(FluxFactors)  :: particle          !< denotes density contribution
            TYPE(FluxFactors)  :: momentum          !< denotes mementum contribution
            TYPE(FluxFactors)  :: thermal           !< denotes energy density contribution
!           CHARACTER(40)       :: get='thermal particle momentum'
        END TYPE TransType
!> Collection of flux factors, as exampled as
!! fluxcoeff%electron%thermal%diffusivity = \f$\chi_e\f$
        TYPE :: TransSpeciesCoeff
            TYPE(TransType), ALLOCATABLE      :: ion(:)    !< transport coefficients associated with ions 
            TYPE(TransType)               :: electron           !< transport coefficients associated with electrons
            TYPE(TransType)               :: AvgHydrogenic      !< transport coefficients associated with averaged hydrogenic ions
            TYPE(TransType)               :: AvgImpurity        !< transport coefficients associated with averaged impurity ions 
            TYPE(TransType)               :: AvgFastIon         !< transport coefficients associated with averaged fast ions
!            CHARACTER(40)      :: get='ion_# electron AvgImpurity AvgFastIon'
        END TYPE TransSpeciesCoeff     
!=================================================================================
 

 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        INTEGER                          :: jcount
        

!---------------------------Type_mod specific routines--------------------- 
        CONTAINS

!--------------------------------------------------------------------------------- 
!           Initialization of Data Structures
!--------------------------------------------------------------------------------- 
!--------------------------------------------------------------------------------- 
!> Used in model specific flag initialization.  Necessary to identify
!! which fluxes are trully calculated from the model.
        SUBROUTINE InitCompFluxes(cflux)
            TYPE(compfluxes)      :: cflux


            cflux%bulkion%particle=0
            cflux%bulkion%thermal=0
            cflux%bulkion%momentum=0
            cflux%electron%particle=0
            cflux%electron%thermal=0
            cflux%electron%momentum=0
            cflux%impurity%particle=0
            cflux%impurity%thermal=0
            cflux%impurity%momentum=0


        END SUBROUTINE InitCompFluxes

!--------------------------------------------------------------------------------- 
!       SUBROUTINE FmInitMagGeom
!-----------------------------------------------------------------------
!> It is necessary to initialize the magentic geometric quantities to zero prior to
!!setting.
        SUBROUTINE FmInitMagGeom(eqMG,error ) !#WRAP:INIT
        TYPE(MagGeom), INTENT(OUT)             :: eqMG
        INTEGER                                :: error
!-----------------------------------------------------------------------    
! mageom
        eqMG%rmin=0_r8
        eqMG%rminor=0_r8
        eqMG%Rmaj=0_r8
        eqMG%DRmajDrho=0_r8
        eqMG%rho=0_r8
        eqMG%DrDrho=0_r8
        eqMG%kappa=0_r8
        eqMG%DkappaDrho=0_r8
        eqMG%delta=0_r8
        eqMG%deltaMiller=0_r8
        eqMG%DdeltaDrho=0_r8
        eqMG%DdeltaMillerDrho=0_r8
        eqMG%q=0_r8
        eqMG%DlnQDlnRho=0_r8
        eqMG%gradRho=0_r8
        eqMG%gradRhoSq=0_r8
        eqMG%arho=0_r8
        eqMG%Bt=0_r8                        ![Tesla]
        eqMG%RMajor=0_r8       
        eqMG%RBtor=0_r8
        eqMG%fluxAvgB=0_r8
        eqMG%BSq=0_r8
        eqMG%invBSq=0_r8
        eqMG%invBSqGradRhoSq=0_r8
        eqMG%coulombLog=0_r8
        eqMG%volumePrime=0_r8
        eqMG%localInvAspect=0_r8
        eqMG%gradRhoSqInvRmajSq=0_r8
        eqMG%InvRmajSq=0_r8
 
        RETURN
        END SUBROUTINE FmInitMagGeom
!--------------------------------------------------------------------------------- 
!--------------------------------------------------------------------------------- 
!       SUBROUTINE InitIon
!-----------------------------------------------------------------------
!> This initializes the ion species' attributes to zero.
!> @note The primary method to use is InitAllSpecies.
        SUBROUTINE InitIon(ion,nspec,error)
        TYPE(species), ALLOCATABLE, INTENT(INOUT)  :: ion(:) 
        INTEGER                       :: nspec !< number of ion species
        INTEGER                       :: icount,error, nspec2
!-----------------------------------------------------------------------  
        error=0
        ALLOCATE(ion(nspec))
        IF(.NOT.ALLOCATED(ion)) THEN
         CALL killFmcfm("InitIon: Must allocate.",error)
        ENDIF

        nspec2=SIZE(ion)
        DO icount=1,nspec2
        
            ion(icount)%charge= 1        
            ion(icount)%nprotons= 1.       
            ion(icount)%amu= 1.00794_r8       
            ion(icount)%temperature= 0_r8       
            ion(icount)%gradTemp= 0_r8      
            ion(icount)%density= 0_r8      
            ion(icount)%gradDen=0_r8       
            ion(icount)%vphi= 0_r8
            ion(icount)%vpara= 0_r8       
            ion(icount)%vperp= 0_r8 
            ion(icount)%speciesType='thermalion'
        END DO
        END SUBROUTINE InitIon
!--------------------------------------------------------------------------------- 
!--------------------------------------------------------------------------------- 
!       SUBROUTINE DelIon
!-----------------------------------------------------------------------
!> This deallocates the ion species
!> @note The primary method to use is DelAllSpecies.
        SUBROUTINE DelIon(ion,error)
        TYPE(species), ALLOCATABLE,INTENT(INOUT)  :: ion(:)
        INTEGER                       :: icount,error, nspec2
!-----------------------------------------------------------------------  
        IF(ALLOCATED(ion)) THEN
         DEALLOCATE(ion)
        ELSE
         CALL killFmcfm("DelIon:Must preallocate.",error)
        ENDIF

        END SUBROUTINE DelIon
!-----------------------------------------------------------------------
!       SUBROUTINE InitElectron
!-----------------------------------------------------------------------
!> This initializes the electron species' attributes to zero.
!> @note The primary method to use is InitAllSpecies.

        SUBROUTINE InitElectron(elec,error)
        TYPE(species), INTENT(OUT)  :: elec
        INTEGER                     :: error
!-----------------------------------------------------------------------  

        elec%charge= -1        
        elec%nprotons= 0_r8      
        elec%amu= 5.4892E-4_r8       
        elec%temperature= -HUGE(1)       
        elec%gradTemp= -HUGE(1)       
        elec%density= -HUGE(1)
        elec%gradDen= -HUGE(1)
        elec%vphi= 0_r8       
        elec%vpara= 0_r8       
        elec%vperp= 0_r8       
        elec%speciesType='electron'

        END SUBROUTINE InitElectron 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> This initializes the average hydrogenic species' attributes to zero.
!> @note The primary method to use is InitAllSpecies.
        SUBROUTINE InitAvgHydrogenic(myspecies,error)
        TYPE(Species), INTENT(INOUT)  :: mySpecies
        INTEGER                       :: error
!-----------------------------------------------------------------------  

        mySpecies%charge= 0_r8        
        mySpecies%nprotons= 0_r8       
        mySpecies%amu= 0_r8       
        mySpecies%temperature= 0_r8       
        mySpecies%gradTemp= 0_r8       
        mySpecies%density= 0_r8
        mySpecies%gradDen= 0_r8
        mySpecies%vphi= 0_r8       
        mySpecies%vpara= 0_r8       
        mySpecies%vperp= 0_r8       
        mySpecies%speciesType='avghydrogenic'

        END SUBROUTINE InitAvgHydrogenic
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> This initializes the average impurity species' attributes to zero.
!> @note The primary method to use is InitAllSpecies.
        SUBROUTINE InitAvgImpurity(myspecies,error)
        TYPE(Species), INTENT(OUT)  :: mySpecies
        INTEGER                       :: error
!-----------------------------------------------------------------------  

        mySpecies%charge= 1        
        mySpecies%nprotons= 0_r8       
        mySpecies%amu= 0_r8       
        mySpecies%temperature= 0_r8       
        mySpecies%gradTemp= 0_r8       
        mySpecies%density= 0_r8
        mySpecies%gradDen= 0_r8
        mySpecies%vphi= 0_r8       
        mySpecies%vpara= 0_r8       
        mySpecies%vperp= 0_r8       
        mySpecies%speciesType='avgimpurity'
        END SUBROUTINE InitAvgImpurity
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> This initializes the average fast ion species' attributes to zero.
!> @note The primary method to use is InitAllSpecies.
        SUBROUTINE InitAvgFastIon(myspecies,error)
        TYPE(Species), INTENT(OUT)  :: mySpecies
        INTEGER                       :: error
!-----------------------------------------------------------------------  

        mySpecies%charge= 0_r8        
        mySpecies%nprotons= 0_r8       
        mySpecies%amu= 0_r8       
        mySpecies%temperature= 0_r8       
        mySpecies%gradTemp= 0_r8       
        mySpecies%density= 0_r8
        mySpecies%gradDen= 0_r8
        mySpecies%vphi= 0_r8       
        mySpecies%vpara= 0_r8       
        mySpecies%vperp= 0_r8       
        mySpecies%speciesType='avgfastion'
        END SUBROUTINE InitAvgFastIon
!!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!       SUBROUTINE FmInitAllSpecies
!-----------------------------------------------------------------------
!> A method that initializes all the species to default values.
!! @note This is the prefered way to initialize all species attributes
        SUBROUTINE FmInitAllSpecies(myspecies,num,error) !#WRAP:INIT
        TYPE(AllSpecies), INTENT(INOUT)       :: myspecies !< The full collection.
        INTEGER, INTENT(IN)                   :: num !< number of intended ion species 
        INTEGER                               :: error
!-----------------------------------------------------------------------
        error=0
        CALL InitIon(myspecies%ion,num,error)
        CALL InitAvgHydrogenic(myspecies%AvgHydrogenic,error)
        CALL InitAvgImpurity(myspecies%AvgImpurity,error)
        CALL InitAvgFastIon(myspecies%AvgFastIon,error)
        CALL InitElectron(myspecies%electron,error)
        myspecies%Zeff=0_r8
        myspecies%numImpurity=0
        myspecies%numFastIon=0
        myspecies%numIon=num


        IF(SIZE(myspecies%ion) /= myspecies%numIon) THEN
        CALL killFmcfm("InitAllSpieces:numSpec.neq.allcoated",error) 
        ENDIF

        END SUBROUTINE FmInitAllSpecies
!-----------------------------------------------------------------------
!> A method that deletes all the species.
!! @note This is the prefered way clean things up.
        SUBROUTINE FmDelAllSpecies(myspecies,error) !#WRAP:DEL
        TYPE(AllSpecies), INTENT(INOUT)       :: myspecies !< The full collection.
        INTEGER                               :: error
        IF(ALLOCATED(mySpecies%ion)) THEN
         CALL DelIon(mySpecies%ion,error)
        ENDIF

        END SUBROUTINE FmDelAllSpecies
!-----------------------------------------------------------------------
!           SUBROUTINE FmInitTransSpecies        
!-----------------------------------------------------------------------
!> Initializes fluxes to default values of a HUGE number.
!! This will allow the user to quickly identify which transport flux
!! is not implemented by the selected transport model.

        SUBROUTINE  FmInitTransSpecies(TranspOut,nspec,error) !#WRAP:INIT
        TYPE(TransSpecies), INTENT(INOUT)         :: TranspOut
        INTEGER                                   :: nspec
        INTEGER                                   :: error
        INTEGER                                   :: nspec2, icount
        DOUBLE PRECISION                        :: thehuge
!-----------------------------------------------------------------------
        error = 0 
        ! already initialized
        IF(ALLOCATED(TranspOut%ion)) THEN
          RETURN
        ENDIF
        ! allocate new arrays
        ALLOCATE(TranspOut%ion(nspec))
        !
        nspec2=SIZE(TranspOut%ion)
        DO icount = 1,nspec2
           TranspOut%ion(icount)%particle = zero
           TranspOut%ion(icount)%thermal  = zero
           TranspOut%ion(icount)%momentum = zero
        END DO
        !
        TranspOut%electron%thermal       = zero
        TranspOut%AvgHydrogenic%thermal  = zero
        TranspOut%AvgImpurity%thermal    = zero
        TranspOut%AvgFastIon%thermal     = zero
        TranspOut%electron%particle      = zero
        TranspOut%AvgHydrogenic%particle = zero
        TranspOut%AvgImpurity%particle   = zero
        TranspOut%AvgFastIon%particle    = zero
        TranspOut%electron%momentum      = zero
        TranspOut%AvgHydrogenic%momentum = zero
        TranspOut%AvgImpurity%momentum   = zero
        TranspOut%AvgFastIon%momentum    = zero
        !  
        END SUBROUTINE FmInitTransSpecies
!-----------------------------------------------------------------------
!> Clean allocated arrays assocaited with fluxes.
        SUBROUTINE  FmDelTransSpecies(TranspOut,error) !#WRAP:DEL
        TYPE(TransSpecies), INTENT(InOut)         :: TranspOut
        INTEGER                                   :: error

        IF(ALLOCATED(TranspOut%ion)) THEN
         DEALLOCATE(TranspOut%ion)
        ELSE
         CALL killFmcfm("DelTransSpecies: Not preallocated",error)
        ENDIF

        END SUBROUTINE  FmDelTransSpecies
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> Initializes transport coefficients to zero.
!! This will allow the user for addition of diffusive and convective parts.
        SUBROUTINE  InitTransSpeciesCoeff(TranspOut,nspec,error)
        TYPE(TransSpeciesCoeff), INTENT(InOut)         :: TranspOut
        INTEGER                                   :: nspec
        INTEGER                                   :: nspec2, icount
        INTEGER                                   :: error
        DOUBLE PRECISION                        :: thehuge
!-----------------------------------------------------------------------
        error=0       
        IF (ALLOCATED(TranspOut%ion)) RETURN
        !
        ALLOCATE(TranspOut%ion(nspec))
        !
        nspec2=SIZE(TranspOut%ion)
        DO icount = 1,nspec2
           TranspOut%ion(icount)%thermal%diffusivity    = 0_r8
           TranspOut%ion(icount)%thermal%convectiveVel  = 0_r8
           TranspOut%ion(icount)%particle%diffusivity = 0_r8
           TranspOut%ion(icount)%particle%convectiveVel = 0_r8
           TranspOut%ion(icount)%momentum%diffusivity = 0_r8
           TranspOut%ion(icount)%momentum%convectiveVel = 0_r8
        END DO
        !
        TranspOut%electron%thermal%diffusivity       = 0_r8
        TranspOut%electron%thermal%convectiveVel       = 0_r8
        TranspOut%AvgHydrogenic%thermal%diffusivity  = 0_r8
        TranspOut%AvgHydrogenic%thermal%convectiveVel  = 0_r8
        TranspOut%AvgImpurity%thermal%diffusivity    = 0_r8
        TranspOut%AvgImpurity%thermal%convectiveVel    = 0_r8
        TranspOut%AvgFastIon%thermal%diffusivity     = 0_r8
        TranspOut%AvgFastIon%thermal%convectiveVel     = 0_r8
        TranspOut%electron%particle%diffusivity      = 0_r8
        TranspOut%electron%particle%convectiveVel      = 0_r8
        TranspOut%AvgHydrogenic%particle%diffusivity = 0_r8
        TranspOut%AvgHydrogenic%particle%convectiveVel = 0_r8
        TranspOut%AvgImpurity%particle%diffusivity   = 0_r8
        TranspOut%AvgImpurity%particle%convectiveVel   = 0_r8
        TranspOut%AvgFastIon%particle%diffusivity    = 0_r8
        TranspOut%AvgFastIon%particle%convectiveVel    = 0_r8
        TranspOut%electron%momentum%diffusivity      = 0_r8
        TranspOut%electron%momentum%convectiveVel      = 0_r8
        TranspOut%AvgHydrogenic%momentum%diffusivity = 0_r8
        TranspOut%AvgHydrogenic%momentum%convectiveVel = 0_r8
        TranspOut%AvgImpurity%momentum%diffusivity   = 0_r8
        TranspOut%AvgImpurity%momentum%convectiveVel   = 0_r8
        TranspOut%AvgFastIon%momentum%diffusivity    = 0_r8
        TranspOut%AvgFastIon%momentum%convectiveVel    = 0_r8
        !  
        END SUBROUTINE InitTransSpeciesCoeff
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> Clean allocated arrays assocaited with fluxes.
        SUBROUTINE  DelTransSpeciesCoeff(TranspOut,error)
        TYPE(TransSpeciesCoeff), INTENT(InOut)         :: TranspOut
        INTEGER                                   :: error

        IF(ALLOCATED(TranspOut%ion)) THEN
         DEALLOCATE(TranspOut%ion)
        ELSE
         CALL killFmcfm("DelTransSpeciesCoeff: must preallocate",error)
        ENDIF

        END SUBROUTINE  DelTransSpeciesCoeff

!-----------------------------------------------------------------------
         ! Calculate averages**************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> This calculates the average hydrodenic attributes from user defined all
!! species. 
!! @note The perfered method for average calculations is CalcAvgAllSpecies.
        SUBROUTINE CalcAvgHydrogenicAs(myions)         
        TYPE(allSpecies)       :: myions
        INTEGER                :: num,nspec        
        REAL(kind=r8)          :: denSum
        REAL(kind=r8)          :: tempCharge,tempNprotons,tempAmu
        REAL(kind=r8)          :: tempGradDen,tempTemp,tempGradTemp
        REAL(kind=r8)          :: tempVphi,tempVpara,tempVperp

        num=0 ; denSum=0_r8 
        !nspec=SIZE(myions%ion)
        nspec=myions%numIon

        IF(myions%AvgHydrogenic%density==0) THEN 
        DO jcount=1,nspec
         IF(myions%ion(jcount)%nprotons==1 .AND. TRIM(myions%ion(jcount)%speciesType)=='thermalion') THEN
            num=num+1
            denSum=denSum + myions%ion(jcount)%density
         ENDIF
        END DO
            IF (num /=0.) THEN
            myions%AvgHydrogenic%density=denSum/num 
            END IF
        ENDIF

     IF(myIons%AvgHydrogenic%density /= 0.) THEN   
!    Density weighted averages.
        tempCharge = 0
        tempNprotons = 0
        tempAmu = 0
        tempGradDen = 0
        tempTemp = 0
        tempGradTemp = 0
        tempVphi = 0
        tempVpara = 0
        tempVperp = 0
        DO jcount=1,nspec
         IF(myions%ion(jcount)%nprotons==1 .AND. TRIM(myions%ion(jcount)%speciesType)=='thermalion') THEN
           
           tempCharge=tempCharge &
     &    + myions%ion(jcount)%charge*myions%ion(jcount)%density

             tempNprotons=tempNprotons &
     &    + myions%ion(jcount)%nprotons*myions%ion(jcount)%density         
        
             tempAmu=tempAmu &
     &     + myions%ion(jcount)%amu * myions%ion(jcount)%density

          tempGradDen=tempGradDen &
     &     + myions%ion(jcount)%gradDen * myions%ion(jcount)%density

          tempTemp=tempTemp &
     &     + myions%ion(jcount)%temperature * myions%ion(jcount)%density

          tempGradTemp=tempGradTemp &
     &     + myions%ion(jcount)%gradTemp * myions%ion(jcount)%density

          tempVphi=tempVphi &
     &     + myions%ion(jcount)%vphi * myions%ion(jcount)%density

          tempVpara=tempVpara &
     &     + myions%ion(jcount)%vpara * myions%ion(jcount)%density

          tempVperp=tempVperp &
     &     + myions%ion(jcount)%vperp * myions%ion(jcount)%density

         ENDIF
        END DO




        IF( myions%AvgHydrogenic%charge ==0_r8) THEN 
        myions%AvgHydrogenic%charge=tempCharge &
     &      /denSum
        ENDIF

        IF( myions%AvgHydrogenic%Nprotons ==0_r8) THEN 
        myions%AvgHydrogenic%Nprotons=tempNprotons &
     &      /denSum
        ENDIF

        IF(myions%AvgHydrogenic%amu ==0_r8) THEN
        myions%AvgHydrogenic%amu  = tempAmu &
     &      /denSum
        ENDIF

        IF(myions%AvgHydrogenic%gradDen ==0_r8) THEN
        myions%AvgHydrogenic%gradDen  = tempGradDen &
     &          /denSum
        ENDIF

        IF(myions%AvgHydrogenic%temperature ==0_r8 ) THEN
        myions%AvgHydrogenic%temperature=tempTemp &
     &          /denSum
        ENDIF

        IF(myions%AvgHydrogenic%gradTemp ==0_r8 ) THEN
        myions%AvgHydrogenic%gradTemp = tempGradTemp &
     &          /denSum
        ENDIF

        IF(myions%AvgHydrogenic%vphi ==0_r8) THEN
        myions%AvgHydrogenic%vphi  = tempVphi &
     &          /denSum
        ENDIF

        IF(myions%AvgHydrogenic%vpara ==0_r8) THEN
        myions%AvgHydrogenic%vpara  = myions%AvgHydrogenic%vpara &
     &          /denSum
        ENDIF

        IF(myions%AvgHydrogenic%vperp ==0_r8) THEN
        myions%AvgHydrogenic%vperp  = tempVperp &
     &          /denSum
        ENDIF

        END IF

        END SUBROUTINE CalcAvgHydrogenicAs
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!> This calculates the average impurity attributes from user defined all
! species. 
! @note The perfered method for average calculations is CalcAvgAllSpecies.
        SUBROUTINE CalcAvgImpurityAs(myions,num)         
        TYPE(allSpecies)       :: myions
        INTEGER ,INTENT(OUT)               :: num      
        INTEGER                :: nspec        
        REAL(kind=r8)          :: denSum
        REAL(kind=r8)          :: tempCharge,tempNprotons,tempAmu
        REAL(kind=r8)          :: tempGradDen,tempTemp,tempGradTemp
        REAL(kind=r8)          :: tempVphi,tempVpara,tempVperp

        num=0 ; denSum=0 
        nspec=myions%numIon
        write(*,*) 
        IF(myions%AvgImpurity%density==0_r8) THEN 
        DO jcount=1,nspec
         IF(TRIM(myions%ion(jcount)%speciesType)=='impurity') THEN
            num=num+1
            denSum=denSum + myions%ion(jcount)%density
         ENDIF
        END DO
            IF (num /=0_r8 ) THEN
            myions%AvgImpurity%density=denSum/num 
            END IF
        ELSE
            denSum=myions%AvgImpurity%density
        ENDIF

        IF ( myIons%AvgImpurity%density /= 0_r8) THEN
!    Density weighted averages.
        tempCharge = 0
        tempNprotons = 0
        tempAmu = 0
        tempGradDen = 0
        tempTemp = 0
        tempGradTemp = 0
        tempVphi = 0
        tempVpara = 0
        tempVperp = 0

        DO jcount=1,nspec
            
         IF(TRIM(myions%ion(jcount)%speciesType)=='impurity') THEN
           
           tempCharge=tempCharge &
     &    + myions%ion(jcount)%charge*myions%ion(jcount)%density

             tempNprotons=tempNprotons &
     &    + myions%ion(jcount)%nprotons*myions%ion(jcount)%density         
        
             tempAmu=tempAmu &
     &     + myions%ion(jcount)%amu * myions%ion(jcount)%density

          tempGradDen=tempGradDen &
     &     + myions%ion(jcount)%gradDen * myions%ion(jcount)%density

          tempTemp=tempTemp &
     &     + myions%ion(jcount)%temperature * myions%ion(jcount)%density

          tempGradTemp=tempGradTemp &
     &     + myions%ion(jcount)%gradTemp * myions%ion(jcount)%density

          tempVphi=tempVphi &
     &     + myions%ion(jcount)%vphi * myions%ion(jcount)%density

          tempVpara=tempVpara &
     &     + myions%ion(jcount)%vpara * myions%ion(jcount)%density

          tempVperp=tempVperp &
     &     + myions%ion(jcount)%vperp * myions%ion(jcount)%density

         ENDIF
        END DO

        IF( myions%AvgImpurity%charge ==0_r8) THEN 
        myions%AvgImpurity%charge=tempCharge &
     &      /denSum
        ENDIF

        IF( myions%AvgImpurity%Nprotons ==0_r8) THEN 
        myions%AvgImpurity%Nprotons=tempNprotons &
     &      /denSum
        ENDIF

        IF(myions%AvgImpurity%amu ==0_r8) THEN
        myions%AvgImpurity%amu  = tempAmu &
     &      /denSum
        ENDIF

        IF(myions%AvgImpurity%gradDen ==0_r8) THEN
        myions%AvgImpurity%gradDen  = tempGradDen &
     &          /denSum
        ENDIF

        IF(myions%AvgImpurity%temperature ==0_r8 ) THEN
        myions%AvgImpurity%temperature=tempTemp &
     &          /denSum
        ENDIF

        IF(myions%AvgImpurity%gradTemp ==0_r8 ) THEN
        myions%AvgImpurity%gradTemp = tempGradTemp &
     &          /denSum
        ENDIF

        IF(myions%AvgImpurity%vphi ==0_r8) THEN
        myions%AvgImpurity%vphi  = tempVphi &
     &          /denSum
        ENDIF

        IF(myions%AvgImpurity%vpara ==0_r8) THEN
        myions%AvgImpurity%vpara  = myions%AvgImpurity%vpara &
     &          /denSum
        ENDIF

        IF(myions%AvgImpurity%vperp ==0_r8) THEN
        myions%AvgImpurity%vperp  = tempVperp &
     &          /denSum
        ENDIF
        END IF

        END SUBROUTINE CalcAvgImpurityAs
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> This calculates the average Fast Ion attributes from user defined all
!! species. 
!! @note The perfered method for average calculations is CalcAvgAllSpecies.
        SUBROUTINE CalcAvgFastIonAs(myions,num)
        TYPE(allSpecies)       :: myions
        INTEGER, INTENT(OUT)   :: num
        INTEGER                ::nspec
        REAL(kind=r8)          :: denSum
        
        

        num=0 ; denSum=0 
        !nspec=SIZE(myions%ion)
        nspec=myions%numIon

 
        DO jcount=1,nspec
         IF  (TRIM(myions%ion(jcount)%speciesType)=='beamion' &
     &   .OR. TRIM(myions%ion(jcount)%speciesType)=='rfminority' &   
     &   .OR. TRIM(myions%ion(jcount)%speciesType)=='fusionion'   ) THEN
            num=num+1
            denSum= denSum + myions%ion(jcount)%density
         ENDIF
        END DO

        IF (num /=0.) THEN 
            myions%AvgFastIon%density=denSum/num 
        END IF

        IF (myIons%AvgFastIon%density /= 0_r8 ) THEN 
!    Density weighted averages.
        DO jcount=1,nspec
         IF  (TRIM(myions%ion(jcount)%speciesType)=='beamion' &
     &   .OR. TRIM(myions%ion(jcount)%speciesType)=='rfminority' &   
     &   .OR. TRIM(myions%ion(jcount)%speciesType)=='fusionion') THEN
           
           myions%AvgFastIon%charge=myions%AvgFastIon%charge  &
     &    + myions%ion(jcount)%charge*myions%ion(jcount)%density
          
           myions%AvgFastIon%Nprotons=myions%AvgFastIon%charge  &
     &    + myions%ion(jcount)%Nprotons*myions%ion(jcount)%density

          myions%AvgFastIon%amu=myions%AvgFastIon%amu &
     &     + myions%ion(jcount)%amu * myions%ion(jcount)%density

          myions%AvgFastIon%gradDen= myions%AvgFastIon%gradDen &
     &     + myions%ion(jcount)%gradDen * myions%ion(jcount)%density

          myions%AvgFastIon%temperature= myions%AvgFastIon%temperature &
     &    + myions%ion(jcount)%temperature * myions%ion(jcount)%density

          myions%AvgFastIon%gradTemp= myions%AvgFastIon%gradTemp &
     &     + myions%ion(jcount)%gradTemp * myions%ion(jcount)%density

          myions%AvgFastIon%vphi= myions%AvgFastIon%vphi &
     &     + myions%ion(jcount)%vphi * myions%ion(jcount)%density

          myions%AvgFastIon%vpara= myions%AvgFastIon%vpara &
     &     + myions%ion(jcount)%vpara * myions%ion(jcount)%density

          myions%AvgFastIon%vperp= myions%AvgFastIon%vperp &
     &     + myions%ion(jcount)%vperp * myions%ion(jcount)%density

         ENDIF
        END DO

        IF(myions%AvgFastIon%charge==0_r8) THEN
        myions%AvgFastIon%charge=myions%AvgFastIon%charge &
     &          /denSum
        END IF
        IF(myions%AvgFastIon%Nprotons==0_r8) THEN
        myions%AvgFastIon%Nprotons=myions%AvgFastIon%Nprotons &
     &          /denSum
        END IF
        IF(myions%AvgFastIon%amu==0_r8) THEN
        myions%AvgFastIon%amu  = myions%AvgFastIon%amu &
     &          /denSum
        END IF
        IF(myions%AvgFastIon%gradDen==0_r8) THEN
        myions%AvgFastIon%gradDen  = myions%AvgFastIon%gradDen &
     &          /denSum
        END IF
        IF(myions%AvgFastIon%temperature==0_r8) THEN
        myions%AvgFastIon%temperature=myions%AvgFastIon%temperature &
     &          /denSum
        END IF
        IF(myions%AvgFastIon%gradTemp==0_r8) THEN
        myions%AvgFastIon%gradTemp = myions%AvgFastIon%gradTemp &
     &          /denSum
        END IF
        IF(myions%AvgFastIon%vphi==0_r8) THEN
        myions%AvgFastIon%vphi  = myions%AvgFastIon%vphi &
     &          /denSum
        END IF
        myions%AvgFastIon%vpara  = myions%AvgFastIon%vpara &
     &          /denSum
        IF(myions%AvgFastIon%vperp ==0_r8) THEN
        myions%AvgFastIon%vperp  = myions%AvgFastIon%vperp &
     &          /denSum
        END IF

        END IF


        END SUBROUTINE CalcAvgFastIonAs
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> This calculates the z-effective for the user defined species.
!! @note The prefered way is to use CalcAvgAllSpecies.
        SUBROUTINE CalcZeff(mySpecies)
        TYPE(allSpecies)       :: mySpecies
        INTEGER                :: num,nspec        
        REAL(kind=r8)          :: denSum, temp_top, temp_bottom

        num=0 ; denSum=0 ; temp_top=0.;temp_bottom=0.
        nspec=mySpecies%numIon


        IF(mySpecies%Zeff==0_r8) THEN
            DO jcount=1,nspec
            temp_top = temp_top+ mySpecies%ion(jcount)%density*mySpecies%ion(jcount)%charge**2
            temp_bottom = temp_bottom +  mySpecies%ion(jcount)%density*mySpecies%ion(jcount)%charge
            ENDDO
          mySpecies%Zeff=temp_top/temp_bottom
        END IF

        END SUBROUTINE CalcZeff
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> Calculates all averages for AllSpecies type and Z effective.
        SUBROUTINE FmCalcAvgAllSpecies(mySpecies,error) !#WRAP:CALCAVG
        TYPE(allSpecies)       :: mySpecies
        INTEGER                :: error
        INTEGER                ::numImp, numFastIon

        IF(.NOT.ALLOCATED(mySpecies%ion)) THEN 
            CALL killFmcfm('FmCalcAvgAllSpecies: Must initialize.',error)
        END IF

        
        CALL CalcAvgImpurityAs(mySpecies,numImp)
        CALL CalcAvgFastIonAs(mySpecies,numFastIon)
        CALL CalcAvgHydrogenicAs(mySpecies)
        CALL CalcZeff(mySpecies)

        mySpecies%numImpurity=numImp
        mySpecies%numFastIon=numFastIon

        END SUBROUTINE FmCalcAvgAllSpecies
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> Caculate the average flux due to Impurity species.
!! @note The prefered way is to use CalcAvgTransSpecies.
        SUBROUTINE CalcAvgImpurityTs(myFlux,mySpecies)
        TYPE(TransSpecies)      :: myFlux
        TYPE(AllSpecies)        :: mySpecies
        INTEGER                 :: jcount
        REAL(kind=r8)           :: tempParticleflux, tempThermalFlux
        REAL(kind=r8)           :: tempMomentumFlux
            
            tempParticleflux=0_r8
            tempThermalflux=0_r8
            tempMomentumflux=0_r8
        DO jcount=1,mySpecies%numIon
         IF(TRIM(mySpecies%ion(jcount)%speciesType)=='impurity' &
     &  .OR. TRIM(mySpecies%ion(jcount)%speciesType)=='avgimpurity') THEN
          tempParticleflux=tempParticleFlux +myFlux%ion(jcount)%particle*mySpecies%ion(jcount)%density
            tempThermalflux=tempThermalFlux + myFlux%ion(jcount)%thermal*mySpecies%ion(jcount)%density
            tempMomentumflux=tempMomentumFlux + myFlux%ion(jcount)%momentum*mySpecies%ion(jcount)%density
         ENDIF
        END DO
     
        

          IF(mySpecies%AvgImpurity%density==0_r8) THEN
          myFlux%AvgImpurity%particle=zero
          myFlux%AvgImpurity%thermal=zero
          myFlux%AvgImpurity%momentum=zero
         ELSE
          myFlux%AvgImpurity%particle=tempParticleflux/mySpecies%AvgImpurity%density
          myFlux%AvgImpurity%thermal=tempthermalflux/mySpecies%AvgImpurity%density
          myFlux%AvgImpurity%momentum=tempMomentumflux/mySpecies%AvgImpurity%density
          END IF 


        END SUBROUTINE CalcAvgImpurityTs
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> Caculate the average flux due to FastIon species.
!! @note The prefered way is to use CalcAvgTransSpecies.
        SUBROUTINE CalcAvgFastIonTs(myFlux,mySpecies)
        TYPE(TransSpecies)      :: myFlux
        TYPE(AllSpecies)        :: mySpecies
        INTEGER                 :: jcount
        REAL(kind=r8)           :: tempParticleflux, tempThermalFlux
        REAL(kind=r8)           :: tempMomentumFlux

            tempParticleflux=0_r8
            tempThermalflux=0_r8
            tempMomentumflux=0_r8
        DO jcount=1,mySpecies%numIon
         IF  (TRIM(mySpecies%ion(jcount)%speciesType)=='beamion' &
     &   .OR. TRIM(mySpecies%ion(jcount)%speciesType)=='rfminority'  &  
     &   .OR. TRIM(mySpecies%ion(jcount)%speciesType)=='fusionion'   &
     &   .OR. TRIM(mySpecies%ion(jcount)%speciesType)=='fastion') THEN
            tempParticleflux=tempParticleFlux + myFlux%ion(jcount)%particle*mySpecies%ion(jcount)%density
            tempThermalflux=tempThermalFlux + myFlux%ion(jcount)%thermal*mySpecies%ion(jcount)%density
            tempMomentumflux=tempMomentumFlux + myFlux%ion(jcount)%momentum*mySpecies%ion(jcount)%density
         ENDIF
        END DO
 
          IF(mySpecies%AvgFastIon%density==0_r8) THEN
          myFlux%AvgFastIon%particle=zero
          myFlux%AvgFastIon%thermal=zero
          myFlux%AvgFastIon%momentum=zero
         ELSE
          myFlux%AvgFastIon%particle=tempParticleflux/mySpecies%AvgFastIon%density
          myFlux%AvgFastIon%thermal=tempthermalflux/mySpecies%AvgFastIon%density
          myFlux%AvgFastIon%momentum=tempMomentumflux/mySpecies%AvgFastIon%density
          END IF
           
        END SUBROUTINE CalcAvgFastIonTs
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> Caculate the average flux due to hydrogenic species
!! @note The prefered way is to use CalcAvgTransSpecies.
        SUBROUTINE CalcAvgHydrogenicTs(myFlux,mySpecies)
        TYPE(TransSpecies)      :: myFlux
        TYPE(AllSpecies)        :: mySpecies
        INTEGER                 :: jcount
        REAL(kind=r8)           :: tempParticleflux, tempThermalFlux
        REAL(kind=r8)           :: tempMomentumFlux

            tempParticleflux=0_r8
            tempThermalflux=0_r8
            tempMomentumflux=0_r8

            DO jcount=1,mySpecies%numIon
         IF(mySpecies%ion(jcount)%nprotons==1 .AND. &
            TRIM(mySpecies%ion(jcount)%speciesType)=='thermalion' &
     &   .OR.TRIM(mySpecies%ion(jcount)%speciesType)=='avghydrogenic') THEN
            tempParticleflux=tempParticleFlux + myFlux%ion(jcount)%particle*mySpecies%ion(jcount)%density
            tempThermalflux=tempThermalFlux + myFlux%ion(jcount)%thermal*mySpecies%ion(jcount)%density
            tempMomentumflux=tempMomentumFlux + myFlux%ion(jcount)%momentum*mySpecies%ion(jcount)%density
         ENDIF
        END DO
        IF(mySpecies%AvgHydrogenic%density==0_r8) THEN
          myFlux%AvgHydrogenic%particle=zero
          myFlux%AvgHydrogenic%thermal=zero
          myFlux%AvgHydrogenic%momentum=zero
         ELSE
         myFlux%AvgHydrogenic%particle=tempParticleflux/mySpecies%AvgHydrogenic%density
          myFlux%AvgHydrogenic%thermal=tempthermalflux/mySpecies%AvgHydrogenic%density
          myFlux%AvgHydrogenic%momentum=tempMomentumflux/mySpecies%AvgHydrogenic%density
        ENDIF

        END SUBROUTINE CalcAvgHydrogenicTs
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> The prefered method to calculates the average flux quantities.
        SUBROUTINE FmCalcAvgTransSpecies(myFlux,mySpecies,error) !#WRAP:CALCAVG
        TYPE(TransSpecies)      :: myFlux
        TYPE(AllSpecies)        :: mySpecies
        INTEGER                 :: error

        IF(.NOT.ALLOCATED(myFlux%ion)) THEN 
            CALL killFmcfm('FmCalcAvgTransSpecies: Must initialize.',error)
        END IF

        CALL CalcAvgImpurityTs(myFlux,mySpecies)
        CALL CalcAvgFastIonTs(myFlux,mySpecies)
        CALL CalcAvgHydrogenicTs(myFlux,mySpecies)




        END SUBROUTINE FmCalcAvgTransSpecies
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!        Dump to standard out or file         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> Write to file or standard output with outunit=6
        SUBROUTINE FmDumpMagGeom(eqMG,outunit,error) !#WRAP:DUMP
        TYPE(MagGeom), INTENT(IN)     :: eqMG
        INTEGER            :: outunit
        INTEGER            :: error
        INTEGER           :: funit, icount
!-----------------------------------------------------------------------        
            funit=outunit

! mageom
        WRITE(funit,FMT=*) "rmin= ",eqMG%rmin
        WRITE(funit,FMT=*) "rminon= ",eqMG%rminor
        WRITE(funit,FMT=*) "Rmaj= ",eqMG%Rmaj
        WRITE(funit,FMT=*) "DRmajDrho= ",eqMG%DRmajDrho
        WRITE(funit,FMT=*) "rho= ",eqMG%rho
        WRITE(funit,FMT=*) "dr/drho= ",eqMG%DrDrho      
        WRITE(funit,FMT=*) "elongation,kappa= ",eqMG%kappa
        WRITE(funit,FMT=*) "Dkappa/Drho= ",eqMG%DkappaDrho
        WRITE(funit,FMT=*) "triangulation, delta= ",eqMG%delta
        WRITE(funit,FMT=*) "triangulation, deltaMiller= ",eqMG%deltaMiller
        WRITE(funit,FMT=*) "Ddelta/Drho= ",eqMG%DdeltaDrho
        WRITE(funit,FMT=*) "DdeltaMiller/Drho= ",eqMG%DdeltaMillerDrho
        WRITE(funit,FMT=*) "q= ",eqMG%q
        WRITE(funit,FMT=*) "DlnQDlnRho= ",eqMG%DlnQDlnRho
        WRITE(funit,FMT=*) "gradRho= ",eqMG%gradRho
        WRITE(funit,FMT=*) "gradRhoSq= ",eqMG%gradRhoSq
        WRITE(funit,FMT=*) "arho= ",eqMG%arho
        WRITE(funit,FMT=*) "Bt= ",eqMG%Bt                        ![Tesla]
        WRITE(funit,FMT=*) "RBtor= ",eqMG%RBtor                        ![Tesla]
        WRITE(funit,FMT=*) "RMajor= ",eqMG%RMajor      
        WRITE(funit,FMT=*) "<B>= ",eqMG%fluxAvgB     
        WRITE(funit,FMT=*) "<BSquared>= ",eqMG%BSq      
        WRITE(funit,FMT=*) "<Inverse BSquared>= ",eqMG%invBSq      
        WRITE(funit,FMT=*) "<Grad Rho Sqrd. / Inverse BSquared>= ",eqMG%invBSq      
        WRITE(funit,FMT=*) "Coulomb Logarithm= ",eqMG%coulombLog      
        WRITE(funit,FMT=*) "Volume Prime= ",eqMG%volumePrime      
        WRITE(funit,FMT=*) "local Inverse Aspect Ratio= ",eqMG%localInvAspect     
        WRITE(funit,FMT=*) "<Grad Rho Sqrd. / Inverse Rmaj Squared>= ",eqMG%gradRhoSqInvRmajSq     
        WRITE(funit,FMT=*) "<1 / Inverse Rmaj Squared>= ",eqMG%invRmajSq     

            RETURN
        END SUBROUTINE FmDumpMagGeom
 
!---------------------------------------------------------------------------------        
!       SUBROUTINE FmDumpTransSpeciesDiff
!       assumed that nspecies is defined        
!---------------------------------------------------------------------------------        
!> Easily write transport diffusivity coefficients to file or standard out.
!!@note The preferred method is DumpTransSpeciesCoeff.
        SUBROUTINE DumpTransSpeciesCoeffDiff(trans,outunit) 
        TYPE(TransSpeciesCoeff),INTENT(IN)                :: trans
        INTEGER              :: outunit
        INTEGER           :: funit, icount,i,nspec
!-----------------------------------------------------------------------
            funit=outunit

        WRITE(UNIT=funit,FMT=*) &
     &     "electron particle diffusivity= ", &
     &     trans%electron%particle%diffusivity
        WRITE(UNIT=funit,FMT=*) &
     &     "electron thermal diffusivity= ", &
     &     trans%electron%thermal%diffusivity
        WRITE(UNIT=funit,FMT=*) &
     &     "electron momentum diffusivity= ", &
     &     trans%electron%momentum%diffusivity



        nspec=SIZE(trans%ion)
          DO i=1,nspec 
          WRITE(UNIT=funit,FMT=*) & 
     &     "ion(",i,") particle diffusivity= ", &
     &      trans%ion(i)%particle%diffusivity
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") thermal diffusivity= ", &
     &      trans%ion(i)%thermal%diffusivity
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") momentum diffusivity= ", &
     &      trans%ion(i)%momentum%diffusivity
          END DO 

          WRITE(UNIT=funit,FMT=*) & 
     &     "AvgHydrogenic particle diffusivity= ", &
     &      trans%AvgHydrogenic%particle%diffusivity
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgHydrogenic thermal diffusivity= ", &
     &      trans%AvgHydrogenic%thermal%diffusivity
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgHydrogenic momentum diffusivity= ", &
     &      trans%AvgHydrogenic%momentum%diffusivity


          WRITE(UNIT=funit,FMT=*) & 
     &     "AvgImpurity particle diffusivity= ", &
     &      trans%AvgImpurity%particle%diffusivity
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgImpurity thermal diffusivity= ", &
     &      trans%AvgImpurity%thermal%diffusivity
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgImpurity momentum diffusivity= ", &
     &      trans%AvgImpurity%momentum%diffusivity


          WRITE(UNIT=funit,FMT=*) & 
     &     "AvgFastIon particle diffusivity= ", &
     &      trans%AvgFastIon%particle%diffusivity
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgFastIon thermal diffusivity= ", &
     &      trans%AvgFastIon%thermal%diffusivity
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgFastIon momentum diffusivity= ", &
     &      trans%AvgFastIon%momentum%diffusivity

      END SUBROUTINE DumpTransSpeciesCoeffDiff
!---------------------------------------------------------------------------------        
!       SUBROUTINE DumpTransSpeciesConVel
!       assumed that nspecies is defined        
!---------------------------------------------------------------------------------        
!> Write to file or standard out transport coefficient for convective velocity
!!term.
!!@note The preferred method is DumpTransSpeciesCoeff.
        SUBROUTINE DumpTransSpeciesCoeffConVel(trans,outunit)
        TYPE(TransSpeciesCoeff),INTENT(IN)            :: trans
        INTEGER             :: outunit
        INTEGER           :: funit, icount,i,nspec
!-----------------------`------------------------------------------------
            funit=outunit

        nspec=SIZE(trans%ion)
        WRITE(UNIT=funit,FMT=*) &
     &     "electron particle convectiveVel= ", &
     &     trans%electron%particle%convectiveVel
        WRITE(UNIT=funit,FMT=*) &
     &     "electron thermal convectiveVel= ", &
     &     trans%electron%thermal%convectiveVel
        WRITE(UNIT=funit,FMT=*) &
     &     "electron momentum convectiveVel= ", &
     &     trans%electron%momentum%convectiveVel

          DO i=1,nspec 
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") particle convectiveVel= ", &
     &      trans%ion(i)%particle%convectiveVel
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") thermal convectiveVel= ", &
     &      trans%ion(i)%thermal%convectiveVel
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") momentum convectiveVel= ", &
     &      trans%ion(i)%momentum%convectiveVel
          END DO 

          WRITE(UNIT=funit,FMT=*) & 
     &     "AvgHydrogenic particle convectiveVel= ", &
     &      trans%AvgHydrogenic%particle%convectiveVel
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgHydrogenic thermal convectiveVel= ", &
     &      trans%AvgHydrogenic%thermal%convectiveVel
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgHydrogenic momentum convectiveVel= ", &
     &      trans%AvgHydrogenic%momentum%convectiveVel


          WRITE(UNIT=funit,FMT=*) & 
     &     "AvgImpurity particle convectiveVel= ", &
     &      trans%AvgImpurity%particle%convectiveVel
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgImpurity thermal convectiveVel= ", &
     &      trans%AvgImpurity%thermal%convectiveVel
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgImpurity momentum convectiveVel= ", &
     &      trans%AvgImpurity%momentum%convectiveVel

          WRITE(UNIT=funit,FMT=*) & 
     &     "AvgFastIon particle convectiveVel= ", &
     &      trans%AvgFastIon%particle%convectiveVel

         WRITE(UNIT=funit,FMT=*) &
     &     "AvgFastIon thermal convectiveVel= ", &
     &      trans%AvgFastIon%thermal%convectiveVel
          WRITE(UNIT=funit,FMT=*) &
     &     "AvgFastIon momentum convectiveVel= ", &
     &      trans%AvgFastIon%momentum%convectiveVel
      
        END SUBROUTINE DumpTransSpeciesCoeffConVel
!!---------------------------------------------------------------------------------        
!!---------------------------------------------------------------------------------        
!> Method for writing transport coefficients to file or standard out.
        SUBROUTINE FmDumpTransSpeciesCoeff(trans,outunit,error) !#WRAP:DUMP
        TYPE(TransSpeciesCoeff)            :: trans
        INTEGER             :: outunit
        INTEGER             :: error
        INTEGER           :: funit, icount,i,nspec

!---------------------------------------------------------------------------------        
         funit=outunit
         
        CALL DumpTransSpeciesCoeffDiff(trans,outunit)
        CALL DumpTransSpeciesCoeffConVel(trans,outunit)

        END SUBROUTINE FmDumpTransSpeciesCoeff
!---------------------------------------------------------------------------------        
!---------------------------------------------------------------------------------        
!       SUBROUTINE DumpTransVarsFlux
!       assumed that nspecies is defined        
!---------------------------------------------------------------------------------        
!> Outputs all fluxes to standard out or a file.
        SUBROUTINE FmDumpTransSpeciesFlux(flux,outunit,error) !#WRAP:DUMP
        TYPE(TransSpecies)                :: flux !< The type of flux
        INTEGER           :: outunit                 !< File outunit number, with out then standard out   
        INTEGER           :: error

        INTEGER           :: funit, icount,i,nspec
!-----------------------------------------------------------------------
            funit=outunit


        WRITE(UNIT=funit,FMT=*) &
     &     "electron particle flux= ", &
     &     flux%electron%particle
        WRITE(UNIT=funit,FMT=*) &
     &     "electron thermal flux= ", &
     &     flux%electron%thermal
        WRITE(UNIT=funit,FMT=*) &
     &     "electron momentum flux= ", &
     &     flux%electron%momentum


        nspec=SIZE(flux%ion)

          DO i=1,nspec 
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") particle flux= ", &
     &      flux%ion(i)%particle
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") thermal flux= ", &
     &      flux%ion(i)%thermal
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") momentum flux= ", &
     &      flux%ion(i)%momentum
          END DO 

          WRITE(UNIT=funit,FMT=*)&
     &   "Average Hydrogenic particle flux= ",&
     &   flux%AvgHydrogenic%particle    
      
        WRITE(UNIT=funit,FMT=*)&
     &   "Average Hydrogenic thermal flux= ",&
     &   flux%AvgHydrogenic%thermal    
      
        WRITE(UNIT=funit,FMT=*)&
     &   "Average Hydrogenic momentum flux= ",&
     &   flux%AvgHydrogenic%momentum    

        WRITE(UNIT=funit,FMT=*)&
     &   "Average impurity particle flux= ",&
     &   flux%AvgImpurity%particle
        WRITE(UNIT=funit,FMT=*)&
     &   "Average impurity thermal flux= ",&
     &   flux%AvgImpurity%thermal    
        WRITE(UNIT=funit,FMT=*)&
     &   "Average impurity momentum flux= ",&
     &   flux%AvgImpurity%momentum    
      
        WRITE(UNIT=funit,FMT=*)&
     &   "Average Fast Ion particle flux= ",&
     &   flux%AvgFastIon%particle    
      
        WRITE(UNIT=funit,FMT=*)&
     &   "Average Fast Ion thermal flux= ",&
     &   flux%AvgFastIon%thermal    
      
        WRITE(UNIT=funit,FMT=*)&
     &   "Average Fast Ion momentum flux= ",&
     &   flux%AvgFastIon%momentum    
      
        END SUBROUTINE FmDumpTransSpeciesFlux
!---------------------------------------------------------------------------------        
!---------------------------------------------------------------------------------        
!> Write ion species attributes to file or standard out.
!!@note The prefered method is DumpAllSpecies.
        SUBROUTINE DumpAsIon(mySpecies,num,outunit)
        TYPE(allSpecies), INTENT(IN)    ::  mySpecies
        INTEGER, INTENT(IN)             ::  num
        INTEGER           :: outunit, error
        INTEGER           :: funit, icount,i,nspec
!---------------------------------------------------------------------------------        
       error=0 

            funit=outunit

           

        IF(num > mySpecies%numIon) THEN
            CALL killFmcfm ("DumpAsIon: Incorrect ion.",error)
        END IF
          
            WRITE(UNIT=funit,FMT=*) &
     &     "ion(",num,") charge= ",  &
     &      mySpecies%ion(num)%charge

            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") Nprotons= ",  &
     &      mySpecies%ion(num)%Nprotons
           
            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") amu= ",  &
     &      mySpecies%ion(num)%amu
                                    
            WRITE(UNIT=funit,FMT=*) & 
     &     "ion(",num,") density= ",  &
     &      mySpecies%ion(num)%density
                                                
            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") gradDen= ", & 
     &      mySpecies%ion(num)%gradDen
                                                            
            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") temperature= ", & 
     &      mySpecies%ion(num)%temperature

            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") gradTemp= ", & 
     &      mySpecies%ion(num)%gradTemp

            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") vphi= ",  &
     &      mySpecies%ion(num)%vphi

            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") vpara= ",  &
     &      mySpecies%ion(num)%vpara
            
            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") vperp= ",  &
     &      mySpecies%ion(num)%vperp

            WRITE(UNIT=funit,FMT=*)  &
     &     "ion(",num,") SpeciesType= ", & 
     &      mySpecies%ion(num)%speciesType
           
        END SUBROUTINE DumpAsIon
!!---------------------------------------------------------------------------------        
!!---------------------------------------------------------------------------------        
!> Write average hydrogenic species attributes to file or standard out.
!!@note The prefered method is DumpAllSpecies.
        SUBROUTINE DumpAsAvgHydrogenic(mySpecies,outunit)
        TYPE(allSpecies), INTENT(IN)    ::  mySpecies
        INTEGER            :: outunit
        INTEGER           :: funit, icount,i,nspec
!---------------------------------------------------------------------------------        

            funit=outunit


          
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic charge= ",  &
     &      mySpecies%AvgHydrogenic%charge

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic Nprotons= ",  &
     &      mySpecies%AvgHydrogenic%Nprotons
           
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic amu= ",  &
     &      mySpecies%AvgHydrogenic%amu
                                    
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic density= ",  &
     &      mySpecies%AvgHydrogenic%density
                                                
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic gradDen= ",  &
     &      mySpecies%AvgHydrogenic%gradDen
                                                            
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic temperature= ",  &
     &      mySpecies%AvgHydrogenic%temperature

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic gradTemp= ",  &
     &      mySpecies%AvgHydrogenic%gradTemp

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic vphi= ",  &
     &      mySpecies%AvgHydrogenic%vphi

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic vpara= ",  &
     &      mySpecies%AvgHydrogenic%vpara
            
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic vperp= ",  &
     &      mySpecies%AvgHydrogenic%vperp

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgHydrogenic SpeciesType= ",  &
     &      mySpecies%AvgHydrogenic%speciesType
           
        END SUBROUTINE DumpAsAvgHydrogenic
!!---------------------------------------------------------------------------------   
!!---------------------------------------------------------------------------------        
!> Write average impurity species attributes to file or standard out.
!!@note The prefered method is DumpAllSpecies.
        SUBROUTINE DumpAsAvgImpurity(mySpecies,outunit)
        TYPE(allSpecies), INTENT(IN)    ::  mySpecies
        INTEGER              :: outunit
        INTEGER           :: funit, icount,i,nspec
!---------------------------------------------------------------------------------        

            funit=outunit


          
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity charge= ",  &
     &      mySpecies%AvgImpurity%charge

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity Nprotons= ",  &
     &      mySpecies%AvgImpurity%Nprotons
           
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity amu= ",  &
     &      mySpecies%AvgImpurity%amu
                                    
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity density= ",  &
     &      mySpecies%AvgImpurity%density
                                                
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity gradDen= ",  &
     &      mySpecies%AvgImpurity%gradDen
                                                            
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity temperature= ", & 
     &      mySpecies%AvgImpurity%temperature

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity gradTemp= ",  &
     &      mySpecies%AvgImpurity%gradTemp

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity vphi= ",  &
     &      mySpecies%AvgImpurity%vphi

            WRITE(UNIT=funit,FMT=*) & 
     &     "AvgImpurity vpara= ",  &
     &      mySpecies%AvgImpurity%vpara
            
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity vperp= ",  &
     &      mySpecies%AvgImpurity%vperp

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgImpurity SpeciesType= ",  &
     &      mySpecies%AvgImpurity%speciesType
           
        END SUBROUTINE DumpAsAvgImpurity
!!---------------------------------------------------------------------------------        
!!---------------------------------------------------------------------------------        
!> Write average fast ion species attributes to file or standard out.
!!@note The prefered method is DumpAllSpecies.
        SUBROUTINE DumpAsAvgFastIon(mySpecies,outunit)
        TYPE(allSpecies), INTENT(IN)    ::  mySpecies
        INTEGER              :: outunit
        INTEGER           :: funit, icount,i,nspec
!---------------------------------------------------------------------------------        

            funit=outunit


          
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon charge= ",  &
     &      mySpecies%AvgFastIon%charge

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon Nprotons= ",  &
     &      mySpecies%AvgFastIon%Nprotons
           
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon amu= ",  &
     &      mySpecies%AvgFastIon%amu
                                    
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon density= ",  &
     &      mySpecies%AvgFastIon%density
                                                
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon gradDen= ",  &
     &      mySpecies%AvgFastIon%gradDen
                                                            
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon temperature= ",  &
     &      mySpecies%AvgFastIon%temperature

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon gradTemp= ",  &
     &      mySpecies%AvgFastIon%gradTemp

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon vphi= ",  &
     &      mySpecies%AvgFastIon%vphi

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon vpara= ",  &
     &      mySpecies%AvgFastIon%vpara
            
            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon vperp= ",  &
     &      mySpecies%AvgFastIon%vperp

            WRITE(UNIT=funit,FMT=*)  &
     &     "AvgFastIon SpeciesType= ",  &
     &      mySpecies%AvgFastIon%speciesType
           
        END SUBROUTINE DumpAsAvgFastIon
!
!!---------------------------------------------------------------------------------        
!!---------------------------------------------------------------------------------        
!> Write z-effective user defined species to file or standard out.
!!@note The prefered method is DumpAllSpecies.
        SUBROUTINE DumpAsZeff(myspecies,outunit)
        TYPE(allSpecies), INTENT(IN)    ::  mySpecies
        INTEGER              :: outunit
        INTEGER           :: funit, icount,i,nspec
!---------------------------------------------------------------------------------        

            funit=outunit
            WRITE(UNIT=funit,FMT=*)  &
     &     "Zeff = ",  &
     &      mySpecies%Zeff       
        END SUBROUTINE DumpAsZeff
!!---------------------------------------------------------------------------------        
!!---------------------------------------------------------------------------------        
!>  Write all species attributes to file or standard out.
        SUBROUTINE FmDumpAllSpecies(myspecies,outunit,error) !#WRAP:DUMP
        TYPE(allSpecies), INTENT(IN)    ::  mySpecies
        INTEGER              :: outunit
        INTEGER              :: error
        INTEGER              :: num, icount

        num=mySpecies%numIon

        DO icount=1,num
        CALL DumpAsIon(mySpecies,icount,outunit)
        END DO
        CALL DumpAsAvgHydrogenic(myspecies,outunit)
        CALL DumpAsAvgImpurity(myspecies,outunit)
        CALL DumpAsAvgFastIon(myspecies,outunit)
        CALL DumpAsZeff(myspecies,outunit)


        END SUBROUTINE FmDumpAllSpecies 
!!---------------------------------------------------------------------------------        
!!       SUBROUTINE killFmcfm
!!---------------------------------------------------------------------------------        
        SUBROUTINE killFmcfm(message,error)    
        CHARACTER(*)     :: message
        INTEGER, INTENT(INOUT)          :: error


        WRITE(*,*) message
        WRITE(*,*) "error =", error
      


        END SUBROUTINE killFmcfm
!        
!!================================================================================
        END MODULE type_mod     
