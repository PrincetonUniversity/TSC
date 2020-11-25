!## ############################################################
!##
!## File:  flow_shear_mod.f90
!##
!## Purpose: Module to compute shearing rates 
!## 
!##
!## $Id: flow_shear_mod.f90 24 2007-07-04 09:39:36Z pankin $
!## ########################################################### 
!> This module is used to calculate flow shearing rates as found
!! in the initial driver of the NTCC GLF23.
MODULE flow_shear
  USE spec_kind_mod
  USE glf23_mod
  USE anomalous_type_mod
  !
  IMPLICIT NONE
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c

  CONTAINS

!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!       SUBROUTINE calc_flow_shear 
!>          Calculate ExB flow shear 
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
   SUBROUTINE calcFlowShearGlf23(gf, mg, aS, anomSv, jmaxm, ierr)
     ! 
     TYPE :: shearingRates
      REAL(kind=r8) :: egamma     !< E x B shearing rate (flow shear)
      REAL(kind=r8) :: vstarp     !< Diamagnetic part of E x B drift velocity shear
      REAL(kind=r8) :: gamma_p    !< Parallel velocity shearing rate
     END TYPE shearingRates
     !
     TYPE(Glf23Flags), INTENT(IN)                     :: gf !< flags need to be set
     TYPE(MagGeom), DIMENSION(:), INTENT(IN)          :: mg
     TYPE(allSpecies), DIMENSION(:), INTENT(INOUT)    :: aS
     TYPE(anomSurfVars), DIMENSION(:), INTENT(INOUT)  :: anomSv
     INTEGER, INTENT(IN)       :: jmaxm               ! array dimensions
     INTEGER, INTENT(OUT)      :: ierr
     ! 
     ! internal variables
     ! 
     INTEGER       :: jm
     REAL(kind=r8), ALLOCATABLE :: drhodrrrho(:)
     TYPE(shearingRates), ALLOCATABLE :: sr(:)
     REAL(kind=r8), PARAMETER   :: epsilon = 1.d-34 
     REAL(kind=r8), ALLOCATABLE :: ve(:)
     REAL(kind=r8), ALLOCATABLE :: vpar(:)
     REAL(kind=r8), ALLOCATABLE :: vper(:)
     REAL(kind=r8), ALLOCATABLE :: csda(:)
     REAL(kind=r8), ALLOCATABLE :: rhosda(:)
     REAL(kind=r8)              :: bteff
     REAL(kind=r8)              :: fc
     REAL(kind=r8)              :: alpha_neo
     REAL(kind=r8), ALLOCATABLE :: zppi(:), zppi_neo(:) ! renormalized pressure gradients
     REAL(kind=r8), ALLOCATABLE :: ztau(:)           ! temperature ratio
     !
     ierr = 0
     !
     IF ( (SIZE(mg)*SIZE(aS) < jmaxm**2) .OR. (jmaxm < 2 ) ) THEN
       ierr = 1
       RETURN
     ENDIF
     !
     IF ( ANY(mg(1:jmaxm)%drdRho.eq.0) ) THEN
       print *, jmaxm, mg(1:jmaxm)%drdRho
       ierr = 2
       RETURN
     ENDIF
     !
     ! allocate internal arrays
     !
     ALLOCATE(drhodrrrho(1:jmaxm),ve(1:jmaxm),vpar(1:jmaxm),vper(1:jmaxm))
     ALLOCATE(csda(1:jmaxm),rhosda(1:jmaxm))
     ALLOCATE(sr(1:jmaxm))
     ALLOCATE(ztau(1:jmaxm), zppi(1:jmaxm), zppi_neo(1:jmaxm))
     !
     DO jm=1,jmaxm
       !
       ! Geometric factor
       !
       drhodrrrho(jm) = mg(jm)%rmin / ((mg(jm)%rho+epsilon) * mg(jm)%drdRho )
       ! use of aRho is conversion from FMCFM to GLF23.
       !
       ! Local rate unit for flow shear
       !
       !csda(jm) = 9.75e5_r8*SQRT(sv(jm)%Te)/(mg(jm)%arho*100.)/SQRT(sv(jm)%amassGas)
       ! Temperature in csda is in [eV]
       csda(jm) = 9.79e3_r8 * SQRT(aS(jm)%electron%temperature) / &
                 mg(jm)%arho / SQRT(aS(jm)%avgHydrogenic%amu) 
       ! assuming working gas is Hydrogenic
       ! Temperature in GLF23 and FMCFM in [Kev] 

       !
       ! Local rho_star
       ! Use effective B-field if bt_flag > 0
       !
       IF (gf%basic%BtFlag .GT. 0) THEN
          bteff = mg(jm)%Bt * mg(jm)%rho * mg(jm)%arho**2 / mg(jm)%rmin / mg(jm)%drdRho
          !no need for the extra arho, because they would cancel
          rhosda(jm) = ( 1.02e-4_r8 * SQRT(aS(jm)%electron%temperature) / bteff ) &
                        * SQRT(aS(jm)%avgHydrogenic%amu) / mg(jm)%arho
       ELSE
          rhosda(jm) = ( 1.02e-4_r8 * SQRT(aS(jm)%electron%temperature) / mg(jm)%Bt ) &
                        * SQRT(aS(jm)%avgHydrogenic%amu) / mg(jm)%arho
       ENDIF
       !
       ! Calculate collisionless limit formulas
       !
       fc = 1.-1.46_r8*SQRT(mg(jm)%rmin/mg(jm)%rmaj)+ &
            0.46*(mg(jm)%rmin/mg(jm)%rmaj)**1.5
       alpha_neo = 1.-0.8839_r8*fc/(0.3477_r8+0.4058_r8*fc)
       ! 
       ! normalized presure gradient
       !
       zppi(jm) = -rhosda(jm) * mg(jm)%arho / mg(jm)%drdRho * & 
                 ( aS(jm)%AvgHydrogenic%gradDen/aS(jm)%avgHydrogenic%density + &
                   aS(jm)%AvgHydrogenic%gradTemp/aS(jm)%avgHydrogenic%temperature)
       !
       zppi_neo(jm) = -rhosda(jm) * mg(jm)%arho / mg(jm)%drdRho * &
                      ( aS(jm)%AvgHydrogenic%gradDen/aS(jm)%avgHydrogenic%density + &
                        alpha_neo*aS(jm)%AvgHydrogenic%gradTemp/aS(jm)%avgHydrogenic%temperature )
       !
       ! temperature ratio
       !
       ztau(jm) = aS(jm)%avgHydrogenic%temperature/aS(jm)%electron%temperature
       !
       ! Compute E x B drift velocity
       !
       ! ...Main component
       !
       !sV (8/25/08) at this point I'm going to assume only hydrogenic species
       ! are to be involved with the flow shear rate calculations.
        
       ve(jm) = ztau(jm)*csda(jm)*mg(jm)%arho* &
                zppi_neo(jm)-mg(jm)%rmin/mg(jm)%q*anomSv(jm)%velAng
       !
       ! ...Parallel velocity
       !
       vpar(jm) = mg(jm)%rmajor*anomSv(jm)%velAng - &
                 ztau(jm)*csda(jm) * mg(jm)%arho**2 * &
                 rhosda(jm) / mg(jm)%drdRho * &
                 (alpha_neo-1.) * &
                 aS(jm)%avgHydrogenic%gradTemp/aS(jm)%avgHydrogenic%temperature* &
                 mg(jm)%rmin/mg(jm)%rmajor/mg(jm)%q
       !
       ! Various optional dependencies on parallel, perpendicular, and toroidal
       !     velocities
       !
       IF ((gf%basic%tranFlags(4).EQ.0) .AND. (gf%basic%tranFlags(5).EQ.0)) THEN
         aS(jm)%avgHydrogenic%vphi  = mg(jm)%rmajor*anomSv(jm)%velAng
         aS(jm)%avgHydrogenic%vpara  = vpar(jm)

         aS(jm)%avgHydrogenic%vperp = ve(jm) + &
           ztau(jm)*csda(jm)*mg(jm)%arho*zppi(jm)
       ENDIF
       IF ((ABS(gf%basic%tranFlags(4)).EQ.1).AND.(gf%basic%tranFlags(5).EQ.1)) THEN
         vpar(jm)=aS(jm)%avgHydrogenic%vpara
       ENDIF
!
       IF ((ABS(gf%basic%tranFlags(4)).EQ.1).AND.(gf%basic%tranFlags(5).EQ.0)) THEN
         !
         ve(jm)=ztau(jm)*csda(jm)*mg(jm)%arho*zppi_neo(jm) - &
                mg(jm)%rmin/mg(jm)%q*aS(jm)%avgHydrogenic%vphi/mg(jm)%rmajor
         !
         as(jm)%avgHydrogenic%vperp=ve(jm)-(aS(jm)%avgHydrogenic%temperature/aS(jm)%electron%temperature)* &
            csda(jm)*mg(jm)%arho*zppi(jm)

         vpar(jm)=aS(jm)%avgHydrogenic%vphi - &
                 ztau(jm)*csda(jm) * mg(jm)%arho**2 * &
                 rhosda(jm) / mg(jm)%drdRho * &
                 (alpha_neo-1.) * &
                 aS(jm)%avgHydrogenic%gradTemp/aS(jm)%avgHydrogenic%temperature* &
                 mg(jm)%rmin/mg(jm)%rmajor/mg(jm)%q

       ENDIF
       !
       IF (gf%basic%tranFlags(5).EQ.1) &
          ve(jm) = aS(jm)%avgHydrogenic%vperp + ztau(jm)*csda(jm)*mg(jm)%arho*zppi(jm)
       !
     ENDDO
     !
     ! Computation of shearing rates
     !
     DO jm=1,jmaxm-1
       !
       ! Compute E x B shearing rate (flow shear)
       !
       sr(jm)%egamma=drhodrrrho(jm) * &
                     ( mg(jm)%rho + mg(jm+1)%rho ) / & 
                     ( mg(jm)%q + mg(jm+1)%q ) * &
                     ( ve(jm+1) * mg(jm+1)%q/mg(jm+1)%rmin - ve(jm)*mg(jm)%q/mg(jm)%rmin ) / &
                     ( mg(jm+1)%rho - mg(jm)%rho + epsilon) / csda(jm)
       !
       ! Diamagnetic part of E x B drift velocity shear
       !
       sr(jm)%vstarp=-( mg(jm+1)%rmin + mg(jm)%rmin ) / 2. * &
                      ( ztau(jm+1) * csda(jm+1) * zppi(jm+1) / mg(jm+1)%rmin &
                       - ztau(jm) * csda(jm) * zppi(jm) / mg(jm)%rmin ) / &
                      ( mg(jm+1)%rho - mg(jm)%rho + epsilon) / csda(jm)
       ! 
       ! Parallel velocity shearing rate
       !
       sr(jm)%gamma_p=-1 / mg(jm)%drdRho * & 
                      ( vpar(jm+1) - vpar(jm) ) / &
                      ( mg(jm+1)%rho - mg(jm)%rho + epsilon) / csda(jm)
       !
     ENDDO
     ! linear extrapolation to the last grid point
     IF (jmaxm > 2 ) THEN
       sr(jmaxm)%egamma = sr(jmaxm-1)%egamma+(sr(jmaxm-1)%egamma-sr(jmaxm-2)%egamma)/ &
                          (mg(jmaxm-1)%rho*mg(jmaxm-1)%arho-mg(jmaxm-2)%rho*mg(jmaxm-2)%arho)* &
                          (mg(jmaxm)%rho*mg(jmaxm)%arho-mg(jmaxm-1)%rho*mg(jmaxm-1)%arho)
       sr(jmaxm)%vstarp = sr(jmaxm-1)%vstarp+(sr(jmaxm-1)%vstarp-sr(jmaxm-2)%vstarp)/ &
                          (mg(jmaxm-1)%rho*mg(jmaxm-1)%arho-mg(jmaxm-2)%rho*mg(jmaxm-2)%arho)* &
                          (mg(jmaxm)%rho*mg(jmaxm)%arho-mg(jmaxm-1)%rho*mg(jmaxm-1)%arho)
       sr(jmaxm)%gamma_p = sr(jmaxm-1)%gamma_p+(sr(jmaxm-1)%gamma_p-sr(jmaxm-2)%gamma_p)/ &
                          (mg(jmaxm-1)%rho*mg(jmaxm-1)%arho-mg(jmaxm-2)%rho*mg(jmaxm-2)%arho)* &
                          (mg(jmaxm)%rho*mg(jmaxm)%arho-mg(jmaxm-1)%rho*mg(jmaxm-1)%arho)
     ELSE
       !
       ! E x B shearing rate (flow shear)
       sr(jmaxm)%egamma = sr(jmaxm-1)%egamma
       !
       ! Diamagnetic part of E x B drift velocity shear
       sr(jmaxm)%vstarp = sr(jmaxm-1)%vstarp
       ! 
       ! Parallel velocity shearing rate
       sr(jmaxm)%gamma_p = sr(jmaxm-1)%gamma_p
     ENDIF
     ! 
     ! Update GLF derived types with flow shear rate information 
     !
     DO jm=1, jmaxm
       anomSv(jm)%wExB = sr(jm)%egamma
       anomSv(jm)%wExBdia = sr(jm)%vstarp
       anomSv(jm)%wpara = sr(jm)%gamma_p
     ENDDO
     DEALLOCATE(drhodrrrho,ve,vpar,vper,csda,rhosda,sr)
     DEALLOCATE(ztau, zppi, zppi_neo)
   END SUBROUTINE calcFlowShearGlf23
END MODULE flow_shear
