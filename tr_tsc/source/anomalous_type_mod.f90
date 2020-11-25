!## File:  anomalous_type_mod.F90
!##
!## Purpose:  Declare the common derived types
!##                          
!##           
!##
!## $Id: anomalous_type_mod.F90 24 2007-02-16 00:17:02Z srinath $
!## ############################################################
!> Attributes associated with anomalous models and may not be associated with
!! other transort models.

      MODULE anomalous_type_mod
        USE spec_kind_mod
        USE type_mod
        IMPLICIT NONE
!> Scalars associated with anomalous transport for a specific flux surface.
        TYPE :: AnomSurfVars !#SETGET
           REAL(kind=r8) :: velAng                    !< angular velocity  
           REAL(kind=r8) :: wExb                     !< Experimental E x B shearing rate \f$[1/s]\f$
           REAL(kind=r8) :: wExbDia                 !< Diamagnetic part of E x B shearing rate \f$[1/s]\f$
           REAL(kind=r8) :: wPara                   !< Parallel velocity shearing rate \f$[1/s]\f$
           REAL(kind=r8) :: alphaMhd                  !< alpha_MHD \f$-q^2R\(d\beta /dr\)\f$
        END TYPE AnomSurfVars

             
!> Diffusivitities from flows          
        TYPE :: MomementumDiffusivities !#GET
           REAL(kind=r8) :: phi        ! toroidal momentum diffusivity
           REAL(kind=r8) :: par        ! parallel component of toroidal momentum diffusivity
           REAL(kind=r8) :: perp       ! perpendicular component of toroidal momentum diffusivity
        END TYPE MomementumDiffusivities

!> Growth rate for instabilities in the anomalous         
        TYPE :: Frequencies
           REAL(kind=r8) :: omega
           REAL(kind=r8) :: gamma
        END TYPE Frequencies
!> Rates associated with flows associated with the anomalous model
        TYPE :: Rates
           REAL(kind=r8) :: exch
           REAL(kind=r8) :: ExbGamma
           REAL(kind=r8) :: DexbGamma
           REAL(kind=r8) :: velParGamma
        END TYPE Rates

        INTEGER         :: nionspecies=2  !*FD typically bulk and fast ions 
!---------------------------------------------------------------------------------
!   Common Output for Anomalous Models
!---------------------------------------------------------------------------------
!> Typical output from anomalous model codes.
        TYPE :: AnomTransDetails 
           REAL(kind=r8), ALLOCATABLE         :: diffMatrix(:,:)  ! full diffusivity matrix
           TYPE(TransSpecies), ALLOCATABLE    :: driftmodes(:)   ! contribution of different drift modes
           TYPE(Frequencies), ALLOCATABLE     :: freq(:)         ! frequencies and growth rates
           TYPE(Frequencies), ALLOCATABLE     :: spectrum(:,:)   ! spectrum (GLF23 only)
           TYPE(Rates)                        :: rate           ! growth rates of unstable mode
           TYPE(MomementumDiffusivities)      :: tmomentum      ! toroidal modementum diffusivities  
        END TYPE AnomTransDetails

!==============================================================================
        CONTAINS
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------        
!> Sets anomalous surface variables to default values.

            SUBROUTINE FmInitAnomSurfVars(sV,error) !#WRAP:INIT
            TYPE(AnomSurfVars)          :: sV
            INTEGER                     :: error
!------------------------------------------------------------------------------        
    
           sv%velAng=0.                    
           sv%wExb=0.                     
           sv%wExbDia=0.                 
           sv%wPara =0.                  
           sv%alphaMhd   =0.               
                        
            END SUBROUTINE FmInitAnomSurfVars

!------------------------------------------------------------------------------        
!------------------------------------------------------------------------------        
           SUBROUTINE FmDumpAnomSurfVars(sV,outunit,error) !#WRAP:DUMP
           TYPE(AnomSurfVars)           ::sV
           INTEGER                      :: outunit
           INTEGER                      :: error
           INTEGER                      :: funit

           funit=outunit 


                     
            WRITE(UNIT=funit,FMT=*)  &
     &     "AnomSurfVar velAng= ",  &
     &     sV%velAng

           WRITE(UNIT=funit,FMT=*)  &
     &     "AnomSurfVar wExb= ",  &
     &     sV%wExb

       WRITE(UNIT=funit,FMT=*)  &
     &     "AnomSurfVar wExbDia= ",  &
     &     sV%wExbDia

       WRITE(UNIT=funit,FMT=*)  &
     &     "AnomSurfVar wPara= ",  &
     &     sV%wPara


       WRITE(UNIT=funit,FMT=*)  &
     &     "AnomSurfVar alphaMhd= ",  &
     &     sV%alphaMhd



            END SUBROUTINE FmDumpAnomSurfVars
!------------------------------------------------------------------------------        
!       SUBROUTINE InitAnomTransDetails        


        SUBROUTINE InitAnomTransDetails(AT,nmodes,nfreq,nroot, &
     &          nspectrum)
        TYPE(AnomTransDetails), INTENT(OUT) :: AT
        INTEGER, INTENT(IN):: nmodes, nfreq, nroot, nspectrum
        INTEGER :: icount,lcount
!---------------------------------------------------
!        CALL GetNumIonSpecies(nionspecies)
        IF (nmodes>0) THEN
         ALLOCATE(AT%driftmodes(nmodes))
         DO icount=1,nmodes
           ALLOCATE(AT%driftmodes(icount)%ion(2))
           DO lcount=1,nionspecies
           AT%driftmodes(icount)%ion(lcount)%thermal=0.
           AT%driftmodes(icount)%ion(lcount)%particle=0.
           AT%driftmodes(icount)%ion(lcount)%momentum=0.
           ENDDO
         ENDDO
          AT%driftmodes%electron%thermal=0.
          AT%driftmodes%electron%particle=0.
          AT%driftmodes%electron%momentum=0.
        ENDIF
!
        IF (nfreq>0)  THEN
          ALLOCATE(AT%freq(nfreq))
          AT%freq%omega=0.
          AT%freq%gamma=0.
        ENDIF
!
        IF (nspectrum>0) THEN
          ALLOCATE(AT%spectrum(nfreq,nspectrum))
          AT%spectrum%omega = 0.
          AT%spectrum%gamma = 0.
        ENDIF
!
        AT%rate%exch=0.
        AT%rate%ExbGamma=0.
        AT%rate%DexbGamma=0.
        AT%rate%velParGamma=0.
!
!       Total diffusivity matrix
!
        IF (nroot>0) THEN
          ALLOCATE(AT%diffMatrix(nroot,nroot))
          AT%diffMatrix = 0.
        ENDIF
!
  
!
! toroidal momentum diffusivity
!
        AT%tmomentum%phi         = 0.
        AT%tmomentum%par         = 0.
        AT%tmomentum%perp        = 0.
! 
!
       
        END SUBROUTINE InitAnomTransDetails       
!--------------------------------------------------------------------------------- 
!       SUBROUTINE DelTransIndepth
!--------------------------------------------------------------------------------- 
        SUBROUTINE DelAnomTransDetails(AT)
        TYPE(AnomTransDetails) :: AT
!---------------------------------------------------

        IF (ALLOCATED(AT%driftmodes))      DEALLOCATE(AT%driftmodes)
        IF (ALLOCATED(AT%freq))            DEALLOCATE(AT%freq)
        IF (ALLOCATED(AT%diffMatrix))      DEALLOCATE(AT%diffMatrix)
        IF (ALLOCATED(AT%spectrum))        DEALLOCATE(AT%spectrum)

        
        END SUBROUTINE DelAnomTransDetails
!--------------------------------------------------------------------------------      
!--------------------------------------------------------------------------------
        SUBROUTINE FmDumpAnomTransDetails(AT,outunit,error) !#WRAP:DUMP
        TYPE(AnomTransDetails)  :: AT
        INTEGER      :: outunit
        INTEGER      :: error
        INTEGER      :: funit
        INTEGER      :: idim(2) , i ,j, icount
        CHARACTER*25 :: sfmt,t1,t2
!--------------------------------------------------------------------------------
!        IF(PRESENT(outunit)) THEN
                funit=outunit
!        ELSE
!                funit = 6
!        END IF

        IF (ALLOCATED(AT%driftmodes)) THEN
          WRITE(UNIT=funit,FMT=*) "--------------------------"
         DO icount=1,SIZE(AT%driftmodes)
          WRITE(UNIT=funit,FMT=*)  &
     &     "electron particle diffusivity for mode ",icount," = ", &
     &      AT%driftmodes(icount)%electron%particle 
          DO i=1,nionspecies
            WRITE(UNIT=funit,FMT=*) &
     &      "ion(",i,") particle diffusivity for mode ",icount," = ", &
     &      AT%driftmodes(icount)%ion(i)%particle
          END DO 
          WRITE(UNIT=funit,FMT=*) &
     &     "electron thermal diffusivity for mode ",icount," = ", &
     &    AT%driftmodes(icount)%electron%thermal
          DO i=1,nionspecies
          WRITE(UNIT=funit,FMT=*) &
     &     "ion(",i,") thermal diffusivity for mode ",icount,"  = ", &
     &     AT%driftmodes(icount)%ion(i)%thermal
          END DO
         END DO
        ELSE
         WRITE(UNIT=funit,FMT=*) "--------------------------"
         WRITE(UNIT=funit,FMT=*)  &
     &   "No driftmodes were allocated."
       END IF


        WRITE(UNIT=funit,FMT=*) &
     &            "eta Phi= ", AT%tmomentum%phi
        WRITE(UNIT=funit,FMT=*) &
     &            "eta parallel= ", AT%tmomentum%par 
        WRITE(UNIT=funit,FMT=*) &
     &            "eta perpendicular= ", AT%tmomentum%perp

        IF (ALLOCATED(AT%diffMatrix)) THEN
          idim=SHAPE(AT%diffMatrix)
          WRITE(UNIT=t1, FMT='(I2)') idim(1)
          WRITE(UNIT=t2, FMT='(I2)') idim(2)
          WRITE(UNIT=funit,FMT='(1A)') 'Full diff.matrix:'
          DO i =1,idim(1)
                  WRITE(UNIT=funit,FMT= '('//TRIM(t2)//'(3x,es10.3))') &
     &                       AT%diffMatrix(i,1:idim(2)) 
          END DO
       ENDIF
!
!       WRITE(*,*) 'Convective velocity', AT%convectiveVel
!
        IF (ALLOCATED(AT%freq)) THEN
          idim=SIZE(AT%freq)
          WRITE(UNIT=t1, FMT='(I2)') idim(1)
          WRITE(UNIT=funit,FMT=*) "--------------------------"
          WRITE(UNIT=funit,FMT='(1A,'//TRIM(t1)//'(3x,es10.3))') &
     &       "Growth rates of most unstable modes= ", AT%freq(:)%gamma
          WRITE(UNIT=funit,FMT='(1A,'//TRIM(t1)//'(3x,es10.3))') &
     &       "Frequencies of most unstable modes= ", AT%freq(:)%omega
        ENDIF
!
        IF (ALLOCATED(AT%spectrum)) THEN
          idim=SHAPE(AT%spectrum)
          WRITE(UNIT=t1, FMT='(I2)') idim(1)
          WRITE(UNIT=t2, FMT='(I2)') idim(2)
          WRITE(UNIT=funit,FMT=*)"--------------------------"
          WRITE(UNIT=funit,FMT=*)"Growth rates:"
          DO i =1,idim(2)
               WRITE(UNIT=funit,FMT= '('//TRIM(t1)//'(3x,es10.3))') &
     &                AT%spectrum(1:idim(1),i:i)%gamma
          ENDDO
          WRITE(UNIT=funit,FMT=*)"Frequences:"
          DO i =1,idim(2)
               WRITE(UNIT=funit,FMT= '('//TRIM(t1)//'(3x,es10.3))') &
     &                AT%spectrum(1:idim(1),i:i)%omega
          ENDDO
        ENDIF
!
        WRITE(UNIT=funit,FMT=*) "--------------------------"
        WRITE(UNIT=funit,FMT=*) &
     &  "turbulent electron to ion ENERGY exchange in MW/mFMT=*FMT=*3= " &
     &  ,AT%rate%exch
        WRITE(UNIT=funit,FMT=*) "ExB shear rate = ", AT%rate%ExbGamma
        WRITE(UNIT=funit,FMT=*) "ExB delayed shear rate = ", &
     &  AT%rate%DexbGamma
        WRITE(UNIT=funit,FMT=*) "parallel velcoity shear rate= ", &
     &  AT%rate%velParGamma

        END SUBROUTINE FmDumpAnomTransDetails



!=================================================================================



      END MODULE anomalous_type_mod
