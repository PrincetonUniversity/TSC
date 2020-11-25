      SUBROUTINE alloc_arrays  
      
      USE ARRAYS
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( PF(18), STAT=istat)
        PF = 0 
        ALLOCATE ( PFE(18), STAT=istat)
        PFE = 0 
        ALLOCATE ( FCOM(18), STAT=istat)
        FCOM = 0 
        ALLOCATE ( vchopper(18), STAT=istat)
        vchopper = 0 
        ALLOCATE ( FCOMSV(18), STAT=istat)
        FCOMSV = 0                                                         
            if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_arrays   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_arrays  
