      SUBROUTINE dealloc_arrays  
      
      USE ARRAYS
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( PF, STAT=istat)
        DEALLOCATE ( PFE, STAT=istat)
        DEALLOCATE ( FCOM, STAT=istat)
        DEALLOCATE ( vchopper, STAT=istat)
        DEALLOCATE ( FCOMSV, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_arrays   ' 
      end if
            
      return
      END SUBROUTINE dealloc_arrays  
