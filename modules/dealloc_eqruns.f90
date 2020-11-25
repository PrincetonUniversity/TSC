      SUBROUTINE dealloc_eqruns  
      
      USE EQRUNS
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_eqruns   ' 
      end if
            
      return
      END SUBROUTINE dealloc_eqruns  
