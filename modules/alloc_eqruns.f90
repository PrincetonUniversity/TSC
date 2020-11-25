      SUBROUTINE alloc_eqruns  
      
      USE EQRUNS
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_eqruns   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_eqruns  
