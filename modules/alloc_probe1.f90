      SUBROUTINE alloc_probe1  
      
      USE PROBE1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_probe1   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_probe1  
