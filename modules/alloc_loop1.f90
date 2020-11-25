      SUBROUTINE alloc_loop1  
      
      USE LOOP1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_loop1   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_loop1  
