      SUBROUTINE alloc_cb  
      
      USE CB_ALL
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_cb   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_cb  
