      SUBROUTINE alloc_scr17  
      
      USE SCR17
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr17   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr17  
