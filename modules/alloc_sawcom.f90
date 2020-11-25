      SUBROUTINE alloc_sawcom  
      
      USE SAWCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_sawcom   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_sawcom  
