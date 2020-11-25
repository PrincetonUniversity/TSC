      SUBROUTINE alloc_decom  
      
      USE DECOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_decom   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_decom  
