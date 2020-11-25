      SUBROUTINE alloc_diagcom  
      
      USE DIAGCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_diagcom   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_diagcom  
