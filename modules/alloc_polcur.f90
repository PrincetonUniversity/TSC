      SUBROUTINE alloc_polcur  
      
      USE POLCUR
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                            
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_polcur   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_polcur  
