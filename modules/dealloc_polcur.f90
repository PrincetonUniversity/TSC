      SUBROUTINE dealloc_polcur  
      
      USE POLCUR
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_polcur   ' 
      end if
            
      return
      END SUBROUTINE dealloc_polcur  
