      SUBROUTINE dealloc_loop1  
      
      USE LOOP1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_loop1   ' 
      end if
            
      return
      END SUBROUTINE dealloc_loop1  
