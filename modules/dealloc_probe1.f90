      SUBROUTINE dealloc_probe1  
      
      USE PROBE1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_probe1   ' 
      end if
            
      return
      END SUBROUTINE dealloc_probe1  
