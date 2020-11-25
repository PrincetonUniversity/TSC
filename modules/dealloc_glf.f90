      SUBROUTINE dealloc_glf  
      
      USE GLF
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_glf   ' 
      end if
            
      return
      END SUBROUTINE dealloc_glf  
