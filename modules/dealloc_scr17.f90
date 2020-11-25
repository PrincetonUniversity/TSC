      SUBROUTINE dealloc_scr17  
      
      USE SCR17
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr17   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr17  
