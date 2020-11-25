      SUBROUTINE dealloc_cb  
      
      USE CB_ALL
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_cb   ' 
      end if
            
      return
      END SUBROUTINE dealloc_cb  
