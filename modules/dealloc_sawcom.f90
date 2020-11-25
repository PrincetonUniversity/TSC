      SUBROUTINE dealloc_sawcom  
      
      USE SAWCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_sawcom   ' 
      end if
            
      return
      END SUBROUTINE dealloc_sawcom  
