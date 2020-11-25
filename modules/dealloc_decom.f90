      SUBROUTINE dealloc_decom  
      
      USE DECOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_decom   ' 
      end if
            
      return
      END SUBROUTINE dealloc_decom  
