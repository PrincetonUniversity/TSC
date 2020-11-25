      SUBROUTINE dealloc_diagcom  
      
      USE DIAGCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_diagcom   ' 
      end if
            
      return
      END SUBROUTINE dealloc_diagcom  
