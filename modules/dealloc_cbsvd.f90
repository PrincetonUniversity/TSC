      SUBROUTINE dealloc_cbsvd  
      
      USE CBSVD
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_cbsvd   ' 
      end if
            
      return
      END SUBROUTINE dealloc_cbsvd  
