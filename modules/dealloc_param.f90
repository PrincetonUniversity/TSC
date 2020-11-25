      SUBROUTINE dealloc_param  
      
      USE PARAM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_param   ' 
      end if
            
      return
      END SUBROUTINE dealloc_param  
