      SUBROUTINE alloc_param  
      
      USE PARAM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_param   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_param  
