      SUBROUTINE alloc_glf  
      
      USE GLF
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_glf   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_glf  
