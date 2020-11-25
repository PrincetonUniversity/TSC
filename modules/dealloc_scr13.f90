      SUBROUTINE dealloc_scr13  
      
      USE SCR13
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( a, STAT=istat)
        DEALLOCATE ( xnf, STAT=istat)
        DEALLOCATE ( ynf, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr13   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr13  
