      SUBROUTINE dealloc_scr7  
      
      USE SCR7
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( aa1, STAT=istat)
        DEALLOCATE ( dert, STAT=istat)
        DEALLOCATE ( derb, STAT=istat)
        DEALLOCATE ( derl, STAT=istat)
        DEALLOCATE ( derr, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr7   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr7  
