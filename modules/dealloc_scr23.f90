      SUBROUTINE dealloc_scr23  
      
      USE SCR23
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( vpolsegpl, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr23   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr23  
