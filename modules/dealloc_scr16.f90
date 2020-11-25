      SUBROUTINE dealloc_scr16  
      
      USE SCR16
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( itap1, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr16   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr16  
