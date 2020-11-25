      SUBROUTINE dealloc_scr11  
      
      USE SCR11
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( grsum, STAT=istat)
        DEALLOCATE ( grsum0, STAT=istat)
        DEALLOCATE ( gcurr, STAT=istat)
        DEALLOCATE ( gvsum, STAT=istat)
        DEALLOCATE ( gvsum0, STAT=istat)
        DEALLOCATE ( gcur0ka, STAT=istat)
        DEALLOCATE ( gcurfka, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr11   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr11  
