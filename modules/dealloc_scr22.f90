      SUBROUTINE dealloc_scr22  
      
      USE SCR22
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( kgrouplw, STAT=istat)
        DEALLOCATE ( jplexv, STAT=istat)
        DEALLOCATE ( xplw1, STAT=istat)
        DEALLOCATE ( zplw1, STAT=istat)
        DEALLOCATE ( xplw2, STAT=istat)
        DEALLOCATE ( zplw2, STAT=istat)
        DEALLOCATE ( polsegpl, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr22   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr22  
