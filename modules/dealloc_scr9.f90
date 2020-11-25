      SUBROUTINE dealloc_scr9  
      
      USE SCR9
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( izer, STAT=istat)
        DEALLOCATE ( jzer, STAT=istat)
        DEALLOCATE ( ans1v, STAT=istat)
        DEALLOCATE ( sns1v, STAT=istat)
        DEALLOCATE ( ans2v, STAT=istat)
        DEALLOCATE ( sns2v, STAT=istat)
        DEALLOCATE ( r1, STAT=istat)
        DEALLOCATE ( z1, STAT=istat)
        DEALLOCATE ( r2, STAT=istat)
        DEALLOCATE ( z2, STAT=istat)
        DEALLOCATE ( g1, STAT=istat)
        DEALLOCATE ( g2, STAT=istat)
        DEALLOCATE ( g3, STAT=istat)
        DEALLOCATE ( g4, STAT=istat)
        DEALLOCATE ( g5, STAT=istat)
        DEALLOCATE ( g6, STAT=istat)
        DEALLOCATE ( plgf1, STAT=istat)
        DEALLOCATE ( plgf2, STAT=istat)
        DEALLOCATE ( signzer, STAT=istat)
        DEALLOCATE ( sumi, STAT=istat)
        DEALLOCATE ( sumj, STAT=istat)
        DEALLOCATE ( gfun4, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr9   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr9  
