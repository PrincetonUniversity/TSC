      SUBROUTINE dealloc_scr6  
      
      USE SCR6
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( jtemp, STAT=istat)
        DEALLOCATE ( wrk1, STAT=istat)
        DEALLOCATE ( wrk3, STAT=istat)
        DEALLOCATE ( wrk2, STAT=istat)
        DEALLOCATE ( xtmp1, STAT=istat)
        DEALLOCATE ( ztmp1, STAT=istat)
        DEALLOCATE ( xtmp2, STAT=istat)
        DEALLOCATE ( ztmp2, STAT=istat)
        DEALLOCATE ( atrn1, STAT=istat)
        DEALLOCATE ( atrn2, STAT=istat)
        DEALLOCATE ( dxtmp1, STAT=istat)
        DEALLOCATE ( dztmp1, STAT=istat)
        DEALLOCATE ( dxtmp2, STAT=istat)
        DEALLOCATE ( dztmp2, STAT=istat)
        DEALLOCATE ( g4, STAT=istat)
        DEALLOCATE ( g5, STAT=istat)
        DEALLOCATE ( g6, STAT=istat)
        DEALLOCATE ( g1, STAT=istat)
        DEALLOCATE ( g2, STAT=istat)
        DEALLOCATE ( g3, STAT=istat)
        DEALLOCATE ( atans, STAT=istat)
        DEALLOCATE ( a1, STAT=istat)
        DEALLOCATE ( a1i, STAT=istat)
        DEALLOCATE ( a2, STAT=istat)
        DEALLOCATE ( a3, STAT=istat)
        DEALLOCATE ( a4, STAT=istat)
        DEALLOCATE ( rr, STAT=istat)
        DEALLOCATE ( vdt, STAT=istat)
        DEALLOCATE ( veg0, STAT=istat)
        DEALLOCATE ( vdiff, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr6   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr6  
