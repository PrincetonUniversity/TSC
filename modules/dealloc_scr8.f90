      SUBROUTINE dealloc_scr8  
      
      USE SCR8
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( afac, STAT=istat)
        DEALLOCATE ( bfac, STAT=istat)
        DEALLOCATE ( cfac, STAT=istat)
        DEALLOCATE ( sn, STAT=istat)
        DEALLOCATE ( cn, STAT=istat)
        DEALLOCATE ( x, STAT=istat)
        DEALLOCATE ( z, STAT=istat)
        DEALLOCATE ( derlt, STAT=istat)
        DEALLOCATE ( derrt, STAT=istat)
        DEALLOCATE ( abnd, STAT=istat)
        DEALLOCATE ( bbnd, STAT=istat)
        DEALLOCATE ( bv, STAT=istat)
        DEALLOCATE ( cv, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr8   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr8  
