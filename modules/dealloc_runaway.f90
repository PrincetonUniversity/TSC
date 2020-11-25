      SUBROUTINE dealloc_runaway  
      
      USE RUNAWAY
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( ajpre, STAT=istat)
        DEALLOCATE ( ajprecc, STAT=istat)
        DEALLOCATE ( ajpresf, STAT=istat)
        DEALLOCATE ( anre, STAT=istat)
        DEALLOCATE ( sresf, STAT=istat)
        DEALLOCATE ( etafac, STAT=istat)
        DEALLOCATE ( adnre, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_runaway   ' 
      end if
            
      return
      END SUBROUTINE dealloc_runaway  
