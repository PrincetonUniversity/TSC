      SUBROUTINE dealloc_scr12a  
      
      USE SCR12A
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( uplot, STAT=istat)
        DEALLOCATE ( eplot, STAT=istat)
        DEALLOCATE ( bplot, STAT=istat)
        DEALLOCATE ( cplot, STAT=istat)
        DEALLOCATE ( xxplot, STAT=istat)
        DEALLOCATE ( tplot, STAT=istat)
        DEALLOCATE ( uplot2, STAT=istat)
        DEALLOCATE ( eplot2, STAT=istat)
        DEALLOCATE ( bplot2, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr12a   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr12a  
