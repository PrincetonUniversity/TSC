      SUBROUTINE dealloc_worka  
      
      USE WORKA
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( aw, STAT=istat)
        DEALLOCATE ( bw, STAT=istat)
        DEALLOCATE ( cw, STAT=istat)
        DEALLOCATE ( d, STAT=istat)
        DEALLOCATE ( e, STAT=istat)
        DEALLOCATE ( f, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_worka   ' 
      end if
            
      return
      END SUBROUTINE dealloc_worka  
