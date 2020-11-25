      SUBROUTINE dealloc_cintg  
      
      USE CINTG
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( ahat, STAT=istat)
        DEALLOCATE ( bhat, STAT=istat)
        DEALLOCATE ( fhat, STAT=istat)
        DEALLOCATE ( chat, STAT=istat)
        DEALLOCATE ( gapdw, STAT=istat)
        DEALLOCATE ( dgapdw, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_cintg   ' 
      end if
            
      return
      END SUBROUTINE dealloc_cintg  
