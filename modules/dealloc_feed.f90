      SUBROUTINE dealloc_feed  
      
      USE FEED
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( viabturn, STAT=istat)
        DEALLOCATE ( viabturno, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_feed   ' 
      end if
            
      return
      END SUBROUTINE dealloc_feed  
