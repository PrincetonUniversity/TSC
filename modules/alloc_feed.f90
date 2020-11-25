      SUBROUTINE alloc_feed  
      
      USE FEED
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( viabturn(18), STAT=istat)
        viabturn = 0 
        ALLOCATE ( viabturno(18), STAT=istat)
        viabturno = 0                                                      
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_feed   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_feed  
