      SUBROUTINE dealloc_scr21  
      
      USE SCR21
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( ipath, STAT=istat)
        DEALLOCATE ( xpc1, STAT=istat)
        DEALLOCATE ( zpc1, STAT=istat)
        DEALLOCATE ( xpc2, STAT=istat)
        DEALLOCATE ( zpc2, STAT=istat)
        DEALLOCATE ( orient, STAT=istat)
        DEALLOCATE ( polsegc, STAT=istat)
        DEALLOCATE ( polwir, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr21   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr21  
