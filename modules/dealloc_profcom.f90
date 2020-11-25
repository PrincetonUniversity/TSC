      SUBROUTINE dealloc_profcom  
      
      USE PROFCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( npsitsv, STAT=istat)
        DEALLOCATE ( ihds, STAT=istat)
        DEALLOCATE ( xplt, STAT=istat)
        DEALLOCATE ( yplt, STAT=istat)
        DEALLOCATE ( scary, STAT=istat)
        DEALLOCATE ( psix, STAT=istat)
        DEALLOCATE ( timesv, STAT=istat)
        DEALLOCATE ( ymaxsf, STAT=istat)
        DEALLOCATE ( yminsf, STAT=istat)
        DEALLOCATE ( xcs, STAT=istat)
        DEALLOCATE ( ycs, STAT=istat)
        DEALLOCATE ( div, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_profcom   ' 
      end if
            
      return
      END SUBROUTINE dealloc_profcom  
