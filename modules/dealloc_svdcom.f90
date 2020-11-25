      SUBROUTINE dealloc_svdcom  
      
      USE SVDCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( am, STAT=istat)
        DEALLOCATE ( uu, STAT=istat)
        DEALLOCATE ( vv, STAT=istat)
        DEALLOCATE ( datav, STAT=istat)
        DEALLOCATE ( sigma, STAT=istat)
        DEALLOCATE ( rv1, STAT=istat)
        DEALLOCATE ( coef, STAT=istat)
        DEALLOCATE ( rcocom, STAT=istat)
        DEALLOCATE ( zcocom, STAT=istat)
        DEALLOCATE ( rnormv, STAT=istat)
        DEALLOCATE ( znormv, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_svdcom   ' 
      end if
            
      return
      END SUBROUTINE dealloc_svdcom  
