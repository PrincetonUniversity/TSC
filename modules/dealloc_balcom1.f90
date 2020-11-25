      SUBROUTINE dealloc_balcom1  
      
      USE BALCOM1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( bsqi, STAT=istat)
        DEALLOCATE ( xbal, STAT=istat)
        DEALLOCATE ( zbal, STAT=istat)
        DEALLOCATE ( delt, STAT=istat)
        DEALLOCATE ( deltp, STAT=istat)
        DEALLOCATE ( deltm, STAT=istat)
        DEALLOCATE ( delp1, STAT=istat)
        DEALLOCATE ( delm1, STAT=istat)
        DEALLOCATE ( delp, STAT=istat)
        DEALLOCATE ( del, STAT=istat)
        DEALLOCATE ( qpmh, STAT=istat)
        DEALLOCATE ( fpe, STAT=istat)
        DEALLOCATE ( fmh, STAT=istat)
        DEALLOCATE ( fph, STAT=istat)
        DEALLOCATE ( fm3, STAT=istat)
        DEALLOCATE ( ppmh, STAT=istat)
        DEALLOCATE ( gmh, STAT=istat)
        DEALLOCATE ( gpmh, STAT=istat)
        DEALLOCATE ( sum1, STAT=istat)
        DEALLOCATE ( sum2, STAT=istat)
        DEALLOCATE ( sum3, STAT=istat)
        DEALLOCATE ( sum4, STAT=istat)
        DEALLOCATE ( sum5, STAT=istat)
        DEALLOCATE ( sum6, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_balcom1   ' 
      end if
            
      return
      END SUBROUTINE dealloc_balcom1  
