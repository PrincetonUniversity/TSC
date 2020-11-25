      SUBROUTINE dealloc_s50com  
      
      USE S50COM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( iv, STAT=istat)
        DEALLOCATE ( itera, STAT=istat)
        DEALLOCATE ( r1, STAT=istat)
        DEALLOCATE ( z1, STAT=istat)
        DEALLOCATE ( r2, STAT=istat)
        DEALLOCATE ( z2, STAT=istat)
        DEALLOCATE ( gr, STAT=istat)
        DEALLOCATE ( gz, STAT=istat)
        DEALLOCATE ( gg, STAT=istat)
        DEALLOCATE ( grz, STAT=istat)
        DEALLOCATE ( gzz, STAT=istat)
        DEALLOCATE ( grr, STAT=istat)
        DEALLOCATE ( betaa, STAT=istat)
        DEALLOCATE ( wk, STAT=istat)
        DEALLOCATE ( alfr, STAT=istat)
        DEALLOCATE ( alfi, STAT=istat)
        DEALLOCATE ( v, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_s50com   ' 
      end if
            
      return
      END SUBROUTINE dealloc_s50com  
