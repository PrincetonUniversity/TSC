      SUBROUTINE dealloc_noncor  
      
      USE NONCOR
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( nzspec, STAT=istat)
        DEALLOCATE ( fracim, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_noncor   ' 
      end if
            
      return
      END SUBROUTINE dealloc_noncor  
