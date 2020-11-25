      SUBROUTINE alloc_noncor  
      
      USE NONCOR
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( nzspec(3), STAT=istat)
        nzspec = 0 
        ALLOCATE ( fracim(3), STAT=istat)
        fracim = 0                                                         
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_noncor   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_noncor  
