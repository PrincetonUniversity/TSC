      SUBROUTINE alloc_worka  
      
      USE WORKA
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( aw(pnthe), STAT=istat)
        aw = 0 
        ALLOCATE ( bw(pnthe), STAT=istat)
        bw = 0 
        ALLOCATE ( cw(pnthe), STAT=istat)
        cw = 0 
        ALLOCATE ( d(pnthe), STAT=istat)
        d = 0 
        ALLOCATE ( e(pnthe), STAT=istat)
        e = 0 
        ALLOCATE ( f(pnthe), STAT=istat)
        f = 0                                                              
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_worka   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_worka  
