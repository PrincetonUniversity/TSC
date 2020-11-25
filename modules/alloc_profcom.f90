      SUBROUTINE alloc_profcom  
      
      USE PROFCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( npsitsv(50), STAT=istat)
        npsitsv = 0 
        ALLOCATE ( ihds(4), STAT=istat)
        ihds = 0 
        ALLOCATE ( xplt(pnx+ppsi+pnseg), STAT=istat)
        xplt = 0 
        ALLOCATE ( yplt(pnx+ppsi+pnseg), STAT=istat)
        yplt = 0 
        ALLOCATE ( scary(pnplts+10,pnx+ppsi+pnseg), STAT=istat)
        scary = 0 
        ALLOCATE ( psix(pnx), STAT=istat)
        psix = 0 
        ALLOCATE ( timesv(50), STAT=istat)
        timesv = 0 
        ALLOCATE ( ymaxsf(500), STAT=istat)
        ymaxsf = 0 
        ALLOCATE ( yminsf(500), STAT=istat)
        yminsf = 0 
        ALLOCATE ( xcs(4), STAT=istat)
        xcs = 0 
        ALLOCATE ( ycs(4), STAT=istat)
        ycs = 0 
        ALLOCATE ( div(2), STAT=istat)
        div = 0                                                            
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_profcom   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_profcom  
