      SUBROUTINE alloc_svdcom  
      
      USE SVDCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( am(2*nptsp,nordp), STAT=istat)
        am = 0 
        ALLOCATE ( uu(2*nptsp,2*nptsp), STAT=istat)
        uu = 0 
        ALLOCATE ( vv(nordp,nordp), STAT=istat)
        vv = 0 
        ALLOCATE ( datav(2*nptsp), STAT=istat)
        datav = 0 
        ALLOCATE ( sigma(nordp), STAT=istat)
        sigma = 0 
        ALLOCATE ( rv1(nordp), STAT=istat)
        rv1 = 0 
        ALLOCATE ( coef(nordp), STAT=istat)
        coef = 0 
        ALLOCATE ( rcocom(nptsp), STAT=istat)
        rcocom = 0 
        ALLOCATE ( zcocom(nptsp), STAT=istat)
        zcocom = 0 
        ALLOCATE ( rnormv(nptsp), STAT=istat)
        rnormv = 0 
        ALLOCATE ( znormv(nptsp), STAT=istat)
        znormv = 0                                                         
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_svdcom   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_svdcom  
