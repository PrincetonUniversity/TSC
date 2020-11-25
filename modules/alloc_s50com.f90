      SUBROUTINE alloc_s50com  
      
      USE S50COM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( iv(pnx), STAT=istat)
        iv = 0 
        ALLOCATE ( itera(pnx), STAT=istat)
        itera = 0 
        ALLOCATE ( r1(pnx), STAT=istat)
        r1 = 0 
        ALLOCATE ( z1(pnx), STAT=istat)
        z1 = 0 
        ALLOCATE ( r2(pnx), STAT=istat)
        r2 = 0 
        ALLOCATE ( z2(pnx), STAT=istat)
        z2 = 0 
        ALLOCATE ( gr(pnx), STAT=istat)
        gr = 0 
        ALLOCATE ( gz(pnx), STAT=istat)
        gz = 0 
        ALLOCATE ( gg(pnx), STAT=istat)
        gg = 0 
        ALLOCATE ( grz(pnx), STAT=istat)
        grz = 0 
        ALLOCATE ( gzz(pnx), STAT=istat)
        gzz = 0 
        ALLOCATE ( grr(pnx), STAT=istat)
        grr = 0 
        ALLOCATE ( betaa(pnx), STAT=istat)
        betaa = 0 
        ALLOCATE ( wk(pnx), STAT=istat)
        wk = 0 
        ALLOCATE ( alfr(pnx), STAT=istat)
        alfr = 0 
        ALLOCATE ( alfi(pnx), STAT=istat)
        alfi = 0 
        ALLOCATE ( v(pnx), STAT=istat)
        v = 0                                                              
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_s50com   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_s50com  
