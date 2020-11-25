      SUBROUTINE alloc_scr6  
      
      USE SCR6
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( jtemp(pncoil), STAT=istat)
        jtemp = 0 
        ALLOCATE ( wrk1(2*pngroup), STAT=istat)
        wrk1 = 0 
        ALLOCATE ( wrk3(2*pngroup), STAT=istat)
        wrk3 = 0 
        ALLOCATE ( wrk2(pngroup), STAT=istat)
        wrk2 = 0 
        ALLOCATE ( xtmp1(pncoil), STAT=istat)
        xtmp1 = 0 
        ALLOCATE ( ztmp1(pncoil), STAT=istat)
        ztmp1 = 0 
        ALLOCATE ( xtmp2(pncoil), STAT=istat)
        xtmp2 = 0 
        ALLOCATE ( ztmp2(pncoil), STAT=istat)
        ztmp2 = 0 
        ALLOCATE ( atrn1(pncoil), STAT=istat)
        atrn1 = 0 
        ALLOCATE ( atrn2(pncoil), STAT=istat)
        atrn2 = 0 
        ALLOCATE ( dxtmp1(pncoil), STAT=istat)
        dxtmp1 = 0 
        ALLOCATE ( dztmp1(pncoil), STAT=istat)
        dztmp1 = 0 
        ALLOCATE ( dxtmp2(pncoil), STAT=istat)
        dxtmp2 = 0 
        ALLOCATE ( dztmp2(pncoil), STAT=istat)
        dztmp2 = 0 
        ALLOCATE ( g4(pncoil), STAT=istat)
        g4 = 0 
        ALLOCATE ( g5(pncoil), STAT=istat)
        g5 = 0 
        ALLOCATE ( g6(pncoil), STAT=istat)
        g6 = 0 
        ALLOCATE ( g1(pncoil), STAT=istat)
        g1 = 0 
        ALLOCATE ( g2(pncoil), STAT=istat)
        g2 = 0 
        ALLOCATE ( g3(pncoil), STAT=istat)
        g3 = 0 
        ALLOCATE ( atans(pngroup,pncoil), STAT=istat)
        atans = 0 
        ALLOCATE ( a1(pngroup,pngroup), STAT=istat)
        a1 = 0 
        ALLOCATE ( a1i(pngroup,pngroup), STAT=istat)
        a1i = 0 
        ALLOCATE ( a2(pngroup,pngroup), STAT=istat)
        a2 = 0 
        ALLOCATE ( a3(pngroup,pngroup), STAT=istat)
        a3 = 0 
        ALLOCATE ( a4(pngroup,pngroup), STAT=istat)
        a4 = 0 
        ALLOCATE ( rr(pngroup), STAT=istat)
        rr = 0 
        ALLOCATE ( vdt(pngroup), STAT=istat)
        vdt = 0 
        ALLOCATE ( veg0(pngroup), STAT=istat)
        veg0 = 0 
        ALLOCATE ( vdiff(pngroup), STAT=istat)
        vdiff = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr6   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr6  
