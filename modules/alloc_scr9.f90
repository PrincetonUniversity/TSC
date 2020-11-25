      SUBROUTINE alloc_scr9  
      
      USE SCR9
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( izer(pboun), STAT=istat)
        izer = 0 
        ALLOCATE ( jzer(pboun), STAT=istat)
        jzer = 0 
        ALLOCATE ( ans1v(penx,2,pnwire,2), STAT=istat)
        ans1v = 0 
        ALLOCATE ( sns1v(penx,2,pnwire,2), STAT=istat)
        sns1v = 0 
        ALLOCATE ( ans2v(penz,2,pnwire,2), STAT=istat)
        ans2v = 0 
        ALLOCATE ( sns2v(penz,2,pnwire,2), STAT=istat)
        sns2v = 0 
        ALLOCATE ( r1(pboun), STAT=istat)
        r1 = 0 
        ALLOCATE ( z1(pboun), STAT=istat)
        z1 = 0 
        ALLOCATE ( r2(pboun), STAT=istat)
        r2 = 0 
        ALLOCATE ( z2(pboun), STAT=istat)
        z2 = 0 
        ALLOCATE ( g1(pboun), STAT=istat)
        g1 = 0 
        ALLOCATE ( g2(pboun), STAT=istat)
        g2 = 0 
        ALLOCATE ( g3(pboun), STAT=istat)
        g3 = 0 
        ALLOCATE ( g4(pboun), STAT=istat)
        g4 = 0 
        ALLOCATE ( g5(pboun), STAT=istat)
        g5 = 0 
        ALLOCATE ( g6(pboun), STAT=istat)
        g6 = 0 
        ALLOCATE ( plgf1(pboun), STAT=istat)
        plgf1 = 0 
        ALLOCATE ( plgf2(pboun), STAT=istat)
        plgf2 = 0 
        ALLOCATE ( signzer(pboun), STAT=istat)
        signzer = 0 
        ALLOCATE ( sumi(pnx), STAT=istat)
        sumi = 0 
        ALLOCATE ( sumj(pnz), STAT=istat)
        sumj = 0 
        ALLOCATE ( gfun4(pboun,pboun), STAT=istat)
        gfun4 = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr9   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr9  
