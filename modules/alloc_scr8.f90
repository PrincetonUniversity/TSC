      SUBROUTINE alloc_scr8  
      
      USE SCR8
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( afac(3,0:penx-1), STAT=istat)
        afac = 0 
        ALLOCATE ( bfac(3,0:penx-1), STAT=istat)
        bfac = 0 
        ALLOCATE ( cfac(3,0:penx-1), STAT=istat)
        cfac = 0 
        ALLOCATE ( sn(3,0:penz-1,0:penz-1), STAT=istat)
        sn = 0 
        ALLOCATE ( cn(3,0:penz-1,0:penz-1), STAT=istat)
        cn = 0 
        ALLOCATE ( x(0:penx-1), STAT=istat)
        x = 0 
        ALLOCATE ( z(0:penz-1), STAT=istat)
        z = 0 
        ALLOCATE ( derlt(0:penz), STAT=istat)
        derlt = 0 
        ALLOCATE ( derrt(0:penz), STAT=istat)
        derrt = 0 
        ALLOCATE ( abnd(-1:penx), STAT=istat)
        abnd = 0 
        ALLOCATE ( bbnd(-1:penx), STAT=istat)
        bbnd = 0 
        ALLOCATE ( bv(penx), STAT=istat)
        bv = 0 
        ALLOCATE ( cv(penx), STAT=istat)
        cv = 0                                                             
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr8   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr8  
