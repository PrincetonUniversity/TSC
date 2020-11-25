      SUBROUTINE alloc_balcom1  
      
      USE BALCOM1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( bsqi(pnthe,ppsi), STAT=istat)
        bsqi = 0 
        ALLOCATE ( xbal(pnthe,ppsi), STAT=istat)
        xbal = 0 
        ALLOCATE ( zbal(pnthe,ppsi), STAT=istat)
        zbal = 0 
        ALLOCATE ( delt(pnthe,ppsi), STAT=istat)
        delt = 0 
        ALLOCATE ( deltp(pnthe,ppsi), STAT=istat)
        deltp = 0 
        ALLOCATE ( deltm(pnthe,ppsi), STAT=istat)
        deltm = 0 
        ALLOCATE ( delp1(pnthe,ppsi), STAT=istat)
        delp1 = 0 
        ALLOCATE ( delm1(pnthe,ppsi), STAT=istat)
        delm1 = 0 
        ALLOCATE ( delp(pnthe,ppsi), STAT=istat)
        delp = 0 
        ALLOCATE ( del(pnthe,ppsi), STAT=istat)
        del = 0 
        ALLOCATE ( qpmh(ppsi), STAT=istat)
        qpmh = 0 
        ALLOCATE ( fpe(ppsi), STAT=istat)
        fpe = 0 
        ALLOCATE ( fmh(ppsi), STAT=istat)
        fmh = 0 
        ALLOCATE ( fph(ppsi), STAT=istat)
        fph = 0 
        ALLOCATE ( fm3(ppsi), STAT=istat)
        fm3 = 0 
        ALLOCATE ( ppmh(ppsi), STAT=istat)
        ppmh = 0 
        ALLOCATE ( gmh(ppsi), STAT=istat)
        gmh = 0 
        ALLOCATE ( gpmh(ppsi), STAT=istat)
        gpmh = 0 
        ALLOCATE ( sum1(ppsi), STAT=istat)
        sum1 = 0 
        ALLOCATE ( sum2(ppsi), STAT=istat)
        sum2 = 0 
        ALLOCATE ( sum3(ppsi), STAT=istat)
        sum3 = 0 
        ALLOCATE ( sum4(ppsi), STAT=istat)
        sum4 = 0 
        ALLOCATE ( sum5(ppsi), STAT=istat)
        sum5 = 0 
        ALLOCATE ( sum6(ppsi), STAT=istat)
        sum6 = 0                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_balcom1   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_balcom1  
