      SUBROUTINE alloc_balcli  
      
      USE BALCLI
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( alfa1(pnthe,ppsi), STAT=istat)
        alfa1 = 0 
        ALLOCATE ( alfa2(pnthe,ppsi), STAT=istat)
        alfa2 = 0 
        ALLOCATE ( alfav(pnthe,ppsi), STAT=istat)
        alfav = 0 
        ALLOCATE ( beta1(pnthe,ppsi), STAT=istat)
        beta1 = 0 
        ALLOCATE ( betav(pnthe,ppsi), STAT=istat)
        betav = 0 
        ALLOCATE ( gamav(pnthe,ppsi), STAT=istat)
        gamav = 0 
        ALLOCATE ( gama1(pnthe,ppsi), STAT=istat)
        gama1 = 0 
        ALLOCATE ( di(ppsi), STAT=istat)
        di = 0 
        ALLOCATE ( idi(ppsi), STAT=istat)
        idi = 0 
        ALLOCATE ( idr(ppsi), STAT=istat)
        idr = 0 
        ALLOCATE ( idn(ppsi), STAT=istat)
        idn = 0 
        ALLOCATE ( node1(ppsi), STAT=istat)
        node1 = 0 
        ALLOCATE ( idf(ppsi), STAT=istat)
        idf = 0                                                            
            if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_balcli   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_balcli  
