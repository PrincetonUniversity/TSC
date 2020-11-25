      SUBROUTINE alloc_scr4  
      
      USE SCR4
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( ama6(pnx), STAT=istat)
        ama6 = 0 
        ALLOCATE ( amc4(pnx), STAT=istat)
        amc4 = 0 
        ALLOCATE ( gze6a(pnx), STAT=istat)
        gze6a = 0 
        ALLOCATE ( gze5a(pnx), STAT=istat)
        gze5a = 0 
        ALLOCATE ( gze4c(pnx), STAT=istat)
        gze4c = 0 
        ALLOCATE ( gze5c(pnx), STAT=istat)
        gze5c = 0 
        ALLOCATE ( savecur(pncoil), STAT=istat)
        savecur = 0 
        ALLOCATE ( t6y(pnz), STAT=istat)
        t6y = 0 
        ALLOCATE ( t4y(pnz), STAT=istat)
        t4y = 0 
        ALLOCATE ( t8y(pnz), STAT=istat)
        t8y = 0 
        ALLOCATE ( v4tv(pnx), STAT=istat)
        v4tv = 0 
        ALLOCATE ( v5tv(pnx), STAT=istat)
        v5tv = 0 
        ALLOCATE ( v6tv(pnx), STAT=istat)
        v6tv = 0                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr4   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr4  
