      SUBROUTINE alloc_scr5  
      
      USE SCR5
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
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
         print *, 'Allocation Error : alloc_scr5   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr5  
