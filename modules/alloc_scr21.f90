      SUBROUTINE alloc_scr21  
      
      USE SCR21
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( ipath(pnwire, 8), STAT=istat)
        ipath = 0 
        ALLOCATE ( xpc1(pncoil), STAT=istat)
        xpc1 = 0 
        ALLOCATE ( zpc1(pncoil), STAT=istat)
        zpc1 = 0 
        ALLOCATE ( xpc2(pncoil), STAT=istat)
        xpc2 = 0 
        ALLOCATE ( zpc2(pncoil), STAT=istat)
        zpc2 = 0 
        ALLOCATE ( orient(pncoil), STAT=istat)
        orient = 0 
        ALLOCATE ( polsegc(pncoil), STAT=istat)
        polsegc = 0 
        ALLOCATE ( polwir(pncoil), STAT=istat)
        polwir = 0                                                         
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr21   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr21  
