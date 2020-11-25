      SUBROUTINE alloc_wallcl  
      
      USE WALLCL
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( xw(pnthe,ppsi), STAT=istat)
        xw = 0 
        ALLOCATE ( zw(pnthe,ppsi), STAT=istat)
        zw = 0 
        ALLOCATE ( bmagfc(pnthe,ppsi), STAT=istat)
        bmagfc = 0 
        ALLOCATE ( sqg(pnthe,ppsi), STAT=istat)
        sqg = 0 
        ALLOCATE ( g33(pnthew,ppsi), STAT=istat)
        g33 = 0 
        ALLOCATE ( g11(pnthew,ppsi), STAT=istat)
        g11 = 0 
        ALLOCATE ( g22(pnthew,ppsi), STAT=istat)
        g22 = 0 
        ALLOCATE ( g21(pnthew,ppsi), STAT=istat)
        g21 = 0 
        ALLOCATE ( g12(pnthew,ppsi), STAT=istat)
        g12 = 0 
        ALLOCATE ( xcentfc(pnthe,ppsi), STAT=istat)
        xcentfc = 0 
        ALLOCATE ( apzone(pnthe,ppsi), STAT=istat)
        apzone = 0 
        ALLOCATE ( aplost(pnthe,ppsi), STAT=istat)
        aplost = 0 
        ALLOCATE ( zcentfc(pnthe,ppsi), STAT=istat)
        zcentfc = 0 
        ALLOCATE ( ripcon(pnthe,ppsi), STAT=istat)
        ripcon = 0 
        ALLOCATE ( zws(pnthe,ppsi), STAT=istat)
        zws = 0 
        ALLOCATE ( aplostr(pnthe,ppsi), STAT=istat)
        aplostr = 0                                                        
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_wallcl   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_wallcl  
