      SUBROUTINE alloc_scr22  
      
      USE SCR22
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( kgrouplw(2*pncoil), STAT=istat)
        kgrouplw = 0 
        ALLOCATE ( jplexv(2*pncoil), STAT=istat)
        jplexv = 0 
        ALLOCATE ( xplw1(2*pncoil), STAT=istat)
        xplw1 = 0 
        ALLOCATE ( zplw1(2*pncoil), STAT=istat)
        zplw1 = 0 
        ALLOCATE ( xplw2(2*pncoil), STAT=istat)
        xplw2 = 0 
        ALLOCATE ( zplw2(2*pncoil), STAT=istat)
        zplw2 = 0 
        ALLOCATE ( polsegpl(2*pncoil), STAT=istat)
        polsegpl = 0                                                       
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr22   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr22  
