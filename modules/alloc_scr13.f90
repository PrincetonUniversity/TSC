      SUBROUTINE alloc_scr13  
      
      USE SCR13
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( a(penx,penz), STAT=istat)
        a = 0 
        ALLOCATE ( xnf(pnframe), STAT=istat)
        xnf = 0 
        ALLOCATE ( ynf(pnframe,3), STAT=istat)
        ynf = 0                                                            
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr13   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr13  
