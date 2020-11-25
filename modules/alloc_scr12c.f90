      SUBROUTINE alloc_scr12c  
      
      USE SCR12C
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( uplot(ppsi), STAT=istat)
        uplot = 0 
        ALLOCATE ( eplot(ppsi), STAT=istat)
        eplot = 0 
        ALLOCATE ( bplot(ppsi), STAT=istat)
        bplot = 0 
        ALLOCATE ( cplot(ppsi), STAT=istat)
        cplot = 0 
        ALLOCATE ( xxplot(ppsi), STAT=istat)
        xxplot = 0 
        ALLOCATE ( tplot(ppsi), STAT=istat)
        tplot = 0 
        ALLOCATE ( uplot2(ppsi), STAT=istat)
        uplot2 = 0 
        ALLOCATE ( eplot2(ppsi), STAT=istat)
        eplot2 = 0 
        ALLOCATE ( bplot2(ppsi), STAT=istat)
        bplot2 = 0                                                         
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr12c   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr12c  
