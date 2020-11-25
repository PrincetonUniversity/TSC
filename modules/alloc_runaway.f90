      SUBROUTINE alloc_runaway  
      
      USE RUNAWAY
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( ajpre(pnx,pnz), STAT=istat)
        ajpre = 0 
        ALLOCATE ( ajprecc(pnx,pnz), STAT=istat)
        ajprecc = 0 
        ALLOCATE ( ajpresf(ppsi), STAT=istat)
        ajpresf = 0 
        ALLOCATE ( anre(ppsi), STAT=istat)
        anre = 0 
        ALLOCATE ( sresf(ppsi), STAT=istat)
        sresf = 0 
        ALLOCATE ( etafac(ppsi), STAT=istat)
        etafac = 0 
        ALLOCATE ( adnre(ppsi), STAT=istat)
        adnre = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_runaway   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_runaway  
