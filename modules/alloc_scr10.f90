      SUBROUTINE alloc_scr10  
      
      USE SCR10
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( prpest2(ppsi), STAT=istat)
        prpest2 = 0 
        ALLOCATE ( pppest2(ppsi), STAT=istat)
        pppest2 = 0 
        ALLOCATE ( ajpest2(ppsi), STAT=istat)
        ajpest2 = 0 
        ALLOCATE ( ppest(pnx), STAT=istat)
        ppest = 0 
        ALLOCATE ( gpest(pnx), STAT=istat)
        gpest = 0 
        ALLOCATE ( opest(2*pnz,pnx), STAT=istat)
        opest = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr10   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr10  
