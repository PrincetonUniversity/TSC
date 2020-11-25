      SUBROUTINE alloc_newplot  
      
      USE NEWPLOT
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( chienca(ppsi), STAT=istat)
        chienca = 0 
        ALLOCATE ( chiinca(ppsi), STAT=istat)
        chiinca = 0 
        ALLOCATE ( chiicopi(ppsi), STAT=istat)
        chiicopi = 0 
        ALLOCATE ( chiecopi(ppsi), STAT=istat)
        chiecopi = 0 
        ALLOCATE ( diffary(ppsi), STAT=istat)
        diffary = 0                                                        
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_newplot   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_newplot  
