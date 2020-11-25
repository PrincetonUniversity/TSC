      SUBROUTINE alloc_dercom  
      
      USE DERCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( factb(4), STAT=istat)
        factb = 0 
        ALLOCATE ( alfafb(4,ppsi), STAT=istat)
        alfafb = 0 
        ALLOCATE ( betafb(4,ppsi), STAT=istat)
        betafb = 0 
        ALLOCATE ( gamafb(4,ppsi), STAT=istat)
        gamafb = 0 
        ALLOCATE ( alfa0(ppsi), STAT=istat)
        alfa0 = 0 
        ALLOCATE ( beta0(ppsi), STAT=istat)
        beta0 = 0 
        ALLOCATE ( gama0(ppsi), STAT=istat)
        gama0 = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_dercom   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_dercom  
