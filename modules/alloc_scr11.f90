      SUBROUTINE alloc_scr11  
      
      USE SCR11
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( grsum(pngroup), STAT=istat)
        grsum = 0 
        ALLOCATE ( grsum0(pngroup), STAT=istat)
        grsum0 = 0 
        ALLOCATE ( gcurr(pngroup), STAT=istat)
        gcurr = 0 
        ALLOCATE ( gvsum(pngroup), STAT=istat)
        gvsum = 0 
        ALLOCATE ( gvsum0(pngroup), STAT=istat)
        gvsum0 = 0 
        ALLOCATE ( gcur0ka(pngroup), STAT=istat)
        gcur0ka = 0 
        ALLOCATE ( gcurfka(pngroup), STAT=istat)
        gcurfka = 0                                                        
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr11   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr11  
