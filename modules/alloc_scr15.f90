      SUBROUTINE alloc_scr15  
      
      USE SCR15
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( jvvlo(pnx), STAT=istat)
        jvvlo = 0 
        ALLOCATE ( jvvhi(pnx), STAT=istat)
        jvvhi = 0 
        ALLOCATE ( ivvlo(pnz), STAT=istat)
        ivvlo = 0 
        ALLOCATE ( ivvhi(pnz), STAT=istat)
        ivvhi = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr15   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr15  
