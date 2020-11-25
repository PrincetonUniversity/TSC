      SUBROUTINE alloc_scr16  
      
      USE SCR16
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( itap1(pnx,pnz), STAT=istat)
        itap1 = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr16   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr16  
