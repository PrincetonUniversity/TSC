      SUBROUTINE alloc_scr14  
      
      USE SCR14
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( big1(2*pnsave), STAT=istat)
        big1 = 0 
        ALLOCATE ( big2(2*pnsave), STAT=istat)
        big2 = 0 
        ALLOCATE ( big3(2*pnsave), STAT=istat)
        big3 = 0 
        ALLOCATE ( big4(2*pnsave), STAT=istat)
        big4 = 0 
        ALLOCATE ( big5(2*pnsave), STAT=istat)
        big5 = 0 
        ALLOCATE ( temp(pglobp), STAT=istat)
        temp = 0                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr14   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr14  
