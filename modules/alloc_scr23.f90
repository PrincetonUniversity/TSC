      SUBROUTINE alloc_scr23  
      
      USE SCR23
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( vpolsegpl(2*pncoil), STAT=istat)
        vpolsegpl = 0                                                      
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr23   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr23  
