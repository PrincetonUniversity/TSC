      SUBROUTINE alloc_scr7  
      
      USE SCR7
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( aa1(penx,penz), STAT=istat)
        aa1 = 0 
        ALLOCATE ( dert(0:penx), STAT=istat)
        dert = 0 
        ALLOCATE ( derb(0:penx), STAT=istat)
        derb = 0 
        ALLOCATE ( derl(0:penz), STAT=istat)
        derl = 0 
        ALLOCATE ( derr(0:penz), STAT=istat)
        derr = 0                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr7   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr7  
