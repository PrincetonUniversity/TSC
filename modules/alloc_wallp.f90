      SUBROUTINE alloc_wallp  
      
      USE WALLP
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( cpw(pnode), STAT=istat)
        cpw = 0 
        ALLOCATE ( tc(pnode), STAT=istat)
        tc = 0 
        ALLOCATE ( tcond(20), STAT=istat)
        tcond = 0 
        ALLOCATE ( csubp(20), STAT=istat)
        csubp = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_wallp   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_wallp  
