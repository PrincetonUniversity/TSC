      SUBROUTINE alloc_cintg  
      
      USE CINTG
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( ahat(201,201), STAT=istat)
        ahat = 0 
        ALLOCATE ( bhat(201,7), STAT=istat)
        bhat = 0 
        ALLOCATE ( fhat(201,7), STAT=istat)
        fhat = 0 
        ALLOCATE ( chat(7,201), STAT=istat)
        chat = 0 
        ALLOCATE ( gapdw(6), STAT=istat)
        gapdw = 0 
        ALLOCATE ( dgapdw(6), STAT=istat)
        dgapdw = 0                                                         
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_cintg   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_cintg  
