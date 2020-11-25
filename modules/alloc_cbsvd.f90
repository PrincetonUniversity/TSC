      SUBROUTINE alloc_cbsvd  
      
      USE CBSVD
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                            
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_cbsvd   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_cbsvd  
