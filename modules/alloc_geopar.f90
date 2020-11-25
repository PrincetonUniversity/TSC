      SUBROUTINE alloc_geopar  
      
      USE GEOPAR
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_geopar   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_geopar  
