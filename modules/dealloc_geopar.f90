      SUBROUTINE dealloc_geopar  
      
      USE GEOPAR
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_geopar   ' 
      end if
            
      return
      END SUBROUTINE dealloc_geopar  
