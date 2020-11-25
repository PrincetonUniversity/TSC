      SUBROUTINE dealloc_wallp  
      
      USE WALLP
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( cpw, STAT=istat)
        DEALLOCATE ( tc, STAT=istat)
        DEALLOCATE ( tcond, STAT=istat)
        DEALLOCATE ( csubp, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_wallp   ' 
      end if
            
      return
      END SUBROUTINE dealloc_wallp  
