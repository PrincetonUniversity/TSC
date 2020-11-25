      SUBROUTINE dealloc_scr15  
      
      USE SCR15
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( jvvlo, STAT=istat)
        DEALLOCATE ( jvvhi, STAT=istat)
        DEALLOCATE ( ivvlo, STAT=istat)
        DEALLOCATE ( ivvhi, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr15   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr15  
