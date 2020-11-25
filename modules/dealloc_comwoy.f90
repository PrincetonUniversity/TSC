      SUBROUTINE dealloc_comwoy  
      
      USE COMWOY
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
      DEALLOCATE ( ipow, STAT=istat)
      DEALLOCATE ( ipowm1, STAT=istat)
      DEALLOCATE ( ipowm2, STAT=istat)
      DEALLOCATE ( ipowp, STAT=istat)
      DEALLOCATE ( flxp, STAT=istat)
      DEALLOCATE ( vegint, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_comwoy   ' 
      end if
            
      return
      END SUBROUTINE dealloc_comwoy  
