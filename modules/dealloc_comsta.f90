      SUBROUTINE dealloc_comsta  
      
      USE COMSTA
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( tpl, STAT=istat)
        DEALLOCATE ( fstabpl, STAT=istat)
        DEALLOCATE ( fdestpl, STAT=istat)
        DEALLOCATE ( fsumpl, STAT=istat)
        DEALLOCATE ( taupl, STAT=istat)
        DEALLOCATE ( cpasopl, STAT=istat)
        DEALLOCATE ( cpasupl, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_comsta   ' 
      end if
            
      return
      END SUBROUTINE dealloc_comsta  
