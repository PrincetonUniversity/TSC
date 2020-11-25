      SUBROUTINE dealloc_newplot  
      
      USE NEWPLOT
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( chienca, STAT=istat)
        DEALLOCATE ( chiinca, STAT=istat)
        DEALLOCATE ( chiicopi, STAT=istat)
        DEALLOCATE ( chiecopi, STAT=istat)
        DEALLOCATE ( diffary, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_newplot   ' 
      end if
            
      return
      END SUBROUTINE dealloc_newplot  
