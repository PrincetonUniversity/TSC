      SUBROUTINE dealloc_scr5  
      
      USE SCR5
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( t6y, STAT=istat)
        DEALLOCATE ( t4y, STAT=istat)
        DEALLOCATE ( t8y, STAT=istat)
        DEALLOCATE ( v4tv, STAT=istat)
        DEALLOCATE ( v5tv, STAT=istat)
        DEALLOCATE ( v6tv, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr5   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr5  
