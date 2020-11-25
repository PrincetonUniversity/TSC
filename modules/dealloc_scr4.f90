      SUBROUTINE dealloc_scr4  
      
      USE SCR4
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( ama6, STAT=istat)
        DEALLOCATE ( amc4, STAT=istat)
        DEALLOCATE ( gze6a, STAT=istat)
        DEALLOCATE ( gze5a, STAT=istat)
        DEALLOCATE ( gze4c, STAT=istat)
        DEALLOCATE ( gze5c, STAT=istat)
        DEALLOCATE ( savecur, STAT=istat)
        DEALLOCATE ( t6y, STAT=istat)
        DEALLOCATE ( t4y, STAT=istat)
        DEALLOCATE ( t8y, STAT=istat)
        DEALLOCATE ( v4tv, STAT=istat)
        DEALLOCATE ( v5tv, STAT=istat)
        DEALLOCATE ( v6tv, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr4   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr4  
