      SUBROUTINE dealloc_scr14  
      
      USE SCR14
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( big1, STAT=istat)
        DEALLOCATE ( big2, STAT=istat)
        DEALLOCATE ( big3, STAT=istat)
        DEALLOCATE ( big4, STAT=istat)
        DEALLOCATE ( big5, STAT=istat)
        DEALLOCATE ( temp, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr14   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr14  
