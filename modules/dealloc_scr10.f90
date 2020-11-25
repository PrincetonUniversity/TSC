      SUBROUTINE dealloc_scr10  
      
      USE SCR10
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( prpest2, STAT=istat)
        DEALLOCATE ( pppest2, STAT=istat)
        DEALLOCATE ( ajpest2, STAT=istat)
        DEALLOCATE ( ppest, STAT=istat)
        DEALLOCATE ( gpest, STAT=istat)
        DEALLOCATE ( opest, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr10   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr10  
