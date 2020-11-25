      SUBROUTINE dealloc_dercom  
      
      USE DERCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( factb, STAT=istat)
        DEALLOCATE ( alfafb, STAT=istat)
        DEALLOCATE ( betafb, STAT=istat)
        DEALLOCATE ( gamafb, STAT=istat)
        DEALLOCATE ( alfa0, STAT=istat)
        DEALLOCATE ( beta0, STAT=istat)
        DEALLOCATE ( gama0, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_dercom   ' 
      end if
            
      return
      END SUBROUTINE dealloc_dercom  
