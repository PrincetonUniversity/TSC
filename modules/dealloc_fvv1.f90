      SUBROUTINE dealloc_fvv1  
      
      USE FVV1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( bpolr, STAT=istat)
        DEALLOCATE ( bpolz, STAT=istat)
        DEALLOCATE ( fwirer, STAT=istat)
        DEALLOCATE ( fwirez, STAT=istat)
        DEALLOCATE ( vvdum, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_fvv1   ' 
      end if
            
      return
      END SUBROUTINE dealloc_fvv1  
