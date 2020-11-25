      SUBROUTINE dealloc_scadvan  
      
      USE SCADVAN
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( iforce, STAT=istat)
        DEALLOCATE ( dxa, STAT=istat)
        DEALLOCATE ( dxc, STAT=istat)
        DEALLOCATE ( dzb, STAT=istat)
        DEALLOCATE ( dzd, STAT=istat)
        DEALLOCATE ( face, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scadvan   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scadvan  
