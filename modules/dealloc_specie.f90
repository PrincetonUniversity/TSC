      SUBROUTINE dealloc_specie  
      
      USE SPECIE
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( nq, STAT=istat)
        DEALLOCATE ( nqo, STAT=istat)
        DEALLOCATE ( dperi, STAT=istat)
        DEALLOCATE ( dpari, STAT=istat)
        DEALLOCATE ( ainz, STAT=istat)
        DEALLOCATE ( rec, STAT=istat)
        DEALLOCATE ( rad, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_specie   ' 
      end if
            
      return
      END SUBROUTINE dealloc_specie  
