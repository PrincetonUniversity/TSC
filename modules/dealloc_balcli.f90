      SUBROUTINE dealloc_balcli  
      
      USE BALCLI
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( alfa1, STAT=istat)
        DEALLOCATE ( alfa2, STAT=istat)
        DEALLOCATE ( alfav, STAT=istat)
        DEALLOCATE ( beta1, STAT=istat)
        DEALLOCATE ( betav, STAT=istat)
        DEALLOCATE ( gamav, STAT=istat)
        DEALLOCATE ( gama1, STAT=istat)
        DEALLOCATE ( di, STAT=istat)
        DEALLOCATE ( idi, STAT=istat)
        DEALLOCATE ( idr, STAT=istat)
        DEALLOCATE ( idn, STAT=istat)
        DEALLOCATE ( node1, STAT=istat)
        DEALLOCATE ( idf, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_balcli   ' 
      end if
            
      return
      END SUBROUTINE dealloc_balcli  
