      SUBROUTINE dealloc_iter1  
      
      USE ITER1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( rz, STAT=istat)
        DEALLOCATE ( zz, STAT=istat)
        DEALLOCATE ( bpolav, STAT=istat)
        DEALLOCATE ( bpolavr, STAT=istat)
        DEALLOCATE ( bpolavz, STAT=istat)
        DEALLOCATE ( bpolmax, STAT=istat)
        DEALLOCATE ( bpolr, STAT=istat)
        DEALLOCATE ( bpolz, STAT=istat)
        DEALLOCATE ( sumr, STAT=istat)
        DEALLOCATE ( sumz, STAT=istat)
        DEALLOCATE ( ansr, STAT=istat)
        DEALLOCATE ( ansz, STAT=istat)
        DEALLOCATE ( rfz, STAT=istat)
        DEALLOCATE ( zfz, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_iter1   ' 
      end if
            
      return
      END SUBROUTINE dealloc_iter1  
