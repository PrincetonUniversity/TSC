      SUBROUTINE dealloc_scratch  
      
      USE SCRATCH
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( bigcom, STAT=istat)
        DEALLOCATE ( vecx, STAT=istat)
        DEALLOCATE ( vecz, STAT=istat)
        DEALLOCATE ( bo2, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scratch   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scratch  
