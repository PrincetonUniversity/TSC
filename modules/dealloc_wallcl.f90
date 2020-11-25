      SUBROUTINE dealloc_wallcl  
      
      USE WALLCL
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( xw, STAT=istat)
        DEALLOCATE ( zw, STAT=istat)
        DEALLOCATE ( bmagfc, STAT=istat)
        DEALLOCATE ( sqg, STAT=istat)
        DEALLOCATE ( g33, STAT=istat)
        DEALLOCATE ( g11, STAT=istat)
        DEALLOCATE ( g22, STAT=istat)
        DEALLOCATE ( g21, STAT=istat)
        DEALLOCATE ( g12, STAT=istat)
        DEALLOCATE ( xcentfc, STAT=istat)
        DEALLOCATE ( apzone, STAT=istat)
        DEALLOCATE ( aplost, STAT=istat)
        DEALLOCATE ( zcentfc, STAT=istat)
        DEALLOCATE ( ripcon, STAT=istat)
        DEALLOCATE ( zws, STAT=istat)
        DEALLOCATE ( aplostr, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_wallcl   ' 
      end if
            
      return
      END SUBROUTINE dealloc_wallcl  
