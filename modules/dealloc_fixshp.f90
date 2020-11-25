      SUBROUTINE dealloc_fixshp  
      
      USE FIXSHP
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( amatsc, STAT=istat)
        DEALLOCATE ( uuscfb, STAT=istat)
        DEALLOCATE ( vvscfb, STAT=istat)
        DEALLOCATE ( bvecsc, STAT=istat)
        DEALLOCATE ( gsumg, STAT=istat)
        DEALLOCATE ( grsumg, STAT=istat)
        DEALLOCATE ( gzsumg, STAT=istat)
        DEALLOCATE ( grrsumg, STAT=istat)
        DEALLOCATE ( grho, STAT=istat)
        DEALLOCATE ( grho2, STAT=istat)
        DEALLOCATE ( sigma, STAT=istat)
        DEALLOCATE ( xvecsc, STAT=istat)
        DEALLOCATE ( xvecsco, STAT=istat)
        DEALLOCATE ( bvecsco, STAT=istat)
        DEALLOCATE ( rv1scfb, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_fixshp   ' 
      end if
            
      return
      END SUBROUTINE dealloc_fixshp  
