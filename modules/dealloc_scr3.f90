      SUBROUTINE dealloc_scr3  
      
      USE SCR3
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( amat, STAT=istat)
        DEALLOCATE ( bmat, STAT=istat)
        DEALLOCATE ( cmat, STAT=istat)
        DEALLOCATE ( emat, STAT=istat)
        DEALLOCATE ( fvec, STAT=istat)
        DEALLOCATE ( dvec, STAT=istat)
        DEALLOCATE ( pvec, STAT=istat)
        DEALLOCATE ( tmp1, STAT=istat)
        DEALLOCATE ( tmp4, STAT=istat)
        DEALLOCATE ( tmp5, STAT=istat)
        DEALLOCATE ( tmp3, STAT=istat)
        DEALLOCATE ( xs6, STAT=istat)
        DEALLOCATE ( as6, STAT=istat)
        DEALLOCATE ( bs6, STAT=istat)
        DEALLOCATE ( cs6, STAT=istat)
        DEALLOCATE ( ds6, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_scr3   ' 
      end if
            
      return
      END SUBROUTINE dealloc_scr3  
