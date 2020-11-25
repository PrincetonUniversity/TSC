      SUBROUTINE dealloc_impcom  
      
      USE IMPCOM
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( nq, STAT=istat)
        DEALLOCATE ( ainz, STAT=istat)
        DEALLOCATE ( rec, STAT=istat)
        DEALLOCATE ( rad, STAT=istat)
        DEALLOCATE ( alinzOx, STAT=istat)
        DEALLOCATE ( alradOx, STAT=istat)
        DEALLOCATE ( alrecOx, STAT=istat)
        DEALLOCATE ( alinzC, STAT=istat)
        DEALLOCATE ( alradC, STAT=istat)
        DEALLOCATE ( alrecC, STAT=istat)
        DEALLOCATE ( alinzFe, STAT=istat)
        DEALLOCATE ( alradFe, STAT=istat)
        DEALLOCATE ( alrecFe, STAT=istat)
        DEALLOCATE ( alinzBe, STAT=istat)
        DEALLOCATE ( alradBe, STAT=istat)
        DEALLOCATE ( alrecBe, STAT=istat)
        DEALLOCATE ( alinzNe, STAT=istat)
        DEALLOCATE ( alradNe, STAT=istat)
        DEALLOCATE ( alrecNe, STAT=istat)
        DEALLOCATE ( alinzKr, STAT=istat)
        DEALLOCATE ( alradKr, STAT=istat)
        DEALLOCATE ( alrecKr, STAT=istat)
        DEALLOCATE ( alinzAr, STAT=istat)
        DEALLOCATE ( alradAr, STAT=istat)
        DEALLOCATE ( alrecAr, STAT=istat)
        DEALLOCATE ( alinzW, STAT=istat)
        DEALLOCATE ( alradW, STAT=istat)
        DEALLOCATE ( alrecW, STAT=istat)
        DEALLOCATE ( altei, STAT=istat)
        DEALLOCATE ( alnei, STAT=istat)
        DEALLOCATE ( nne, STAT=istat)
        DEALLOCATE ( nte, STAT=istat)
        DEALLOCATE ( nchrgsr, STAT=istat)
        DEALLOCATE ( snp, STAT=istat)
        DEALLOCATE ( aimp, STAT=istat)
        DEALLOCATE ( bimp, STAT=istat)
        DEALLOCATE ( cimp, STAT=istat)
        DEALLOCATE ( dimp, STAT=istat)
        DEALLOCATE ( ework, STAT=istat)
        DEALLOCATE ( fwork, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_impcom   ' 
      end if
            
      return
      END SUBROUTINE dealloc_impcom  
