      SUBROUTINE dealloc_saprop  
      
      USE SAPROP
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        DEALLOCATE ( ane, STAT=istat)
        DEALLOCATE ( zeffa, STAT=istat)
        DEALLOCATE ( avez, STAT=istat)
        DEALLOCATE ( simpe, STAT=istat)
        DEALLOCATE ( te, STAT=istat)
        DEALLOCATE ( ti, STAT=istat)
        DEALLOCATE ( sradion, STAT=istat)
        DEALLOCATE ( zeffa2, STAT=istat)
        DEALLOCATE ( ajavlh, STAT=istat)
        DEALLOCATE ( ajpary, STAT=istat)
        DEALLOCATE ( ajavfw, STAT=istat)
        DEALLOCATE ( ajavec, STAT=istat)
        DEALLOCATE ( ajavlh2, STAT=istat)
        DEALLOCATE ( djdelh2, STAT=istat)
        DEALLOCATE ( djdetsc2, STAT=istat)
        DEALLOCATE ( voltlp, STAT=istat)
        DEALLOCATE ( anne, STAT=istat)
        DEALLOCATE ( animp, STAT=istat)
        DEALLOCATE ( sradiono, STAT=istat)
        DEALLOCATE ( aneimp, STAT=istat)
        DEALLOCATE ( avezimpa, STAT=istat)
        DEALLOCATE ( amassimpa, STAT=istat)
        DEALLOCATE ( amasshyda, STAT=istat)
        DEALLOCATE ( aimassa, STAT=istat)
        DEALLOCATE ( sumnia, STAT=istat)
        DEALLOCATE ( sumnqa, STAT=istat)
        DEALLOCATE ( pden, STAT=istat)
        DEALLOCATE ( ajavbs, STAT=istat)
        DEALLOCATE ( ajavcd, STAT=istat)
        DEALLOCATE ( anhy, STAT=istat)
        DEALLOCATE ( anox, STAT=istat)
        DEALLOCATE ( anca, STAT=istat)
        DEALLOCATE ( anfe, STAT=istat)
        DEALLOCATE ( anbe, STAT=istat)
        DEALLOCATE ( powtsc, STAT=istat)
        DEALLOCATE ( currtsc, STAT=istat)
        DEALLOCATE ( djdetsc, STAT=istat)
        DEALLOCATE ( djdelh, STAT=istat)
        DEALLOCATE ( sravejet, STAT=istat)
        DEALLOCATE ( sinzrec, STAT=istat)
        DEALLOCATE ( powtsci, STAT=istat)
     
      if (istat .ne. 0) then 
         print *, 'Deallocation Error : alloc_saprop   ' 
      end if
            
      return
      END SUBROUTINE dealloc_saprop  
