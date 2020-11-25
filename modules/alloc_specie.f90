      SUBROUTINE alloc_specie  
      
      USE SPECIE
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( nq(pchrgmx,pimp,ppsi), STAT=istat)
        nq = 0 
        ALLOCATE ( nqo(pchrgmx,pimp,ppsi), STAT=istat)
        nqo = 0 
        ALLOCATE ( dperi(pchrgmx,ppsi), STAT=istat)
        dperi = 0 
        ALLOCATE ( dpari(pchrgmx,ppsi), STAT=istat)
        dpari = 0 
        ALLOCATE ( ainz(pchrgmx,ppsi), STAT=istat)
        ainz = 0 
        ALLOCATE ( rec(pchrgmx,ppsi), STAT=istat)
        rec = 0 
        ALLOCATE ( rad(pimp,pchrgmx,ppsi), STAT=istat)
        rad = 0                                                            
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_specie   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_specie  
