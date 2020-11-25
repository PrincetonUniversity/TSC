      SUBROUTINE alloc_radtab  
      
      USE RADTAB
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( alinzOx(pcgOx,pte), STAT=istat)
        alinzOx = 0 
        ALLOCATE ( alradOx(pcgOx,pte,pne), STAT=istat)
        alradOx = 0 
        ALLOCATE ( alrecOx(pcgOx,pte,pne), STAT=istat)
        alrecOx = 0 
        ALLOCATE ( alinzC(pcgC ,pte), STAT=istat)
        alinzC = 0 
        ALLOCATE ( alradC(pcgC ,pte,pne), STAT=istat)
        alradC = 0 
        ALLOCATE ( alrecC(pcgC ,pte,pne), STAT=istat)
        alrecC = 0 
        ALLOCATE ( alinzFe(pcgFe,pte), STAT=istat)
        alinzFe = 0 
        ALLOCATE ( alradFe(pcgFe,pte,pne), STAT=istat)
        alradFe = 0 
        ALLOCATE ( alrecFe(pcgFe,pte,pne), STAT=istat)
        alrecFe = 0 
        ALLOCATE ( alinzBe(pcgBe,pte), STAT=istat)
        alinzBe = 0 
        ALLOCATE ( alradBe(pcgBe,pte,pne), STAT=istat)
        alradBe = 0 
        ALLOCATE ( alrecBe(pcgBe,pte,pne), STAT=istat)
        alrecBe = 0 
        ALLOCATE ( alinzNe(pcgNe,pte), STAT=istat)
        alinzNe = 0 
        ALLOCATE ( alradNe(pcgNe,pte,pne), STAT=istat)
        alradNe = 0 
        ALLOCATE ( alrecNe(pcgNe,pte,pne), STAT=istat)
        alrecNe = 0 
        ALLOCATE ( alinzKr(pcgKr,pte), STAT=istat)
        alinzKr = 0 
        ALLOCATE ( alradKr(pcgKr,pte,pne), STAT=istat)
        alradKr = 0 
        ALLOCATE ( alrecKr(pcgKr,pte,pne), STAT=istat)
        alrecKr = 0 
        ALLOCATE ( alinzAr(pcgAr,pte), STAT=istat)
        alinzAr = 0 
        ALLOCATE ( alradAr(pcgAr,pte,pne), STAT=istat)
        alradAr = 0 
        ALLOCATE ( alrecAr(pcgAr,pte,pne), STAT=istat)
        alrecAr = 0 
        ALLOCATE ( alinzW(pcgW ,pte), STAT=istat)
        alinzW = 0 
        ALLOCATE ( alradW(pcgW ,pte,pne), STAT=istat)
        alradW = 0 
        ALLOCATE ( alrecW(pcgW ,pte,pne), STAT=istat)
        alrecW = 0 
        ALLOCATE ( altei(pte,pimp), STAT=istat)
        altei = 0 
        ALLOCATE ( alnei(pne,pimp), STAT=istat)
        alnei = 0 
        ALLOCATE ( nne(pimp), STAT=istat)
        nne = 0 
        ALLOCATE ( nte(pimp), STAT=istat)
        nte = 0 
        ALLOCATE ( nchrgsr(pimp), STAT=istat)
        nchrgsr = 0                                                        
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_radtab   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_radtab  
