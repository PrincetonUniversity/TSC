      SUBROUTINE alloc_saprop  
      
      USE SAPROP
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( ane(ppsi), STAT=istat)
        ane = 0 
        ALLOCATE ( zeffa(ppsi), STAT=istat)
        zeffa = 0 
        ALLOCATE ( avez(ppsi), STAT=istat)
        avez = 0 
        ALLOCATE ( simpe(ppsi), STAT=istat)
        simpe = 0 
        ALLOCATE ( te(ppsi), STAT=istat)
        te = 0 
        ALLOCATE ( ti(ppsi), STAT=istat)
        ti = 0 
        ALLOCATE ( sradion(ppsi), STAT=istat)
        sradion = 0 
        ALLOCATE ( zeffa2(ppsi), STAT=istat)
        zeffa2 = 0 
        ALLOCATE ( ajavlh(ppsi), STAT=istat)
        ajavlh = 0 
        ALLOCATE ( ajpary(ppsi), STAT=istat)
        ajpary = 0
        ALLOCATE ( ajavfw(ppsi), STAT=istat)
        ajavfw = 0 
        ALLOCATE ( ajavec(ppsi), STAT=istat)
        ajavec = 0 
        ALLOCATE ( ajavlh2(ppsi), STAT=istat)
        ajavlh2 = 0 
        ALLOCATE ( djdelh2(ppsi), STAT=istat)
        djdelh2 = 0 
        ALLOCATE ( djdetsc2(ppsi), STAT=istat)
        djdetsc2 = 0 
        ALLOCATE ( voltlp(ppsi), STAT=istat)
        voltlp = 0 
        ALLOCATE ( anne(ppsi), STAT=istat)
        anne = 0 
        ALLOCATE ( animp(ppsi), STAT=istat)
        animp = 0 
        ALLOCATE ( sradiono(ppsi), STAT=istat)
        sradiono = 0 
        ALLOCATE ( aneimp(ppsi), STAT=istat)
        aneimp = 0 
        ALLOCATE ( avezimpa(ppsi), STAT=istat)
        avezimpa = 0 
        ALLOCATE ( amassimpa(ppsi), STAT=istat)
        amassimpa = 0 
        ALLOCATE ( amasshyda(ppsi), STAT=istat)
        amasshyda = 0 
        ALLOCATE ( aimassa(ppsi), STAT=istat)
        aimassa = 0 
        ALLOCATE ( sumnia(ppsi), STAT=istat)
        sumnia = 0 
        ALLOCATE ( sumnqa(ppsi), STAT=istat)
        sumnqa = 0 
        ALLOCATE ( pden(ppsi), STAT=istat)
        pden = 0 
        ALLOCATE ( ajavbs(ppsi), STAT=istat)
        ajavbs = 0 
        ALLOCATE ( ajavcd(ppsi), STAT=istat)
        ajavcd = 0 
        ALLOCATE ( anhy(ppsi), STAT=istat)
        anhy = 0 
        ALLOCATE ( anox(ppsi), STAT=istat)
        anox = 0 
        ALLOCATE ( anca(ppsi), STAT=istat)
        anca = 0 
        ALLOCATE ( anfe(ppsi), STAT=istat)
        anfe = 0 
        ALLOCATE ( anbe(ppsi), STAT=istat)
        anbe = 0 
        ALLOCATE ( powtsc(ppsi), STAT=istat)
        powtsc = 0 
        ALLOCATE ( currtsc(ppsi), STAT=istat)
        currtsc = 0 
        ALLOCATE ( djdetsc(ppsi), STAT=istat)
        djdetsc = 0 
        ALLOCATE ( djdelh(ppsi), STAT=istat)
        djdelh = 0 
        ALLOCATE ( sravejet(ppsi), STAT=istat)
        sravejet = 0 
        ALLOCATE ( sinzrec(pnthe), STAT=istat)
        sinzrec = 0 
        ALLOCATE ( powtsci(pimp+2,ppsi), STAT=istat)
        powtsci = 0                                                        
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_saprop   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_saprop  
