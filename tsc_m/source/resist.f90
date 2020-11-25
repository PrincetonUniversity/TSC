      subroutine resist
!.....3.41 xptcalc
!
      USE CLINAM
      USE SAPROP
      USE SCR1
      USE SCR15
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     dimension difcur(pnx,pnz), volt2d(pnx,pnz)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER l,lsav,i,j,iimax,ii,kk,k,ll
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 denom,ffact,acfw,dcfw,etabase,facold,psitil,formf
      REAL*8 allam,rval,dens,tmid,etac,thalo,facsym,finterp,ajfac
      REAL*8 btfac,bsfac,fac1,ajbs,ajcd,ajlh,finterp2,curkm,curk
      REAL*8 ajtot,ajfw,ajec,pval,ppave,x2aves,facdia
!============
      call limpsi
      call magaxis
      call delpdef
      if(lrswtch.eq.0 .and. igone.eq.0) go to 10
!-> Note:
      do 490 l=1,ntpts-1
      lsav = l
      if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 497
  490 continue
  497 continue
!.....linear time point interpolation
      denom = tpro(lsav+1)-tpro(lsav)
      if(denom.eq.0) denom = 1._R8
      ffact = (time-tpro(lsav))/denom
      acfw = acfwd(lsav) + ffact * (acfwd(lsav+1) - acfwd(lsav))
      dcfw = dcfwd(lsav) + ffact * (dcfwd(lsav+1) - dcfwd(lsav))
!
      do 6 i=2,nxp
      do 5 j=2,nzp
      etabase = etav
      if(iexv(i,j).eq.0 .and. lrswtch.eq.0) then
      facold = (psi(i,j) - psilim)/(psimin - psilim)
      psitil = abs(1._R8-facold)
      formf = dcfw**2/( (psitil-acfw)**2 + dcfw**2)
      etabase = etah / formf
      if(etabase.gt.etav) etabase = etav
!
!          if(iexv(i+1,j) .ne. 0) etabase=2.*etabase
!          if(iexv(i-1,j) .ne. 0) etabase=2.*etabase
           if(iexv(i,j+1) .ne. 0) etabase=2._R8*etabase
!          if(iexv(i,j-1) .ne. 0) etabase=2.*etabase
      endif
      if(etay(i,j).ge.etabase) go to 4
      etay(i,j) = etay(i,j)+acoef(70)*(etabase-etay(i,j))
      go to 5
    4 etay(i,j) = etabase
    5 continue
      etay(i,1) = etay(i,2)
      if(isym.eq.1) etay(i,1) = etay(i,3)
    6 continue
      tevmax = 100._R8
      nqflag = nqflag+1
      if(nqflag.lt.40) go to 7
      nqflag = 0
      npts = nx/2
      if(npsi.gt.1) npts = npsi-1
      call qcalc
    7 continue
      return
   10 continue
!
      if(psilim.le.psimin .and. igone.eq.0) ineg=34
      if(ineg.ne.0) return
      if(isvd.ge.1) call xptcalc
!
!.................................................................
!  part 2:  define resistivity array
!.................................................................
!
      if(isurf.ne.0) go to 304
!.....for (isurf.eq.0) calculate etamin every 40 cycles to keep q>acoef(
!
      nqflag = nqflag+1
      if(nqflag.lt.40) go to 340
      nqflag = 0
      npts = nx/2
      if(npsi.gt.1) npts = npsi-1
      call qcalc
      if(times.gt.acoef(95)) qsaw = acoef(96)
!
!.....max temp for resistivity calculation
      iimax = acoef(63)*npts + 1
      if(iimax.lt.1) iimax = 1
      do 310 ii=2,npts
      l = npts+1-ii
      if(qprof(l).lt.qsaw.and.qprof(l+1).gt.qsaw)                        &  
     &go to 320
  310 continue
      tevmax = tprof(1)
      if(qprof(npts).lt.qsaw) tevmax = tprof(iimax)
      go to 330
  320 continue
      tevmax = ((qprof(l+1)-qsaw)*tprof(l)                               &  
     & +(qsaw-qprof(l))*tprof(l+1))                                      &  
     &       /(qprof(l+1)-qprof(l))
  330 continue
      if(tevmax.lt.tprof(iimax)) tevmax = tprof(iimax)
      allam = 24._R8-log(1.E-3_R8*sqrt(udsd*r0)/tevmax)
      etamin = (0.5_R8*1.03E-4_R8*allam)*tevmax**(-1.5_R8)*usdr*zeff
  340 continue
!
!.....define resistivity array
      do 551 i=2,nxp
      do 550 j=2,nzp
!.....check if outside of vacuum vessel
      if(iexv(i,j).eq.1) go to 500
!.....check if on wrong side of separatrix
      if(iexs(i,j).eq.1) go to 532
      go to 498
  532 continue
      if(abs(psi(i,j)-psilim) .lt. abs(phalo-psilim)) go to 499
      go to 500
  498 if(psi(i,j).gt.phalo) go to 500
      if(psi(i,j).gt.psilim) go to 499
      if(igone.eq.1) go to 499
!
!.....point lies in plasma
      rval = .25_R8*(roj(i,j)+roj(i+1,j)+roj(i,j+1)+roj(i+1,j+1))
      dens = udsd*rval
      tmid=.25_R8*(.5_R8*udsh)/dens*(pr(i,j)+pr(i+1,j)+pr(i,j+1)+pr(i+1,  &  
     & j+1))
      if(tmid.lt.tevv)      tmid = tevv
      if(tmid.gt.acoef(71)) tmid = acoef(71)
      allam = 24._R8-log(1.E-3_R8*sqrt(dens)/tmid)
      etac = 0.5_R8*1.03E-4_R8*allam*tmid**(-1.5_R8)*usdr*zeff
      if(etac.lt.etamin) etac = etamin
      if(etac.gt.etah) etac=etah
      etay(i,j) = etac
      go to 550
  499 etay(i,j) = etah
      if (acoef(99).gt.0.0_R8.and. phalo.gt.psilim)   then
      thalo =thalos-(thalos-tevv)*(psi(i,j)-psilim)/(phalo-psilim)
        if (thalo.gt.thalos)   thalo = thalos
        if (thalo.lt.tevv)     thalo = tevv
      etay(i,j) = etah * (thalos/thalo)**1.5_R8
                                                    endif
      go to 550
  500 continue
!
!.....point lies in vacuum
      etay(i,j) = etav
  550 continue
      etay(i,1) = etay(i,2)
      if(isym.eq.1) etay(i,1) = etay(i,3)
  551 continue
!
      go to 569
!
!.....for(isurf.ne.0) interpolate from etpara array
  304 continue
      bsitot = 0._R8
      cditot = 0._R8
      fwitot = 0._R8
      ecitot = 0._R8
      alhitot = 0._R8
      diatot = 0._R8
      do 651 i=2,nxp
      do 650 j=2,nzp
      facsym = 1._R8
      if(isym.eq.1 .and. j.gt.2) facsym=2._R8
      rjcd(i,j) = 0._R8
!.....check if outside of vacuum vessel
      if(iexv(i,j).eq.1) go to 600
!.....check if on wrong side of separatrix
      if(iexs(i,j).eq.1) go to 632
      go to 598
  632 continue
      if(abs(psi(i,j)-psilim) .lt. abs(phalo-psilim)) go to 599
      go to 600
  598 if(psi(i,j).gt.phalo) go to 600
      if(psi(i,j) .gt. psilim) go to 599
      if(igone.eq.1) go to 599
      do 640 kk=2,npsit
      k = kk
  640 if(psi(i,j) .lt. xsv2(k)) go to 630
      if(phalo .gt. psilim) go to 599
      go to 600
!
  630 finterp = (psi(i,j)-xsv2(k-1))/(xsv2(k)-xsv2(k-1))
      etay(i,j) = etpara(k-1)  + finterp*(etpara(k) - etpara(k-1))
      ajfac =      ajave(k-1)  + finterp*(ajave(k)  - ajave(k-1))
      btfac =      ajbtsq(k-1) + finterp*(ajbtsq(k) - ajbtsq(k-1))
      bsfac =      ajbsq(k-1)  + finterp*(ajbsq(k)  - ajbsq(k-1))
      fac1 = gs(i,j)*btfac/bsfac
!
      ajbs = (ajavbs(k-1) +finterp*(ajavbs(k)  - ajavbs(k-1)))*fac1
      ajcd = (ajavcd(k-1) +finterp*(ajavcd(k)  - ajavcd(k-1)))*fac1
      ajlh = (ajavlh2(k-1)+finterp*(ajavlh2(k) - ajavlh2(k-1)))*fac1
!
      do 635 ll=k,k+1
      l = ll
      if(psi(i,j).lt.xsv(l)) go to 625
  635 continue
  625 continue
      if(l.le.2) go to 645
      finterp2 = (psi(i,j)-xsv(l-1))/(xsv(l)-xsv(l-1))
      curkm = (gxmja2(l-1)-gxmja2(l-2))*rdpsi/tpi
      curk  = (gxmja2(l) - gxmja2(l-1))*rdpsi/tpi
      ajtot =(curkm + finterp2*(curk - curkm))*gs(i,j)
      ajfw = (ajavfw(l-1) + finterp2*(ajavfw(l) - ajavfw(l-1)))*fac1
      ajec = (ajavec(l-1) + finterp2*(ajavec(l) - ajavec(l-1)))*fac1
      go to 655
  645 ajtot =((gxmja2(2) - gxmja2(1))*rdpsi/tpi)*gs(i,j)
      ajfw = ajavfw(2)*fac1
      ajec = ajavec(2)*fac1
  655 continue
      call peval(psi(i,j),1+isurf,pval,ppave,i,j)
      x2aves = gs(i,j)**2*ajfac/bsfac
      facdia = (xary(i)**2 - x2aves)*ppave
!
!
      rjcd(i,j)  =  ajcd + ajfw + ajec + ajbs + ajlh  - facdia
      rjcdg(i,j) = (ajcd + ajfw + ajec + ajbs + ajlh  + x2aves*ppave)/gs(i,j)
      bsitot  = bsitot  + ajbs*deex*deez*udsi*facsym/xary(i)
      cditot  = cditot  + ajcd*deex*deez*udsi*facsym/xary(i)
      fwitot  = fwitot  + ajfw*deex*deez*udsi*facsym/xary(i)
      ecitot  = ecitot  + ajec*deex*deez*udsi*facsym/xary(i)
      alhitot = alhitot + ajlh*deex*deez*udsi*facsym/xary(i)
      diatot = diatot - facdia*deex*deez*udsi*facsym/xary(i)
!
!
      go to 649
  599 etay(i,j) = etah
      if (acoef(99).gt.0.0_R8.and. phalo.gt.psilim)   then
      thalo =thalos-(thalos-tevv)*(psi(i,j)-psilim)/(phalo-psilim)
        if (thalo.gt.thalos)   thalo = thalos
        if (thalo.lt.tevv)     thalo = tevv
      etay(i,j) = etah * (thalos/thalo)**1.5_R8
                                                    endif
      go to 650
  600 etay(i,j) = etav
      go to 650
  649 continue
  650 continue
      etay(i,1) = etay(i,2)
      if(isym.eq.1) etay(i,1) = etay(i,3)
      rjcdg(i,1) = rjcdg(i,2)
      if(isym.eq.1) rjcdg(i,1) = rjcdg(i,3)
  651 continue
  569 continue
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
