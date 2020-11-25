      subroutine efield
!
      USE CLINAM
      USE COMWOY
      USE SVDCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER l,ll,llsav,li,ii,ngr,iabs,ngrvwl,ngrvcl,ig,iigr,igr,i
      INTEGER icnt8,icnt9,iw,iswtch,n,nl,ie,np,iem10,iem20,iii
      INTEGER index, lsav
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 coiturn,denom,denom1,term1,term1d,denom2,term2
      REAL*8 term2d,denom3,term3,term3d,sumturn,sumcur,fbchix
      REAL*8 fbmult,rlin,rdis,pval,rval,rpval,zden,zr,za,zip,zk,zd
      REAL*8 zp,zb,zq,ztohmic,ztau,taueka,fbchixold,excess,skal
      REAL*8 dtramp,aplfb,fluxdif,psinn,pinterp,psinp,ellip,delta
      REAL*8 x1,x2,z1,z2,psi1,psi2,psi3,psi4,psisn1,alphac,term
      REAL*8 fluxdif1,formf,fac,facrxw,fluxoff,termd,atuv,ppcur
      REAL*8 aplav,plcur,rebmax,rebmin
      REAL*8 fbchiix, fbchiixold
      REAL*8 tped_corr,fhmodeix,fhmodeixold,a3208
!============
!     dimension factx(ptpts)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: factx
!============      
      IF(.not.ALLOCATED(factx)) ALLOCATE( factx(ptpts), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : efield  ' 
!============      
!
!...number of turns of coi-coil
      coiturn = 5._R8
!.....................................................................
!.....part 1...compute desired currents and pressures
!.....................................................................
!
      do 750 l=1,ntpts
      fact(l) = 0._R8
  750 facd(l) = 0
      if(kcycle .gt. 0) go to 706
      fact(istart) = 1._R8
      denom = tpro(istart+1)-tpro(istart)
      if(denom.eq.0) denom = 1._R8
      facd(istart) = -1._R8/denom
      facd(istart+1) = 1._R8/denom
      go to 720
  706 continue
      do 700 l=1,ntpts-1
      if(tpro(l).le.time .and.                                           &  
     &  tpro(l+1).gt.time) go to 710
  700 continue
      fact(ntpts) = 1._R8
      go to 720
  710 continue
      if(icube .ge. 1) go to 719
!
!.....linear time point interpolation
      denom = tpro(l+1)-tpro(l)
      if(denom.eq.0) denom = 1._R8
      fact(l) = (tpro(l+1)-time)/denom
      fact(l+1) = (time-tpro(l))/denom
      facd(l) = -1._R8/denom
      facd(l+1) = 1._R8/denom
      go to 720
!
!     cubic time point interpolation
  719 denom1 = (tpro(l+1)-tpro(l))**3
      term1 = (time-tpro(l))**2*(3._R8*tpro(l+1)-tpro(l)-2._R8*time)/    &  
     & denom1
      fact(l+1) = term1
      fact(l) = 1._R8-term1
      term1d = 6._R8*(time-tpro(l+1))*(time-tpro(l))/denom1
      facd(l) = term1d
      facd(l+1) = -term1d
      if(l.le.1) go to 718
      denom2 = (tpro(l+1)-tpro(l))**2*(tpro(l+1)-tpro(l-1))
      term2 = (time-tpro(l+1))**2*(time-tpro(l))/denom2
      fact(l-1) = -term2
      fact(l+1) = fact(l+1) + term2
      term2d = (time-tpro(l+1))*(3._R8*time-2._R8*tpro(l)-tpro(l+1))/    &  
     & denom2
      facd(l-1) = -term2d
      facd(l+1) = facd(l+1) + term2d
      do ll=l-1,l+2
      enddo
  718 continue
      if(l+1 .ge. ntpts) go to 720
      denom3 = (tpro(l+1)-tpro(l))**2*(tpro(l+2)-tpro(l))
      term3 = (time-tpro(l))**2*(time-tpro(l+1))/denom3
      fact(l) = fact(l) - term3
      fact(l+2) = term3
      term3d = (time-tpro(l))*(3._R8*time-2._R8*tpro(l+1)-tpro(l))/      &  
     & denom3
      facd(l) = facd(l) - term3d
      facd(l+2) = term3d
      do ll=l,l+2
      enddo
!
  720 continue
      llsav = l
      do 721 l=1,ntpts
      if(fact(l) .le. 2) go to 721
      write(nout,1721) llsav,time,denom1,term1,denom2,term2,denom3,      &  
     & term3
      write(nout,1722) (tpro(li),li=1,ntpts)
      write(nout,1723) (fact(li),li=1,ntpts)
 1721 format(" error, llsave,time,denom,term",i3,1p7e12.4)
 1722 format(" tpro ",1p10e12.4)
 1723 format(" fact",1p10e12.4)
      ineg=22
  721 continue
!
!
!......redesigned 2/18/10    scj
      do l=1,ntpts
        factx(l) = fact(l)
      enddo
!
!.....for idata=9 or 11, call spdnstx and change fact array (but not factx)
      if(idata.eq.9 .or. idata.eq.11) call spdnstx
!
!
!......define coil current arrays
      do 595 ii=1,nwire
      ngr = iabs(igroupw(ii))
      cwire0(ii) = 0._R8
      if(kcycle.le.0) cwire0(ii) = cwics(ii)
      dcoil(ii) = 0._R8
      do 730 l=1,ntpts
      cwire0(ii)=cwire0(ii)+fact(l)*aturnsw(ii)*gcur(l,ngr)*usdi
!
      ngrvwl = ngrvw1(ii)
      if(ngrvwl.le.0) go to 829
      cwire0(ii) = cwire0(ii)+atnvw1(l,ii)*fact(l)*gcur(l,ngrvwl)*usdi
!
  829 continue
      ngrvwl = ngrvw2(ii)
      if(ngrvwl.le.0) go to 839
      cwire0(ii) = cwire0(ii)+atnvw2(l,ii)*fact(l)*gcur(l,ngrvwl)*usdi
!
  839 continue
      ngrvwl = ngrvw3(ii)
      if(ngrvwl.le.0) go to 849
      cwire0(ii) = cwire0(ii)+atnvw3(l,ii)*fact(l)*gcur(l,ngrvwl)*usdi
!
  849 continue
      ngrvwl = ngrvw4(ii)
      if(ngrvwl.le.0) go to 859
      cwire0(ii) = cwire0(ii)+atnvw4(l,ii)*fact(l)*gcur(l,ngrvwl)*usdi
!
  859 continue
  730 continue
!
!.....superimpose sinesoidal oscillation
      if(ngr.ne.int(acoef(80))) go to 594
      cwire0(ii) = cwire0(ii)+1000._R8*acoef(81)*usdi*aturnsw(ii)        &  
     &             *(sin(times*tpi/acoef(82)))
  594 continue
      if(ngr.ne.int(acoef(83))) go to 604
      cwire0(ii) = cwire0(ii)+1000._R8*acoef(84)*usdi*aturnsw(ii)        &  
     &             *(sin(times*tpi/acoef(82)))
  604 continue
      if(ngr.ne.int(acoef(85))) go to 614
      cwire0(ii) = cwire0(ii)+1000._R8*acoef(86)*usdi*aturnsw(ii)        &  
     &             *(sin(times*tpi/acoef(82)))
  614 continue
      if(ngr.ne.int(acoef(87))) go to 624
      cwire0(ii) = cwire0(ii)+1000._R8*acoef(88)*usdi*aturnsw(ii)        &  
     &             *(sin(times*tpi/acoef(82)))
  624 continue
  595 continue
!
      if(ncoil.eq.nwire) go to 598
      do 596 ii=1,ncoil-nwire
      ngr = iabs(igroupc(ii))
      ccoil0(ii) = 0._R8
      if(kcycle.le.0) ccoil0(ii) = ccics(ii)
      do 731 l=1,ntpts
      ccoil0(ii) = ccoil0(ii)+fact(l)*aturnsc(ii)*gcur(l,ngr)*usdi
!
      ngrvcl = ngrvc1(ii)
      if(ngrvcl.le.0) go to 929
      ccoil0(ii) = ccoil0(ii)+atnvc1(l,ii)*fact(l)*gcur(l,ngrvcl)*usdi
!
  929 continue
      ngrvcl = ngrvc2(ii)
      if(ngrvcl.le.0) go to 939
      ccoil0(ii) = ccoil0(ii)+atnvc2(l,ii)*fact(l)*gcur(l,ngrvcl)*usdi
!
  939 continue
      ngrvcl = ngrvc3(ii)
      if(ngrvcl.le.0) go to 949
      ccoil0(ii) = ccoil0(ii)+atnvc3(l,ii)*fact(l)*gcur(l,ngrvcl)*usdi
!
  949 continue
      ngrvcl = ngrvc4(ii)
      if(ngrvcl.le.0) go to 959
      ccoil0(ii) = ccoil0(ii)+atnvc4(l,ii)*fact(l)*gcur(l,ngrvcl)*usdi
!
  959 continue
  731 continue
!
!.....superimpose sinesoidal oscillation
      if(ngr.ne.int(acoef(80))) go to 597
      ccoil0(ii) = ccoil0(ii)+1000._R8*acoef(81)*usdi*aturnsc(ii)        &  
     &             *(sin(times*tpi/acoef(82)))
  597 continue
      if(ngr.ne.int(acoef(83))) go to 607
      ccoil0(ii) = ccoil0(ii)+1000._R8*acoef(84)*usdi*aturnsc(ii)        &  
     &             *(sin(times*tpi/acoef(82)))
  607 continue
      if(ngr.ne.int(acoef(85))) go to 617
      ccoil0(ii) = ccoil0(ii)+1000._R8*acoef(86)*usdi*aturnsc(ii)        &  
     &             *(sin(times*tpi/acoef(82)))
  617 continue
      if(ngr.ne.int(acoef(87))) go to 627
      ccoil0(ii) = ccoil0(ii)+1000._R8*acoef(88)*usdi*aturnsc(ii)        &  
     &             *(sin(times*tpi/acoef(82)))
  627 continue
  596 continue
!..rxw/18/12/87
!======================================================================
      if(acoef(290).ne.4._R8) go to 598
      do 10010 ig = 1,nreg
       sumturn = 0._R8
       sumcur  = 0._R8
       do 10020 ii = 1,ncoil-nwire
        if(igroupc(ii).ne.jnreg(ig)) goto 10020
        sumturn = sumturn + aturnsc(ii)
        sumcur  = sumcur  + ccoil0(ii)
10020  continue
!      preprogrammed external coil group current:
!@@@a
       if(sumturn .eq. 0) then
      write(nout,*)'***** efield: number of turns of ext. group=0 *****  &  
     & '
       endif
!
       ipow(ig)  = sumcur/sumturn * udsi
!
       if (nreg .ne. 7) then
      write(nout,*)'** external coils not usable in woyke controller **  &  
     & '
       endif
!@@@e
10010 continue
  598 continue
!
      apld0 = 0._R8
      gzero = 0._R8
      p0 = 0._R8
      r0 = 0._R8
      vloopp = 0._R8
      fbchix = 0._R8
      fbchiix = 0._R8
      xmagw=0._R8
      rzerw=0._R8
      azerw=0._R8
      ezerw=0._R8
      dzerw=0._R8
      zmagw=0._R8
      thalos = 0._R8
      heact = 0._R8
      whalos = 0._R8
      if(itevv .eq. 0) tevv = 0._R8
      if(iffac .eq. 0) ffac = 0._R8
      if(iimp .eq. 0) zeff = 0._R8
      alphar = 0._R8
      betar  = 0._R8
      fracpar = 0._R8
      qadd = 0.
      fhmodeix = 0._R8
      fhmodei = 0._R8
      pwidthc = 0.
      chiped = 0.
      tped = 0.
      nflag = 0.
      expn1 = 0.
      expn2 = 0.
      firitb = 0.
      secitb = 0.
      fracn0 = 0.
      newden = 0.
      do i=1,8
        fraci(i) = 0.
      enddo
!
      do 739 iigr=1,ngroupt
      igr=nogroupt(iigr)
  739 gvolt0(igr) = 0._R8
      do 740 l=1,ntpts
!
!.....total plasma current
      apld0 = apld0 + fact(l)*pcur(l)*usdi
!
!.....toroidal field constant gzero
      gzero = gzero + fact(l)*gzerov(l)
!
!.....plasma pressure
      p0 = p0 + factx(l)*ppres(l)*usdp
!
!.....normalized density
      r0 = r0 + factx(l)*rnorm(l)
!
!.....preprogrammed loop voltage
      vloopp = vloopp + factx(l)*vloopv(l)
!
!.....enhanced transport
      fbchix = fbchix + factx(l)*fbchia(l)
      fbchiix = fbchiix + factx(l)*fbchiia(l)
      xmagw = xmagw + factx(l)*xmagz(l)
      zmagw = zmagw + factx(l)*zmagz(l)
!
!.....other shape parameters
      rzerw = rzerw + factx(l)*rzerv(l)
      azerw = azerw + factx(l)*azerv(l)
      ezerw = ezerw + factx(l)*ezerv(l)
      dzerw = dzerw + factx(l)*dzerv(l)
!.....vacuum temperature
      if(itevv .eq. 0) tevv = tevv + factx(l)*tevv0(l)
!
!.....mass enhancement parameter
      if(iffac .eq. 0) ffac = ffac + factx(l)*ffac0(l)
!
!.....z-effective for resistivity calculation
      if(iimp .eq. 0) zeff = zeff + factx(l)*zeffv(l)
!
!.....exponents for density function
      alphar = alphar + factx(l)*alpharv(l)
      betar  = betar  + factx(l)*betarv (l)
!.....fraction parallel for current drive
      fracpar = fracpar + factx(l)*frcparv(l)
!.....halo parameters
      thalos = thalos + factx(l)*thalov(l)
      whalos = whalos + factx(l)*whalov(l)
!...helium ash confinement time
      heact  = heact  + factx(l)*heactv(l)
!
!     factor added to q in for itrmod=2  .. variable acoef(123)
      qadd = qadd + factx(l)*qaddv(l)
!
!     fhmode .... variable acoef(3003)
!      fhmodei = fhmodei + factx(l)*fhmodeiv(l)
!
!     (1 - pedestal width) .... variable acoef(3011)
!      pwidthc = pwidthc + factx(l)*pwidthcv(l)
!
!     chi_ped multiplier  .... variable acoef(3006)
!      chiped = chiped + factx(l)*chipedv(l)
!
!     pedestal temperature .... variable acoef(3102)
      tped = tped + factx(l)*tpedv(l)
!
!.....impurity fraction
      do i=1,8
        fraci(i) = fraci(i) + factx(l)*fraciv(i,l)
      enddo
!
!     first transport coefficient...variable acoef(3004)
!      firitb = firitb + factx(l)*firitbv(l)
!
!     second transport coefficient...variable acoef(3005)
!      secitb = secitb + factx(l)*secitbv(l)
!
!     edge density of n0 .... variable acoef(881)
      fracn0 = fracn0 + factx(l)*fracn0v(l)
!
!     new density parameter ... variable acoef(889)
      newden = newden + factx(l)*newdenv(l)
!
!.....preprogrammed group voltage
      do 741 iigr=1,ngroupt
      igr=nogroupt(iigr)
  741 gvolt0(igr) = gvolt0(igr) + fact(l)*gvolt(l,igr)
  740 continue
!
!......set those variable arrays that change abruptly
       do 400 l=1,ntpts-1
      lsav = l
      if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 410
  400 continue
  410 continue
      nflag = nflagv(lsav)
      expn1 = expn1v(lsav)
      expn2 = expn2v(lsav)
      firitb= firitbv(lsav)
      secitb= secitbv(lsav)
      chiped= chipedv(lsav)
      pwidthc=pwidthcv(lsav) 
      fhmodeix=fhmodeiv(lsav) 
      fhmodei=fhmodeiv(lsav) 
!
!.....make sure halo temp is above vacuum temp
      if(thalos.le.tevv) thalos = tevv
!
!.....modify preprogrammed voltages and define source terms
!.....for mission controller when acoef(296) = 5
      if(acoef(296) .eq. 5._R8) call missionc
      if(acoef(297) .eq. 1._R8)  call missionw
!
      e0 = acoef(2)*p0
!
!
!>>>NOTE:  The next 50 lines should be made into a seperate subroutine
!..feedback on chi multiplier to provide a global energy confinement
!..scaling
      fbmult = 1._R8
      if(acoef(3007) .gt. 0._R8) then
!.......line average density
        rlin = 0
        rdis = 0._R8
        do i=iminn,imaxx
          if(iexv(i,nh).eq.1 .or. iexs(i,nh).eq.1) cycle
          pval = psi(i,nh)
          if(pval.gt.psilim) cycle
          call reval(pval,idens,isurf,rval,rpval,i,nh)
          rlin = rlin + rval*deex*udsd
          rdis = rdis + deex
        enddo
        if(rdis.le.deex) rdis = deex
        zden = (rlin/rdis)*1.E-20_R8
        zr = rmajor
        za = rminor
        zip = tcurdtp*tpi*udsi*1.E-6_R8
        zk = shape5
        zd = shape3
!       zp = (pohmic+paux+prad+palpha)*1.e-6
        zp = (pohmic+paux+prad+palpha)*1.E-6_R8- wdotmw
        zb = gzero/zr
        zq = 2.5_R8*za**2*zb*(1._R8+zk**2)/(zip*zr)
        ztohmic = 0.07_R8*zden*zr**2*za*zq
        ztau = ztohmic
        if(zp.gt.0) then
          taueka = .172*zk**.78*zip**0.93*zden**.41*zb**.15             &
     &             *za**.58*zr**1.39/zp**.69
!          taueka = .0578_R8*zk**.64_R8*zip**0.96_R8*zden**.40_R8*zb**.03_R8      &
!                             *zr**1.89_R8*(2.5_R8)**0.2_R8/(zp**.73_R8*za**0.06_R8)
          ztau = 1._R8/sqrt((1._R8/taueka**2)+(1._R8/ztohmic**2))
        endif
        fbmult = enerst2/((hfluxav-wdotmw*1.E6_R8)*acoef(3007)*taueka)
!       fbchiax = fbchia(l)*fbmult
!       write(nterm,*) zden,zr,za,zip,zk,zd,zp,zb,zq
!       write(nterm,*) acoef(3007),taueka,enerst2,hfluxav
!       write(nterm,*) fbmult,fbchiax
      endif
!
      fbchi = 1._R8
      fbchii = 1._R8
      if(fbchix .gt. 0) fbchi = (fbchix*fbmult*acoef(3008)               &  
     & +(1._R8-acoef(3008))*fbchixold)
      fbchixold=fbchi
      if(fbchiix .gt. 0) fbchii = (fbchiix*fbmult*acoef(3008)            &  
     & +(1._R8-acoef(3008))*fbchiixold)
      fbchiixold=fbchii
!
!------------------------------------------------------------------
!	pedestal temperature feedback
!	acoef(3208)=1 uses EPED1 grid of values (ITER only) 
!	acoef(3208)=2 Sugihara
!	acoef(3208)=3 Maget model (ITER only)
!	acoef(3208)>10 uses this as target pedestal 
      tped_corr = 1._R8
      if (fhmodeixold .le. 0.) fhmodeixold = fhmodeix
      if (acoef(3208) .gt. 0.) then 
         call tped_fb(acoef(3208),tped_corr)
         if (fhmodeix .gt. 0.) fhmodei= fhmodeixold*tped_corr
      endif
      fhmodeixold = fhmodei
!------------------------------------------------------------------
!
!......try using plasma current for burn control
      excess = (bpowers+pohmic+palpha) - acoef(47)*1.E6_R8
      if(excess .lt. 0) excess = 0._R8
      apld0 = apld0 - acoef(54)*excess*usdi
!
!
!.....rescale units and reset vacuum temperature
      call rscale
      if(acoef(290).ne.4._R8) go to 10013
!
!**** warning: coi can be used only for woyke-controller
!              if modelled by single coils and group 8,9 !! *****
      icnt8 = 0
      icnt9 = 0
      if (ncoil .eq. nwire) goto 10012
      do 10011 iw = 1,nwire
      if (igroupw(iw) .eq. 8) then
      icnt8 = iw
      ipow (8)  = cwire0 (icnt8) / coiturn * udsi
      endif
      if (igroupw(iw) .eq. 9) then
      icnt9 = iw
      ipow (9)  = cwire0 (icnt9) / coiturn * udsi
      endif
10011 continue
10012 continue
       if (icnt8 .eq. 0 .or. icnt9 .eq. 0) then
      write(nout,*)'** efield: coi not usable in woyke controller **'
       endif
!     sollwerte fuer gleichgewichtsparameter und plasmastrom
      icgso(1) = ipow (8) + ipow (9)
      icgso(2) = ipow (8) - ipow (9)
      gsoll(1) = rzerw
      gsoll(2) = azerw
      gsoll(3) = ezerw
      gsoll(4) = dzerw
      ipsoll   = apld0 * udsi
 
!    skalierungsfaktor fuer kontrollspuleninduktivitaet
      skal = 22.0E-06_R8/ (rswire (1) * udsv / udsi)
 
!   ableitungen der fuehrungsgroessen fuer tracking
!   a) initialisierung bei neustart und restart
        if ((kcycle .le. 1) .or. (iswtch.eq.0)) then
           iswtch = 1
           do 101 i = 1,6
              gsollm1 (i) = gsoll (i)
              gsollp  (i) = 0.0_R8
101        continue
           ipsollm1 = ipsoll
           ipsollp  = 0.0_R8
 
!    b) oder aktualisierung
        else
           do 102 i = 1,6
              gsollp  (i) = (gsoll (i) - gsollm1 (i))/(dt*udst)
              gsollm1 (i) = gsoll (i)
102        continue
           ipsollp  = (ipsoll - ipsollm1)/(dt*udst)
           ipsollm1 = ipsoll
        endif
10013 continue
!@@@e
!
!.......................................................................
!.....part 2: add feedback current to plasma and coil currents
!.......................................................................
!
!     if(acoef(511).le.0 .or. acoef(510).eq.0) go to 7778
!  section added to generate coefficients required in random
!  feedback signal  CK  1/93
!     if(ckstrt .gt. 0.0) go to 7778
!     ckstrt=1.0
!     call g05cce
!     do 7777 kk=1,ifix(acoef(510)/acoef(511))
!     wmega(kk)=float(kk)*2.0*pi/acoef(510)
!     phin(kk)=g05dae(-1.0,1.0)
!     phin(kk)=phin(kk)*pi
!777  continue
!778  continue
      if(idata.eq.7) call tcv2
      if(numfb.eq.0) go to 811
      dtramp = acoef(3)*usdt
!.....accumulate feedback group current in aif5b and gcurfb
!
      aplfb = 0._R8
      if(acoef(901).gt.0._R8) goto 8050
      do 805 ig=1,pngroup
  805 gcurfb(ig) = 0._R8
8050  continue
      do 803 ig=1,pngroup
      do 804 n=1,pncoil
      fadj(ig,n)=1.0_R8
804   continue
803   continue
      do 810 l=1,numfb
!
      ig = nrfb(l)
      nl = l
      ie = ipext(l)
      if(ie.gt.1000) go to 1816
      if(ie .eq. 25) go to 1825
      if(ie.gt.20) go to 1820
      if(ie.gt.10) go to 1910
      go to(806,806,806,1802,1810,1811,1812,1813,1814,1815),ie
  806 continue
      fluxdif = 0._R8
      do 808 ll=1,ntpts
      if(factx(ll).eq.0._R8) go to 808
      n = 2*nfeedv(ll,l)-1
      if(nrecord.eq.0 .or. ie.eq.1) go to 807
      fluxdif = fluxdif + factx(ll)*(fluxu(n,nrecord)                    &  
     &                             -fluxl(n,nrecord))
      go to 808
  807 np = n+1
      psinn = pinterp(xobs(n ),zobs(n ),iobs(n ),jobs(n ))
      psinp = pinterp(xobs(np),zobs(np),iobs(np),jobs(np))
      fluxdif = fluxdif + factx(ll)*(psinn-psinp)
  808 continue
      go to 1803
!
!.....sensing signal is plasma current for ipext(l)=4
 1802 continue
      fluxdif = (tcurdtp*tpi - apld0)*udsi
      go to 1803
 1810 continue
      fluxdif=xmag-xmagw
      go to 1803
 1811 continue
!     sumsig=0.0
!     do 7779 kk=1,ifix(acoef(510)/acoef(511))
!     sumsig=sumsig+exp(-(wmega(kk)**2*acoef(511)**2/8.0))
!    +*cos(wmega(kk)*times-phin(kk))
!779  continue
!     fzmag=2.0*acoef(512)*sqrt(sqrt(pi)*acoef(511)/acoef(510))
!    +*sumsig
!     zmagr=zmagw*fzmag
!     fluxdif=zmag-zmagr
      fluxdif=zmag-zmagw
      go to 1803
 1812 fluxdif = eps1c-eps10
      go to 1803
 1813 fluxdif = eps2c-eps20
      go to 1803
 1814 fluxdif = eps3c-eps30
      go to 1803
 1815 fluxdif = eps4c-eps40
      go to 1803
 1910 continue
      call shape(ellip,delta,x1,x2,z1,z2)
      iem10 = ie - 10
      go to(1911,1912,1913,1914),iem10
 1911 fluxdif = .5_R8*(x2+x1) - rzerw
      go to 1803
 1912 fluxdif = .5_R8*(x2-x1) - azerw
      go to 1803
 1913 fluxdif = ellip - ezerw
      go to 1803
 1914 fluxdif = delta - dzerw
      go to 1803
 1820 continue
      call fluxmod(psi1,psi2,psi3,psi4,psisn1,2.40_R8,0.40_R8)
      iem20 = ie - 20
      go to(1821,1822,1823,1824),iem20
 1821 fluxdif = psi1
      go to 1803
 1822 fluxdif = psi2
      go to 1803
 1823 fluxdif = psi3
      go to 1803
 1824 fluxdif = psi4
      go to 1803
 1816 continue
      index = ie-1000 + ncoil-nwire
      fluxdif = ccoil(index)*udsi
 1825 continue
      alphac=acoef(510)*sqrt((psilim-psimin)/(apl)**2)
      call shpcntl(ig,nl,alphac,term)
      go to 1888
 1803 continue
!
      if(idelay(l) .le. 0) go to 1850
 1801 continue
      indxd1(l) = indxd1(l) + 1
      if(indxd1(l) .gt. idelay(l)) indxd1(l) = indxd1(l) - idelay(l)
      fluxdif1 = 0._R8
      if(kcycle.gt.idelay(l)) fluxdif1 = fsave1(l,indxd1(l))
      fsave1(l,indxd1(l)) = fluxdif
      fluxdif = fluxdif1
 1850 continue
!
      if(tfbon(l).ge.0._R8 .and. time.lt.tfbon(l)) go to 809
      if(tfbons(l).lt.0._R8 .and. kcycle .lt. int(-tfbons(l)) )          &  
     &        go to 809
      if(tfbof(l).ge.0._R8 .and. time.gt.tfbof(l)) go to 809
      if(tfbofs(l).lt.0._R8 .and. kcycle .gt. int(-tfbofs(l)) )          &  
     &       go to 809
!
!.....define linear ramp form-factor for turning on system slowly
      formf = 1._R8
      if(tfbon(l).ge.0._R8 .and. dtramp.gt.0 .and. time.lt.              &  
     &                                     tfbon(l)+dtramp)              &  
     &formf = (time-tfbon(l))/dtramp
!
      fac = 1._R8
      if(fbfac1(l).gt.0) fac = fac*apld0*udsi/pcur(ntpts)
      if(kcycle.le.0.and.irst1.ne.2) savfeed(l) = fluxdif - fbcon(l)*    &  
     & fac
!
!@@@a rxw/16/1/89
!     multiply feedback constants by the preprogrammed plasma
!     current normalized to 1 ma, if feedback on xmag or zmag
      facrxw = 1._R8
      if( (ipext(l).eq.5) .or. (ipext(l).eq.6) )                         &  
     &  facrxw = apld0*udsi*1.E-6_R8
!
      if(kcycle.le.0.and.irst1.ne.2)                                     &  
     &    savfeed(l) = (fluxdif - fbcon(l)*fac) * facrxw
!
!.....average flux signal over two cycles
!
      fluxoff = (fluxdif-fbcon(l)*fac ) * facrxw
!@@@e
      term = fbfac(l)*0.5_R8*(fluxoff+savfeed(l))
!
      sumfeed(l) = sumfeed(l) + fluxoff*dts
      term = term + fbfaci(l)*sumfeed(l)
!
      if((kcycle.eq.0.and.irst1.ne.2).or.dts.eq.0) go to 1804
      termd = (fluxoff - savfeed(l))/dts
      term = term + fbfacd(l)*termd*formf
 1804 savfeed(l) = fluxoff
!
      if(nrfb(l).eq.0) go to 820
      ig = nrfb(l)
!
 1888 continue
      gcurfb(ig) = gcurfb(ig)+term*usdi
      go to 809
  820 aplfb = aplfb + term*usdi
  809 continue
!>>>>>debug
!     write(6,1809) l,ig,ie,fluxdif,term,gcurfb(ig)
!     write(6,1808) tfbon(l),tfbof(l),fbfac1(l),fluxoff,savfeed(l),sumfeed(l),termd,facrxw
!1808 format(12x,1p10e12.4)
!1809 format(" from efield",3i5,1p3e12.4)
  810 continue
  811 continue
!
!
!....add feedback plasma current to preprogrammed value
      apld0 = apld0 + aplfb
!
!....multiply each feedback group current by no of turns
!....and add to preprogrammed current
!
!....note factx multiplies no of turns, but fact multiplies currents
      do 830 ii=1,nwire
      iii=ii+(ncoil-nwire)
      ngr = iabs(igroupw(ii))
      cwire0(ii) = cwire0(ii)+fadj(ngr,iii)*aturnsw(ii)*gcurfb(ngr)
!
      ngrvwl = ngrvw1(ii)
      if(ngrvwl.le.0) go to 840
      atuv = 0.0_R8
      ppcur = 0._R8
      do 835 l=1,ntpts
      ppcur = ppcur + fact(l)*gcur(l,ngrvwl)*usdi
  835 atuv = atuv + atnvw1(l,ii)*factx(l)
      cwire0(ii) = cwire0(ii)                                            &  
     &     +atuv*(fadj(ngrvwl,iii)*gcurfb(ngrvwl)+ppcur)
!
  840 continue
      ngrvwl = ngrvw2(ii)
      if(ngrvwl.le.0) go to 850
      atuv = 0.0_R8
      ppcur = 0._R8
      do 845 l=1,ntpts
      ppcur = ppcur + fact(l)*gcur(l,ngrvwl)*usdi
  845 atuv = atuv + atnvw2(l,ii)*factx(l)
      cwire0(ii) = cwire0(ii)                                            &  
     &     +atuv*(fadj(ngrvwl,iii)*gcurfb(ngrvwl)+ppcur)
!
  850 continue
      ngrvwl = ngrvw3(ii)
      if(ngrvwl.le.0) go to 860
      atuv = 0.0_R8
      ppcur = 0._R8
      do 855 l=1,ntpts
      ppcur = ppcur + fact(l)*gcur(l,ngrvwl)*usdi
  855 atuv = atuv + atnvw3(l,ii)*factx(l)
      cwire0(ii) = cwire0(ii)                                            &  
     &     +atuv*(fadj(ngrvwl,iii)*gcurfb(ngrvwl)+ppcur)
!
  860 continue
      ngrvwl = ngrvw4(ii)
      if(ngrvwl.le.0) go to 870
      atuv = 0.0_R8
      ppcur = 0._R8
      do 865 l=1,ntpts
      ppcur = ppcur + fact(l)*gcur(l,ngrvwl)*usdi
  865 atuv = atuv + atnvw4(l,ii)*factx(l)
      cwire0(ii) = cwire0(ii)                                            &  
     &     +atuv*(fadj(ngrvwl,iii)*gcurfb(ngrvwl)+ppcur)
!
  870 continue
  830 continue
!
      if(ncoil.eq.nwire) go to 890
      do 996 ii=1,ncoil-nwire
      ngr = iabs(igroupc(ii))
      ccoil0(ii) = ccoil0(ii) + fadj(ngr,ii)*aturnsc(ii)*gcurfb(ngr)
!
      ngrvcl = ngrvc1(ii)
      if(ngrvcl.le.0) go to 940
      atuv = 0.0_R8
      ppcur = 0._R8
      do 935 l=1,ntpts
      ppcur = ppcur + fact(l)*gcur(l,ngrvcl)*usdi
  935 atuv = atuv + atnvc1(l,ii)*factx(l)
      ccoil0(ii) = ccoil0(ii)                                            &  
     &     +atuv*(fadj(ngrvcl,ii)*gcurfb(ngrvcl)+ppcur)
!
  940 continue
      ngrvcl = ngrvc2(ii)
      if(ngrvcl.le.0) go to 950
      atuv = 0.0_R8
      ppcur = 0._R8
      do 945 l=1,ntpts
      ppcur = ppcur + fact(l)*gcur(l,ngrvcl)*usdi
  945 atuv = atuv + atnvc2(l,ii)*factx(l)
      ccoil0(ii) = ccoil0(ii)                                            &  
     &     +atuv*(fadj(ngrvcl,ii)*gcurfb(ngrvcl)+ppcur)
!
  950 continue
      ngrvcl = ngrvc3(ii)
      if(ngrvcl.le.0) go to 960
      atuv = 0.0_R8
      ppcur = 0._R8
      do 955 l=1,ntpts
      ppcur = ppcur + fact(l)*gcur(l,ngrvcl)*usdi
  955 atuv = atuv + atnvc3(l,ii)*factx(l)
      ccoil0(ii) = ccoil0(ii)                                            &  
     &     +atuv*(fadj(ngrvcl,ii)*gcurfb(ngrvcl)+ppcur)
!
  960 continue
      ngrvcl = ngrvc4(ii)
      if(ngrvcl.le.0) go to 970
      atuv = 0.0_R8
      ppcur = 0._R8
      do 965 l=1,ntpts
      ppcur = ppcur + fact(l)*gcur(l,ngrvcl)*usdi
  965 atuv = atuv + atnvc4(l,ii)*factx(l)
      ccoil0(ii) = ccoil0(ii)                                            &  
     &     +atuv*(fadj(ngrvcl,ii)*gcurfb(ngrvcl)+ppcur)
  970 continue
!
  996 continue
  890 continue
!
!.....solve circuit equations to define ccoil
      call cirequ
!
!.......................................................................
!.....part 3: compute oh electric field to keep plasma current at preprescribed
!.......................................................................
!
      aplav = .5_R8*(apl+aplo)
      plcur = tcurdtp*tpi
      rebouno = rebounn
      fac = sqrt(nx*(nz-1)/((2-isym)*4._R8*(alx-ccon)*alz))
!.....turn off loop voltage feedback during disruption
      if(times .gt. acoef(95)) fac = 0._R8
      rebounn = (-etav/ndiv*fac*(plcur-apld0)*acoef(11)                  &  
     &        + vloopp*usdv/tpi)
      reboun = .5_R8*(rebouno+rebounn)
!
      if(kcycle .le. 3 .and. irst1.ne.2) reboun = vloopp*usdv/tpi
      if(lrswtch.ne.0) reboun = vloopp*usdv/tpi
      rebmax = 1.E20_R8
      rebmin = -1.E20_R8
      if(acoef(15) .ne. 0) rebmin = acoef(15)*usdv/tpi
      if(acoef(16) .ne. 0) rebmax = acoef(16)*usdv/tpi
      if(reboun .gt. rebmax) reboun = rebmax
      if(reboun .lt. rebmin) reboun = rebmin
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
