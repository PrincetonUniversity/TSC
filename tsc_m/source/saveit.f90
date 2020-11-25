!#include "f77_dcomplx.h"
      subroutine saveit
!......6.91 saveit
!
!.....save data for time history plots
!
      USE CLINAM
      USE SCR14 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.....note that nlab is dimensioned pglob and must be defined
!.....in a data statement in this subroutine
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,k,mwpl2,idisk,l,np,n,iimax,ii,ic,nopl1,nopl2,n10
      INTEGER icc,nmax,indx,ip,iip,iimpgls,llab,iig,ig,il,iabs
      INTEGER izero,i1,i2,i3,iadscr,irec
      INTEGER len, index, if, maxdakoda, ios
      INTEGER, allocatable :: idakoda(:)
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xlist,ylist,xlist2,ylist2,ss,psinn,pinterp,psinp
      REAL*8 tmaxx,tmin,ymin,ymax,tcut,ymax2,diffy,xptmin
      REAL*8 riendsv,ymaxrie,yminrie,rendsv,rieendsv,tplot,dtplot
      REAL*8 xmin,xmax,diffx,ristrsv,rstrsv,smin,smax,aindx,xl,xr
      REAL*8 zl,zr,dis,disncm,ymaxp,gapmin,gapmax,fbmin,fbmax
      REAL*8 fbdiff,amax,amin,sum2
      REAL*8 AREAL, sum
      REAL*8, allocatable :: adakoda(:)
!============
      character*1 lab2(35)
!     CHARACTER*8, DIMENSION(pglobs) :: nlab
      CHARACTER*10, allocatable :: ndakoda(:)
!     data (nlab(i),i=1,40)/
!    1'  amach ','   ekin ','   ekinp','   ekint','  iplim ',
!    2'  zmag  ','  xmag  ','   cur  ','   cur  ','delp*tpi',
!    3'pmin*tpi',' diamag ','surfvolt','  qzero ','  qedge ',
!    3'    dt  ','beta-tor','<n>/nmur','1/q nr/b','li/2+bp ',
!    3'betapol ',' li/2   ',' li vs q','taue-kg ','taue(ms)',
!    3' taue2  ',' te(0)  ',' ti(0)  ','te/te-av','ti/ti-av',
!    3' chiauxs','chiohms ','hflux-mw','   r0   ','minorrad',
!    3'delt-tri','ellip   ',' xsep(1)',' zsep(1)','vsec-int'/
!     data (nlab(k),k=41,96)/
!    3'resv-sec','vsec-abs','loopv-oh','vsec-oh ','ptot(mw)',
!    3'palpha  ','pohmic  ','paux    ','nullapol','nullapol',
!    3'  dipole','  dipole','quadrupl','quadrupl','hexapole',
!    3'hexapole','octapole','octapole','decapole','decapole',
!    4'  q95   ','  qstar ','  qcyl  ',' delta95',' ellip95',
!    5' xsepcal',' zsepcal','  psep  ',' psepcal',' n/ncrit',
!    6' dt1s   ',' dt2s   ',' dt3s   ','1/qcylin',' beta-b0',
!    7' beta-0 ','n-lineav','n-vol av','density ','int ener',
!    8'ellip-90',' delt-90',' volume ',' no part','  eps1a ',
!    9'  eps1c ','  eps10 ','  eps2a ','  eps2c ','  eps20 ',
!    .'  eps3a ','  eps3c ','  eps30 ','  eps4a ','  eps4c ',
!    1'  eps40 '/
!     data (nlab(k),k=97,109)/
!    1'thm ener','pf ener ','tf ener ','hflux en','aux ener',
!    2'int c en','pf bn en','tf bn en','total en','  tevv  ',
!    3'  ffac  ',' npsit  ',' resid  '/
!
!     data (nlab(k),k=125,126)/'cpu(min)','cpu-rate'/
!     data (nlab(k),k=129,130)/'alphdens','alphbeta'/
!     data (nlab(k),k=134,139)/'alphfrac',' zeff   ',
!    1                         ' ni-ave ',' nhe-ave',
!    2                         ' W-dot  ','  H(98) '/
!     data (nlab(k),k=147,148)/' xsep(2)',' zsep(2)'/
!     data (nlab(k),k=151,152) / ' r(q=1) ',' li-(ga)'/
!     data (nlab(k),k=164,165) / ' thalo  ',' whalo  '/
!     data (nlab(k),k=166,167) / ' Pohm-H ',' Pohm-P '/
!     data (nlab(k),k=168,170) / 'VS-Poynt','        ','Ejima Co'/
!     data (nlab(k),k=173,174) /'   li   ','   Ct   '/
!     data (nlab(k),k=175,177) / 'Ptot    ',' Palpha ','frac-T  '/
!     data (nlab(k),k=200,204) / 'dndtpel ','xpel    ',
!    1 'tepel   ','totradp ','Radpel  '/
!     data (nlab(k),k=205,206) /'pellatom','tot atom'/
!     data (nlab(k),k=207,212) / 'iqtrubmx','nskipsf ',
!    1 'runawayc','plas_cur','avalanch','resource' /
!     data (nlab(k),k=218,220) / 'Te(ped) ','Ti(ped) ','Inj Curr'/
!     data (nlab(k),k=221,221) / 'Injcur 2'/
!     data (nlab(k),k=222,222) / '-dW core'/
!     data (nlab(k),k=224,224) / '-dW tot '/
!
      data lab2/"1","2","3","4","5","6","7","8","9",                     &  
     &          "a","b","c","d","e","f","g","h","i","j",                 &  
     &          "k","l","m","n","o","p","q","r","s","t",                 &  
     &          "u","v","w","x","y","z"/
!
!
!
      dimension xlist(20),ylist(20),xlist2(20),ylist2(20)
      data xlist2/                                                       &  
     &  8.0_R8, 8.0_R8, 1.0_R8, 1.0_R8, 8.0_R8,                          &  
     &  6.0_R8, 5.0_R8, 4.0_R8, 3.0_R8, 2.0_R8,                          &  
     &  2.0_R8, 3.0_R8, 3.0_R8, 4.0_R8, 4.0_R8,                          &  
     &  5.0_R8, 5.0_R8, 6.0_R8, 6.0_R8, 8.0_R8/
      data ylist2/                                                       &  
     &  .46_R8, 0.1_R8, 0.1_R8, 1.1_R8, 1.1_R8,                          &  
     &  1.0_R8, 0.9_R8, 0.8_R8, 0.7_R8, .54_R8,                          &  
     &  .34_R8, 0.5_R8, .35_R8, .48_R8, .35_R8,                          &  
     &  .45_R8, .42_R8, .44_R8, .43_R8, .46_R8/
      if(ifrst(5) .ne. 0) go to 50
!
      if(irfp.eq.1) then
      nlab(12) = 'helicity'
      nlab(23) = 'tflux   '
      nlab(36) = ' psi*1th'
      nlab(37) = ' psi*2th'
      nlab(38) = ' psi*3th'
      nlab(39) = ' psi*4th'
      nlab(80) = ' nloop  '
      endif
      if(acoef(504).ne.0) nlab(12) = 'helicity'
      if(acoef(504).ne.0) nlab(80) = 'mag ener'
      if(isym.eq.1 .and. acoef(61).eq.0 .and. lrswtch.gt.0)              &
     & nlab(5) = 'Ivv(pol)'
      if(isym.eq.1 .and. acoef(61).eq.0) nlab(6) = 'R * Btor'
      if(acoef(79) .eq. 1.) nlab(7) = 'Gapvolt'
!
!
!.....open scratch file for global time history plots
      ss = 0
      mwpl2 = 14
      idisk = 1
      len = ((ncycle)/(nskip2+1)+2)*(pglobp)
      ifrst(5) = 1
   50 continue
!
      do 60 l=2,pgls+1+3*pncoil+10*pngroup+6
   60 temp(l) = global(l-1)
      temp(1) = time*udst
!
!.....write global variables into pltsav array
!
      nrecord = nrecord + 1
      if(nrecord.gt.2*pnsave) then
        ineg=63
        write(*,*) "nrecord exceeds pnsave", nrecord
        return
      endif
      lenscr = pgls+1+3*pncoil+10*pngroup+6
      do 12 l=1,lenscr
      index = (nrecord-1)*lenscr+l
      if(index.gt.2*pnsave*pglobp) then
        ineg=63
        write(*,*) "index exceeds array bounds", index
        return
      endif
   12 pltsav(index) = temp(l)
!
!.....save current trajectory
      apld0sv(nrecord) = apld0*udsi
!
!.....save quantities needed by iter
      if(acoef(1) .eq. 5._R8)  call iterbp(1)
!
!.....save flux observation data
!
      do 80 nn=1,nobs-1,2
      np = nn+1
      psinn = pinterp(xobs(nn),zobs(nn),iobs(nn),jobs(nn))
      psinp = pinterp(xobs(np),zobs(np),iobs(np),jobs(np))
      fluxu(nn,nrecord) = psinn
   80 fluxl(nn,nrecord) = psinp
!
!.....save data for ufiles for acoef(3001).eq.1.0
      if(acoef(3001).eq.1) call ufstore
!
!.....save initial values
      if(kcycle.gt.0) return
      do 83 n=1,nobs-1,2
      fluxu0(n) = fluxu(n,nrecord)
      fluxl0(n) = fluxl(n,nrecord)
   83 continue
!
!
      return
!
      entry plotit2
!
!.....produce global time history plots
!
      if(ifrst(5) .eq. 0) return
!
!.....call Neil's subroutine to dump plot data file
      call dumpdat
      tmaxx = time*udst
      tmin = pltsav(1)
      if(tmaxx .le. tmin .or. nrecord.le.1) return
!
      lenscr = pgls+1+3*pncoil+10*pngroup+6
!
!.....note that wire and group plots deferred until later
      iimax = pgls+ncoil
      do 438 i=1,pgls
      if(globmax(i).eq.-1.E30_R8) globmax(i) = 0._R8
      if(globmin(i).eq. 1.E30_R8) globmin(i) = 0._R8
  438 continue
      do 400 ii=1,iimax
!
!
!.....first determine ymin and ymax and store in big1
!
      if(ii.le.pgls) ymin = globmin(ii)
      if(ii.le.pgls) ymax = globmax(ii)
!     KSTAR
!     ymin = globmin(ii)
!     ymax = globmax(ii)
      if(ii.eq.2) ymin = 0._R8
      if(ii.eq.6 .and. isym.eq.1 .and. acoef(61).eq.0)                   &  
     &                 ymax = max(globmax(6),globmax(150))
      if(ii.eq.6 .and. isym.eq.1 .and. acoef(61).eq.0)                   &  
     &                 ymin = min(globmin(6),globmin(150))
      if(ii.eq.8) ymax = max(globmax(8),globmax(9),apld0sv(nrecord),     &
     &                         globmax(235))
      if(ii.eq.8) ymin = min(globmin(8),globmin(9),apld0sv(1),           &
     &                         globmin(235))
      if(ii.eq.10) ymax = max(globmax(10),globmax(234),globmax(236))
      if(ii.eq.10) ymin = min(globmin(10),globmin(234),globmin(236),0.)
!
      if(ii.eq.14.and.isurf.ne.0.and.globmax(171).gt.globmin(171)) then
!     KSTAR
!     ymin = 0.0
!     ymax = 5.0
      ymax = max(globmax(14),globmax(171),globmax(172))
      ymin = min(globmin(14),globmin(171),globmin(172))
      endif
      if(ii.eq.15.and.irfp.ne.1) ymax = max(globmax(15),globmax(61),     &  
     &                          globmax(62),globmax(63))
      if(ii.eq.15.and.irfp.ne.1) ymin = min(globmin(15),globmin(61),     &  
     &                          globmin(62),globmin(63))
      if(ii.eq.15 .and. irfp.ne.1) ymax = min(10._R8,ymax)
      if(ii.eq.15 .and. irfp.ne.1) ymin = max(1._R8,ymin )
      if(ii.eq.16) ymax = max(globmax(16),globmax(71),                   &  
     &                          globmax(72),globmax(73))
      if(ii.eq.16) ymin = min(globmin(16),globmin(71),                   &  
     &                          globmin(72),globmin(73))
      if(ii.eq.17) ymin = min(globmin(17),globmin(75),globmin(76))
      if(ii.eq.17) ymax = max(globmax(17),globmax(75),globmax(76))
      if(ii.eq.18) ymax = max(globmax(18),globmax(70),globmax(149))
      if(ii.eq.18) ymin = min(globmin(18),globmin(70),globmin(149))
      if(ii.eq.20) ymin = min(globmin(20),globmin(21),globmin(22))
      if(ii.eq.20) then
          call disk(ii+1,1,big1,big2)
          tcut = tmin + 0.1_R8*(tmaxx-tmin)
          ymax2 = -1.E60_R8
          do 8011 i=1,nrecord
          if(big2(i).lt.tcut) go to 8011
          ymax2 = max(ymax2,1.5_R8*big1(i))
 8011     continue
          ymax = min(ymax,ymax2)
                   endif
      if(ii.eq.24    .and. irfp.ne.1)                                    &  
     & ymax = max(globmax(24),globmax(131),globmax(132),globmax(133))
      do 452 ic=25,31,2
      if(ii.eq.ic.or.ii.eq.ic+1) ymin=min(globmin(ic),globmin(ic+1))
      if(ii.eq.ic.or.ii.eq.ic+1) ymax=max(globmax(ic),globmax(ic+1))
  452 continue
      if(irfp.eq.1) go to 453
      if(ii.eq.36) ymax = max(globmax(36),globmax(64),globmax(82))
      if(ii.eq.36) ymin = min(globmin(36),globmin(64),globmin(82))
      if(ii.eq.37) ymax = max(globmax(37),globmax(65),globmax(81))
      if(ii.eq.37) ymin = min(globmin(37),globmin(65),globmin(81))
      if(ii.eq.38 .and. isvd.eq.1) ymax = max(globmax(38),globmax(66))
      if(ii.eq.38 .and. isvd.eq.1) ymin = min(globmin(38),globmin(66))
      if(ii.eq.39 .and. isvd.eq.1) ymax = max(globmax(39),globmax(67))
      if(ii.eq.39 .and. isvd.eq.1) ymin = min(globmin(39),globmin(67))
  453 continue
      if(ii.eq.40) ymin = min(globmin(40),globmin(41))
      if(ii.eq.40) ymax = max(globmax(40),globmax(41),                   &  
     &                          globmax(42)-globmin(42)+globmax(44))
      if(ii.eq.45) ymax = max(globmax(45),globmax(46),acoef(706))
      if(ii.eq.45) ymin = min(globmin(119),globmin(120),globmin(121)     &  
     &                        ,acoef(705))
      if(ii.lt.49 .or. ii.gt.59) go to 454
      ymin = min(globmin(ii),globmin(ii+1))
      ymax = max(globmax(ii),globmax(ii+1))
  454 continue
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy .le. 0) ymax = ymin + max(1._R8,abs(ymin))
!
!.....store variable to be plotted in big1  .. time is in big2
      call disk(ii+1,1,big1,big2)
      tmin = big2(1)
!
!.....store time in big3
      do 455 i=1,nrecord
  455 big3(i)=big2(i)
!
!.....special logic for coils only option
      if(lrswtch.ge.1) then
        if(ii.eq.48) go to 903
        if(ii.eq.8) then
          ymax = globmax(8)
          ymin = globmin(8)
          go to 20
                    endif
        if(ii.eq.16) then
          ymax = globmax(16)
          ymin = globmin(16)
          go to 20
                    endif
      if(ii.eq.5 .and. isym.eq.1 .and. acoef(61).eq.0) go to 10
      if(ii.eq.6 .and. isym.eq.1 .and. acoef(61).eq.0) go to 20
      if(ii.eq.12.and. isym.eq.1 .and. acoef(61).eq.0) go to 20
        go to 400
                       endif
!
!.....special logic for RFP option
      if(irfp.eq.1) then
        if(ii.eq.24) go to 810
        if(ii.eq.25) go to 820
        if(ii.eq.18) go to 400
        if(ii.eq.19) go to 782
                    endif
!
!....special logic for non-surface averaged option
      if(isurf.eq.0) then
        if(ii.eq.18) go to 400
        if(ii.eq.19) go to 782
        if(ii.ge.24.and.ii.le.33) go to 400
        if(ii.ge.45.and.ii.lt.48) go to 400
        if(ii.eq.48) go to 903
                     endif
!
!.....special additional logic for CHI modeling
      if(acoef(504).ne.0) then
        if(ii.eq.5  .or. ii.eq.6)  go to 400
        if(ii.eq.14 .or. ii.eq.15) go to 400
!       if(ii.ge.20 .and. ii.le.23)go to 400
        if(ii.eq.34 .or. ii.eq.35) go to 400
        if(ii.eq.36 .or. ii.eq.37) go to 400
                          endif
!
!......check if divertor option is not used
      if(idiv.eq.0 .and. irfp.eq.0) then
        if(ii.eq.38) go to 400
        if(ii.eq.39) go to 400
                    endif
!
!.....check if automatic OH system is turned off
      if(acoef(11).eq.0) then
        if(ii.eq.43) go to 400
        if(ii.eq.44) go to 400
                          endif
!
!.....normal case
      if(ii.ge.61          .and. ii.le.pglobs)      go to 400
      if(ii.gt.pglobs      .and. ii.le.pglobs+ppsi) go to 950
      if(ii.gt.pglobs+ppsi .and. ii.le.pgls)        go to 400
      if(ii.gt.pgls)                                go to 500
!
!                           note: 1. 950 processes balloning info
!                                 2. 500 plots currents
!..debug
!     write(nterm,8833) ii, ymin,ymax
!     write(nout, 8833) ii, ymin,ymax
 8833 format(" debug ii,ymin,ymax =",i4,1p2e10.2)
!..debug
!
!
!
!.....produce a plot for ii .le. 60
!
!
      go to(10,20,30,40,10,20,10,20,40,10,                               &  
     &      20,10,20,10,20,10,20,10,20,10,                               &  
     &      31,41,20,10,20,51,10,53,20,52,                               &  
     &      10,56,20,10,20,10,20,10,20,10,                               &  
     &      42,20,10,20,10,55,56,53,10,56,                               &  
     &      20,54,10,56,20,54,10,56,20,54),ii
   10 continue
!..rxw:
      nopl1 = nopl1+1
      nopl2 = 1
      if(noplot(50+nopl1).gt.0) goto 400
      nopl2 = 0
!
      if(ii.eq.16 .and. ymin.gt.0) then
      call mapgsl(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1._R8)
      else
      xptmin = 0.7_R8
      if(ii.eq.45 .and. isvd.le.0) xptmin = .290_R8
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,xptmin,1._R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      endif
      call tracec(1h*,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
!
      if(ii.eq.10) then
      call disk(235,1,big1,big2)
      call tracec(1hS,big2(1),big1(1),nrecord,-1,-1,0.,0.)
        if(whalos .gt. 0.) then
          call disk(237,1,big1,big2)
          call tracec(1hH,big2(1),big1(1),nrecord,-1,-1,0.,0.)
        endif
      endif
!
!     if(ii.eq.12 .and. (idata.eq.9 .or. idata.eq.11)) then
!     call disk(217,1,big1,big2)
!     call tracec(1hE,big2,big1,nrecord,-1,-1,0._R8,0._R8)
!     endif
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(ii))
!
      if(ii.eq.14 .and. isurf.ne.0) then
      call disk(172,173,big1,big3)
      call tracec(1h2,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call tracec(1h3,big2,big3,nrecord,-1,-1,0._R8,0._R8)
      call disk(75,1,big1,big2)
      call tracec(1hm,big2,big1,nrecord,-1,-1,0.,0.)
                                    endif
      call setold(tmaxx,ymax,1,0,1,0)
      if(ii.eq.5.and.(acoef(61).ne.0 .or. isym.ne.1 .or. lrswtch.eq.0))  &  
     &                      go to 21
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
      if(ii.ge.49 .and. ii.le.60) then
      write(s100,2024)
      call gtextm(s100,80,0,1,5)
 2024 format(//," *...actual  ",/,                                       &  
     &          " 0...preprogrammed",/,                                  &  
     &         " (kA equivalent) ")
      endif
!
      if(ii.eq.36 .and. irfp.ne.1)go to 731
      if(ii.eq.38 .and. irfp.ne.1 .and. isvd.ne.0)go to 741
      if(ii.eq.18) go to 751
      if(ii.eq.16) go to 761
      if(ii.eq.24) go to 762
!
      if(ii.eq.40) then
      write(s100,2022)
      call gtextm(s100,80,0,1,4)
 2022 format(/," R..resistive",/,                                        &  
     &       " *..res+internal",/,                                       &  
     &       " T..res+int+ext ")
      riendsv = big1(nrecord)
      ymaxrie = ymax
      yminrie = ymin
      endif
!
      if(ii.ne.45) go to 400
      write(s100,2021)
      call gtextm(s100,80,0,1,10)
 2021 format(" *..total heating",/,                                      &  
     &       " A..alpha's",/,                                            &  
     &       " O..ohmic",/,                                              &  
     &       " F..ICRH",/,                                               &  
     &       " E..ECRH",/,                                               &  
     &       " L..PLOSS",/,                                              &  
     &       " B..beams",/,                                              &  
     &       " R..brem-r",/,                                             &  
     &       " C..cycl-r",/,                                             &  
     &       " m..impur-r",/,                                            &  
     &       " H..lower-hybrid")
!
      go to 400
   21 continue
      write(s100,1021)
      call gtextm(s100,80,0,1,3)
 1021 format(/,                                                          &  
     &"  0 initializaton" ,/,                                            &  
     &" -n  nth x-point" ,/,                                             &  
     &" >0 limiter number" )
      go to 400
   20 continue
      if(nopl2.eq.1) goto 400
      if(ii.eq.19) go to 781
      if(ii.eq.23 .and. irfp.ne.1) go to 791
      if(ii.ne.15 .or. irfp.eq.1) then
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      else
      ymax = ymin + 10._R8
      call mapgsl(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,        &
     & .584_R8)
      endif
      call tracec(1h*,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(ii))
      write(nsc1,1088) nlab(ii-1),nlab(ii)
 1088 format(3x,a8," and ",a8," vs time")
 1089 format(3x,a8," vs time (2) for CHI ")
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      if(ii.ge.49.and.ii.le.60) then
      write(s100,2024)
      call gtextm(s100,80,0,1,5)
                                endif
 1023 format(/," min ",1pe12.4,/,                                        &  
     &         " max ",1pe12.4)
      if(ii.eq.17) go to 811
      if(ii.eq.15 .and. irfp.ne.1) go to 701
      if(ii.eq.37 .and. irfp.ne.1) go to 711
      if(ii.eq.39 .and. irfp.ne.1 .and. isvd.ne.0) go to 721
!
      if(ii.eq.42) then
      write(s100,2023) xplas,zplas
      call gtextm(s100,80,0,1,5)
 2023 format(//," *..total from all ",/,                                 &  
     &       "    coils at x,z=  ",/,                                    &  
     &       5x,0p2f6.2 )
      endif
      if(ii.eq.6 .and. isym.eq.1 .and. acoef(61).eq.0) then
      write(s100,3023) xplas,zplas
      call gtextm(s100,80,0,1,4)
 3023 format(//," *...grid boundary",/,                                  &  
     &          " x...",0p2f6.2)
        call disk(151,1,big1,big2)
        call tracec(1hx,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
                                                       endif
!
      if(ii.ne.2 .and. ii.ne.8.and.ii.ne.25.and.ii.ne.29                 &  
     &       .and. ii.ne.51.and.ii.ne.55 .and. ii.ne.59) call frscj(6)
      if(ii.eq.8) then
      write(s100,1024)
      call gtextm(s100,80,0,1,3)
                  endif
 1024 format("  *..line integral",/,                                     &  
     &       "  T..area integral",/,                                     &  
     &       "  0..desired value")
      if(ii.eq.8 .and. lrswtch.ne.0) call frscj(6)
      if(ii.eq.6 .and. isym.eq.0) go to 981
      if(ii.eq.39 .and. irfp.ne.1) go to 921
      if(ii.ne.2) go to 400
      write(s100,1022)
      call gtextm(s100,80,0,1,3)
 1022 format(                                                            &  
     &" p.poloidal energy",/,                                            &  
     &" t.toroidal energy",/,                                            &  
     &" *.total energy")
      go to 400
   30 if(nopl2.eq.1) goto 400
      call tracec(1hp,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      go to 400
   31 if(nopl2.eq.1) goto 400
      call tracec(1hb,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      go to 400
   40 if(nopl2.eq.1) goto 400
      call tracec(1ht,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      if(ii.eq.9) then
      call tracec(1h0,big2(1),apld0sv(1),nrecord,                        &
     &                   -1,-1,0._R8,0._R8)
!
!.....new for CHI plots  5/26/05
      call disk(236,1,big1,big2)
      call tracec(1hS,big2(1),big1(1),nrecord,-1,-1,0.,0.)
      endif
      call frscj(6)
!
!.....special coding to make plot for CHI simulations
      if(ii.eq.9 .and. globmax(220) .gt. globmin(220)) then
!
      ymin = globmin(220)
      ymax = globmax(220)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.00_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(221,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(220))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
      ymin = globmin(221)
      ymax = globmax(221)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(222,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(221))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
      write(nsc1,1089) nlab(220)
      call frscj(6)
      endif
!.....end of special coding for CHI
!
!
!......special plots for porcelli sawtooth model
      if(ii.eq.9 .and. isaw.eq.3) then
!
      ymin = min(globmin(222),globmin(223))
      ymax = max(globmax(222),globmax(223))
      if(ymax .gt. .05_R8) ymax = .05_R8
      if(ymin .lt. -.10_R8) ymin = -.10_R8
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.00_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(223,1,big1,big2)
      call tracec(1hW,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(224,1,big1,big2)
      call tracec(1hc,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(222))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
      ymin = min(globmin(224),globmin(226))
      ymax = max(globmax(224),globmax(226))
      if(ymax .gt. .05_R8) ymax = .05_R8
      if(ymin .lt. -.10_R8) ymin = -.10_R8
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
!     call setold(tmaxx,ymin,1,0,1,0)
!     write(s100,6666)
!     call gtext(s100,80,0)
      call sawdraw(s100,tmin,tmaxx,ymin,ymax)
      call disk(225,1,big1,big2)
      call tracec(1hW,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(227,1,big1,big2)
      call tracec(1hC,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(224))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
      write(nsc1,1090)
 1090 format(3x,"delta-W for Porcelli model")
      call frscj(6)
      endif
!......end of special plots for porcelli sawtooth model
!
!
      go to 400
   41 if(nopl2.eq.1) goto 400
      call tracec(1hl,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      go to 400
   51 if(nopl2.eq.1) goto 400
      call tracec(1h2,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      write(s100,1151)
      call setold(tmaxx,ymin+0.4_R8*(ymax-ymin),1,0,1,0)
      call gtextm(s100,80,0,1,4)
 1151 format(//," *..ener/power",/,                                      &  
     &          " 2..ener/hflux")
      call frscj(6)
      go to 400
   52 if(nopl2.eq.1) goto 400
      call tracec(1hi,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      call frscj(6)
      go to 400
   53 if(nopl2.eq.1) goto 400
      call tracec(1hi,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
!     if(idata.eq.9 .or. idata.eq.11) then
!     call disk(218,1,big1,big2)
!     call tracec(1hE,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
!     endif
      if(ii.eq.48) go to 771
      go to 400
   42 if(nopl2.eq.1) goto 400
      call tracec(1hr,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      rendsv = big1(nrecord)
      call disk(43,45,big1,big3)
      do 2242 ic=1,nrecord
 2242 big4(ic) = big1(ic) - big1(1) + big3(ic)
      call tracec(1ht,big2(1),big4(1),nrecord,-1,-1,0._R8,0._R8)
      rieendsv = big4(nrecord)
      write(s100,2243) rendsv,riendsv,rieendsv
      call setold(tmaxx,yminrie + 0.4_R8*(ymaxrie-yminrie),1,0,1,0)
      call gtextm(s100,80,0,1,4)
 2243 format(/," R    ",1pe11.3,/,                                       &  
     &         " R+I  ",1pe11.3,/,                                       &  
     &         " R+I+E",1pe11.3   )
      go to 400
   54 if(nopl2.eq.1) goto 400
      call tracec(1h0,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      call frscj(6)
      go to 400
!
!.....special for lings f-theta plots
  810 call disk(ii+1,1,big1,big2)
      ymin = min(globmin(24),globmin(25))
      ymax = max(globmax(24),globmax(25))
      call maps(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1._R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call tracec(1hf,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd('f and th')
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) globmin(24),globmax(24)
      call gtextm(s100,80,0,1,3)
      write(s100,1023) globmin(25),globmax(25)
      call gtextm(s100,80,0,1,3)
      go to 400
  820 call disk(ii+1,1,big1,big2)
      call tracec(1ht,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
!
      call maps(globmin(25),globmax(25),globmin(24),globmax(24),         &  
     &          .142_R8,.800_R8,.290_R8,.584_R8)
      call disk(26,25,big1,big2)
      call trace(big1,big2,nrecord,-1,-1,0._R8,0._R8)
      n10 = nrecord/10 + 1
      icc = 0
      tplot = big3(1)
      dtplot = (big3(nrecord)-big3(1))/10._R8
      do 2810 ic=1,nrecord
      if(big3(ic).lt.tplot) go to 2810
      tplot = tplot + dtplot
      icc = icc+1
      call setold(big1(ic),big2(ic),1,0,1,0)
      call gtext(lab2(icc),1,0)
 2810 continue
      call tracec(1h*,big1(1),big3(1),nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(' f vs th ')
      write(nsc1,1820)
 1820 format(3x,45h f and th vs time and f vs theta              )
      call frscj(6)
      go to 400
!
55    if(nopl2.eq.1) goto 400
      call tracec(1ha,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      call disk(242,1,big1,big2)
      call tracec(1he,big2,big1,nrecord,-1,-1,0._R8,0._R8)      
      call disk(120,1,big1,big2)
      call tracec(1hr,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(121,1,big1,big2)
      call tracec(1hc,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(122,1,big1,big2)
      call tracec(1hm,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(123,1,big1,big2)
      call tracec(1hf,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(124,1,big1,big2)
      call tracec(1hb,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(125,1,big1,big2)
      call tracec(1hH,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      goto 400
56    if(nopl2.eq.1) goto 400
      call tracec(1ho,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      goto 400
  701 continue
      write(s100,1701)
      call gtextm(s100,80,0,1,4)
      write(98,1998) (big1(i),i=1,nrecord)    ! debug
      call disk(62,1,big1,big2)
      write(98,1998) (big1(i),i=1,nrecord)    ! debug
      call tracec(1h5,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      if(isurf.ne.0) then
      call disk(63,1,big1,big2)
      write(98,1998) (big1(i),i=1,nrecord)    ! debug
      call tracec(1hr,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(64,1,big1,big2)
      write(98,1998) (big1(i),i=1,nrecord)    ! debug
 1998 format(1p10e12.4)     ! debug
      call tracec(1hc,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      endif
      call frscj(6)
      if(isurf.eq.0) go to 400
      ymin = globmin(151)
      ymax = globmax(151)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.00_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(152,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(151))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
      ymin = globmin(152)
      ymax = globmax(152)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(153,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(152))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(nsc1,1088) nlab(151),nlab(152)
      call frscj(6)
      go to 400
 1701 format(" *..edge q",/,                                             &  
     &       " 5..95%  q",/,                                             &  
     &       " R..q-star",/,                                             &  
     &       " C..q-cylin")
  711 continue
      write(s100,1711)
      call gtextm(s100,80,0,1,3)
      call disk(66,1,big1,big2)
      call tracec(1h5,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(82,1,big1,big2)
      call tracec(1h9,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call frscj(6)
      go to 400
 1711 format(" *..edge",/,                                               &  
     &       " 5..95% ",/," 9..90%  ")
  721 continue
      write(s100,1721)
      call gtextm(s100,80,0,1,2)
      call disk(68,1,big1,big2)
      call tracec(1hC,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call frscj(6)
      go to 400
 1721 format(" *..true value",/,                                         &  
     &       " C..calculated")
  921 continue
      if(isym.eq.1) go to 1601
      if(globmax(147) .le. globmin(147)) go to 1601
!.....ii=147   xsep(2)
!
      ymin = globmin(147)
      ymax = globmax(147)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.00_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(148,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(147))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
!.....ii=148  zsep(2)
      ymin = globmin(148)
      ymax = globmax(148)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(149,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(148))
!
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(nsc1,1088) nlab(147),nlab(148)
      call frscj(6)
 1601 continue
      nmax = min(10,ncnt)
      if(ncnt.le.0) go to 400
      do 1602 i=1,nmax
      indx = 178 + 2*(i-1)
      if(globmax(indx) .le. globmin(indx)) go to 1601
!.....ii=indx   delta-psi
!
      ymin = globmin(indx)
      ymax = globmax(indx)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.00_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(indx+1,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
!     call crtbcd(nlab(indx))
      if(acoef(296).eq.6 .or. acoef(296).eq.7) then
      write(s100,1613) i
 1613 format("gap dist",i2)
      else
      write(s100,1603) i
 1603 format("flux diff:",i2)
      endif
      call gtextm(s100,80,0,1,1)
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
!.....ii=indx+1  delta-X
      ymin = max(globmin(indx+1),-0.2_R8)
      ymax = min(globmax(indx+1), 0.2_R8)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(indx+2,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
!     call crtbcd(nlab(indx+1))
      if(acoef(296).eq.6 .or. acoef(296).eq.7) then
      write(s100,1614) i
 1614 format(" gap derivative ",i3)
      else
      write(s100,1604) i
 1604 format(" distance(m)",i2)
      endif
      call gtextm(s100,80,0,1,1)
!
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!     write(nsc1,1088) nlab(indx),nlab(indx+1)
      write(nsc1,1605) i
 1605 format(" flux difference and distance",i2)
      call frscj(6)
 1602 continue
      go to 400
!
  731 continue
      write(s100,1731)
      call gtextm(s100,80,0,1,3)
      call disk(65,1,big1,big2)
      call tracec(1h5,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(83,1,big1,big2)
      call tracec(1h9,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      go to 400
 1731 format(" *..edge",/,                                               &  
     &       " 5..95% ",/," 9..90% ")
  741 continue
      write(s100,1741)
      call gtextm(s100,80,0,1,2)
      call disk(67,1,big1,big2)
      call tracec(1hc,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      go to 400
 1741 format(" *..true value",/,                                         &  
     &       " C..calculated")
  751 continue
      write(s100,1751)
      call gtextm(s100,80,0,1,3)
      call disk(71,1,big1,big2)
      call tracec(1hJ,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(150,1,big1,big2)
      call tracec(1hG,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      go to 400
 1751 format(" *..Murakami ",/,                                          &  
     &       " J..Jet-crit ",/,                                          &  
     &       " G..Greenwald"  )
  761 continue
      write(s100,1761)
      call gtextm(s100,80,0,1,4)
      call disk(72,1,big1,big2)
      call tracec(1h1,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(73,1,big1,big2)
      call tracec(1h2,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(74,1,big1,big2)
      call tracec(1h3,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      go to 400
 1761 format(" *..actual ",/,                                            &  
     &       " 1..resistive",/,                                          &  
     &       " 2..shear wave",/,                                         &  
     &       " 3..fast wave")
!
  762 continue
      write(s100,1762)
      call gtextm(s100,80,0,1,4)
      call disk(132,1,big1,big2)
      call tracec(1hg,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(133,1,big1,big2)
      call tracec(1ha,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(134,1,big1,big2)
      call tracec(1hr,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      go to 400
 1762 format(" *..Kaye-Gold",/,                                          &  
     &       " G..Goldston ",/,                                          &  
     &       " A.ITER98(y,2",/,                                          &  
     &       " R..Rebut-Lal" )
!
  771 continue
      if(isvd.gt.0) go to 772
      write(nsc1,1772)
 1772 format(" input & radiated power vs time")
      call frscj(6)
      go to 897
  772 continue
      xmin = globmin(68)
      xmax = globmax(68)
      ymin = globmin(69)
      ymax = globmax(69)
      diffx = xmax-xmin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffx.le.0) xmax = xmin+max(1._R8,abs(xmin))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin+max(1._R8,abs(ymin))
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call disk(69,70,big1,big2)
      call trace(big1,big2,nrecord,-1,-1,0._R8,0._R8)
      n10 = nrecord/10 + 1
      icc = 0
      tplot = big3(1)
      dtplot = (big3(nrecord)-big3(1))/10._R8
      do 2771 ic=1,nrecord
      if(big3(ic).lt.tplot) go to 2771
      tplot = tplot + dtplot
      icc = icc+1
      call setold(big1(ic),big2(ic),1,0,1,0)
      call gtext(lab2(icc),1,0)
 2771 continue
      call setld(2._R8,18._R8,0,0,2,1)
      write(s100,1771)
      call gtext(s100,80,0)
 1771 format(" psepcal (y) vs psisep (x)")
      write(nsc1,1088) nlab(45),nlab(69)
      call frscj(6)
!
      if(isvd.le.0) go to 897
      ymin = min(globmin(85),globmin(86),globmin(87),-.10_R8)
      ymax = max(globmax(85),globmax(86),globmax(87),0._R8)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin+max(1._R8,abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.1_R8,.47_R8,.72_R8,1._R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(86,1,big1,big2)
      call tracec(1ha,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(87,1,big1,big2)
      call tracec(1hc,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(88,1,big1,big2)
      call tracec(1h0,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setold(tmin-(tmaxx-tmin)*.20_R8,ymin,1,0,1,1)
      write(s100,3301)
      call gtext(s100,80,0)
 3301 format(" outer midplane ")
      ymin = min(globmin(88),globmin(89),globmin(90),0._R8)
      ymax = max(globmax(88),globmax(89),globmax(90),.10_R8)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin+max(1._R8,abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.59_R8,.96_R8,.72_R8,1._R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(89,1,big1,big2)
      call tracec(1ha,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(90,1,big1,big2)
      call tracec(1hc,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(91,1,big1,big2)
      call tracec(1h0,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setold(tmaxx+(tmaxx-tmin)*.05_R8,ymin,1,0,1,1)
      write(s100,3302)
      call gtext(s100,80,0)
 3302 format(" inner midplane ")
      ymin = min(globmin(91),globmin(92),globmin(93))
      ymax = max(globmax(91),globmax(92),globmax(93))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin+max(1._R8,abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.1_R8,.47_R8,.35_R8,.62_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(92,1,big1,big2)
      call tracec(1ha,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(93,1,big1,big2)
      call tracec(1hc,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(94,1,big1,big2)
      call tracec(1h0,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setold(tmin-(tmaxx-tmin)*.20_R8,ymin,1,0,1,1)
      write(s100,3303)
      call gtext(s100,80,0)
 3303 format(" inner strikepoint ")
      call setold(tmin,ymin-.3_R8*(ymax-ymin),1,0,1,0)
      write(s100,3403)
      call gtext(s100,80,0)
 3403 format(" critical surface points:  A-actual ,  C-calculated ",     &  
     &       ",  0-desired " )
      ymin = min(globmin(94),globmin(95),globmin(96))
      ymax = max(globmax(94),globmax(95),globmax(96))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin+max(1._R8,abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.59_R8,.96_R8,.35_R8,.62_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(95,1,big1,big2)
      call tracec(1ha,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(96,1,big1,big2)
      call tracec(1hc,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(97,1,big1,big2)
      call tracec(1h0,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setold(tmaxx+(tmaxx-tmin)*.05_R8,ymin,1,0,1,1)
      write(s100,3304)
      call gtext(s100,80,0)
 3304 format(" outer strikepoint ")
      write(nsc1,3305)
 3305 format(" critical errors vs time")
      call frscj(6)
!
  897 continue
      xmin = 0._R8
      xmax = globmax(80)
      ymin = 0._R8
      ymax = globmax(138)
      diffx = xmax-xmin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffx.le.0) xmax = xmin + max(1._R8,abs(xmin))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymin))
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call disk(81,139,big1,big2)
      call trace(big1,big2,nrecord,-1,-1,0._R8,0._R8)
      icc = 0
      tplot = big3(1)
      dtplot = (big3(nrecord)-big3(1))/10._R8
      do 3791 ic=1,nrecord
      if(big3(ic).lt.tplot) go to 3791
      tplot = tplot + dtplot
      icc = icc+1
      call setold(big1(ic),big2(ic),1,0,1,0)
      call gtext(lab2(icc),1,0)
 3791 continue
      call setld(2._R8,18._R8,0,0,2,1)
      write(s100,3782)
      call gtext(s100,80,0)
 3782 format(" W-dot (MWATT) vs W (MJOULE) ")
      call setold(xmax,ymax,1,0,1,0)
      write(s100,3792) xmin,xmax,ymin,ymax
      call gtextm(s100,80,0,1,5)
 3792 format(/," xmin",1pe12.4,/,                                        &  
     &         " xmax",1pe12.4,/,                                        &  
     &         " ymin",1pe12.4,/,                                        &  
     &         " ymax",1pe12.4)
!.....ii=139  taue/H_98
      ymin = min(globmin(139),globmin(140))
      ymax = max(globmax(139),globmax(140))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1151)
      call gtextm(s100,80,0,1,4)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(140,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(141,1,big1,big2)
      call tracec(1h2,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(139))
!
      write(nsc1,4028)
 4028 format(" W-dot vs W and   H_98  vs time")
      call frscj(6)
!
!......power flow plot
!
!
      ymin = min(globmin(97),globmin(98),globmin(99),globmin(100),       &  
     &             globmin(101),globmin(102),globmin(103),globmin(104),  &  
     &             globmin(105))
      ymax = max(globmax(97),globmax(98),globmax(99),globmax(100),       &  
     &             globmax(101),globmax(102),globmax(103),globmax(104),  &  
     &             globmax(105))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) go to 903
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
 6666 format(" time(sec)")
      do 898 ip=97,105
      call disk(ip+1,1,big1,big2)
      nn = ip-96
  898 call tracec(lab2(nn),big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,40._R8,0,0,2,1)
      write(s100,1898)
      call gtext(s100,80,0)
 1898 format("   power flow .. MWATT")
      write(nsc1,1897)
 1897 format(" power and energy flow-Poynting")
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1899)
      call gtextm(s100,80,0,1,10)
      ymin = min(globmin(110),globmin(111),globmin(112),globmin(113),    &  
     &             globmin(114),globmin(115),globmin(116),globmin(117),  &  
     &             globmin(118))
      ymax = max(globmax(110),globmax(111),globmax(112),globmax(113),    &  
     &             globmax(114),globmax(115),globmax(116),globmax(117),  &  
     &             globmax(118))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      do 899 ip=110,118
      call disk(ip+1,1,big1,big2)
      nn = ip-109
  899 call tracec(lab2(nn),big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,15._R8,0,0,2,1)
      write(s100,1901)
      call gtext(s100,80,0)
 1901 format("energy flow..Mjoules")
 1899 format(/,                                                          &  
     &" 1 ... thermal  ",/,                                              &  
     &" 2 ... pol field",/,                                              &  
     &" 3 ... tor field",/,                                              &  
     &" 4 ... heat flux",/,                                              &  
     &" 5 ... heat source",/,                                            &  
     &" 6 ... int coil   ",/,                                            &  
     &" 7 ... pf boundary",/,                                            &  
     &" 8 ... tf boundary",/,                                            &  
     &" 9 ... total      " )
      call frscj(6)
!
!.....ii=106   tevv
!
  903 continue
      iip = 106
      ymin = globmin(106)
      ymax = globmax(106)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
!..debug
!     write(nterm,8833) iip, ymin,ymax
!     write(nout, 8833) iip, ymin,ymax
      call mapgsl(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call disk(107,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(106))
!
!.....ii=107  ffac
      iip = 107
      ymin = globmin(107)
      ymax = globmax(107)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
!..debug
!     write(nterm,8833) iip, ymin,ymax
!     write(nout, 8833) iip, ymin,ymax
      call mapgsl(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call disk(108,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(107))
!
      write(nsc1,3017)
 3017 format(" ffac and tevv vs time ")
      call frscj(6)
!
!.....ii=165  whalo
      iip = 165
      ymin = globmin(165)
      ymax = globmax(165)
!     if(ymax.le.ymin) go to 3218
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
!..debug
!     write(nterm,8833) iip, ymin,ymax
!     write(nout, 8833) iip, ymin,ymax
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call disk(166,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(165))
!
!.....ii=164   thalo
!
      iip = 164
      ymin = globmin(164)
      ymax = globmax(164)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
!..debug
!     write(nterm,8833) iip, ymin,ymax
!     write(nout, 8833) iip, ymin,ymax
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call disk(165,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(164))
!
      write(nsc1,3117)
 3117 format(" whalo and thalo vs time ")
      call frscj(6)
!
!.....ii=167  pohm-H
      iip = 167
      ymin = globmin(167)
      ymax = globmax(167)
      if(ymax.le.ymin) go to 3218
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
!..debug
!     write(nterm,8833) iip, ymin,ymax
!     write(nout, 8833) iip, ymin,ymax
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call disk(168,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(167))
!
!.....ii=166   Pohm-P
!
      iip = 166
      ymin = globmin(166)
      ymax = globmax(166)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
!..debug
!     write(nterm,8833) iip, ymin,ymax
!     write(nout, 8833) iip, ymin,ymax
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call disk(167,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(166))
!
      write(nsc1,3118)
 3118 format(" Pohm-H and Pohm-P vs time ")
      call frscj(6)
 3218 continue
!
!.....ii=237   Plasma and halo
!
      ymin = min(globmin(237),globmin(238))
      ymax = max(globmax(237),globmax(238))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(238,1,big1,big2)
      call tracec(1hP,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(239,1,big1,big2)
      call tracec(1hH,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(237))
!
!.....ii=238  Plasma - halo
      ymin = globmin(239)
      ymax = globmax(239)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(240,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(239))
!
      write(nsc1,3020)
 3020 format(" Plasma and Halo Area vs time ")
      call frscj(6)
      if(isurf.eq.0 .or. lrswtch.ne.0) go to 979
!
!.....ii=108   npsit
!
      iip = 108
      ymin = globmin(108)
      ymax = globmax(108)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
!     write(nterm,8833) iip, ymin,ymax
!     write(nout, 8833) iip, ymin,ymax
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(109,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(108))
!
!.....ii=109  resid
      ymin = 0._R8
      ymax = globmax(109)
      if(ymax.gt.2) ymax = 2._R8
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(110,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(109))
!
      write(nsc1,3018)
 3018 format(" npsit and resid vs time ")
      call frscj(6)
!
!.....ii=78   volume averaged density
!
      ymin = min(globmin(78),globmin(136),globmin(137))
      ymax = max(globmax(78),globmax(136),globmax(137))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(79,1,big1,big2)
      call tracec(1hE,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(137,1,big1,big2)
      call tracec(1hI,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(138,1,big1,big2)
      call tracec(1hH,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(78))
!
!.....ii=135  zeff
      ymin = globmin(135)
      ymax = globmax(135)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(136,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(135))
!
      write(nsc1,3028)
 3028 format(" vol-dens & zeff vs time ")
      call frscj(6)
      if(ialpha.eq.0) go to 979
!
!.....ii=129   n-alpha
!
      ymin = globmin(129)
      ymax = globmax(129)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(130,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(129))
!
!.....ii=130  beta-alpha
      ymin = globmin(130)
      ymax = globmax(130)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(131,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(130))
!
      write(nsc1,3029)
 3029 format(" n-alpha,b-alpha vs time ")
      call frscj(6)
      if(acoef(114).eq.0) go to 979
!
!.....ii=175   Total and alpha power
!
      ymin = globmin(175)
      ymax = globmax(175)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(176,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(177,1,big1,big2)
      call tracec(1hA,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(176))
!
!.....ii=177  acoef(113)...fraction of tritium
      ymin = globmin(177) - 0.01_R8
      ymax = globmax(177) + 0.01_R8
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(178,1,big1,big2)
      call tracec(1hT,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(177))
!
      write(nsc1,3030)
 3030 format(" Palpha,T-fract  vs time ")
      call frscj(6)
  979 continue
      if(isurf.eq.0 .or. lrswtch.ne.0) go to 989
      if(ibootst.eq.0._R8.and. fracpar.eq.0._R8.and. acoef(1).ne.1._R8   &  
     &  .and. ilhcd .eq. 0 .and. acoef(296).ne.5._R8)  go to 989
!
!.....current drive and bootstrap current
      ymin = min(globmin(141),globmin(142),globmin(146),globmin(198),    &  
     &  globmin(213),globmin(9), globmin(242))
      ymax = max(globmax(141),globmax(142),globmax(146),globmax(198),    &  
     &  globmax(213),globmax(9), globmin(242))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy .le. 0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(142,1,big1,big2)
      call tracec(1hB,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(143,1,big1,big2)
      call tracec(1hN,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(147,1,big1,big2)
      call tracec(1hL,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(199,1,big1,big2)
      call tracec(1hF,big2,big1,nrecord,-1,-1,0._R8,0._R8)
!...added 01/28/11..cj
      call disk(243,1,big1,big2)
      call tracec(1hE,big2,big1,nrecord,-1,-1,0._R8,0._R8)
!...added 10/17/98..scj
      call disk(9,1,big1,big2)
      call tracec(1hP,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(214,1,big1,big2)
      call tracec(1hS,big2,big1,nrecord,-1,-1,0._R8,0._R8)
!...added 02/14/03
      call disk(234,1,big1,big2)
      call tracec(1hG,big2,big1,nrecord,-1,-1,0._R8,0._R8)
!
      call setld(2._R8,22._R8,0,0,2,1)
      write(s100,1985)
      call gtext(s100,80,0)
 1985 format(" current drive ")
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(s100,1986)
      call gtextm(s100,80,0,1,8)
 1986 format(/," B...BOOTSTRAP",/,                                       &  
     &         " N...N. BEAM CUR",/,                                     &  
     &         " F..... FWCD CUR",/,                                     &  
     &         " L...LOWER HYBRID",/,                                    &  
     &         " E..... ECCD CUR",/,                                     &  
     &         " S...SUM of CD   ",/,                                    &  
     &         " P...Plasma Current",/,                                  &  
     &         " G..pressure driven")
!
!.....pedistal T
      ymin = 0._R8
      ymax = max(globmax(218),globmax(219))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(219,1,big1,big2)
      call tracec(1hE,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(220,1,big1,big2)
      call tracec(1hI,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      write(s100,1987)
      call gtext(s100,80,0)
 1987 format(" pedistal T ")
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(s100,1731)
      call gtextm(s100,80,0,1,3)
      write(nsc1,1989)
 1989 format(" current drive and pedistal T")
      call frscj(6)
  989 continue
!
!
!...II=125  CPU time(min)
      ymin = globmin(125)
      ymax = globmax(125)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + 1
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(126,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(125))
!
!...II=126  CPU rate(min/sec)
!
      ymin = globmin(126)
      ymax = globmax(126)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      if(ymin.le.0) ymin = 1.E-5_R8*ymax
      iip = 126
!..debug
!     write(nterm,8833) iip, ymin,ymax
!     write(nout, 8833) iip, ymin,ymax
      call mapgsl(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call disk(127,1,big1,big2)
      do 8981 i=1,nrecord
      if(big1(i).lt.ymin) big1(i) = ymin
      if(big1(i).gt.ymax) big1(i) = ymax
 8981 continue
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(126))
      write(nsc1,3019)
 3019 format(" cpu time and cpu rate vs time")
      call frscj(6)
!
      if(isurf.eq.0.or.lrswtch.ne.0) go to 400
!
!.....ii=207 iqtrubmax....indicator of trouble in flxvol
      ymin = globmin(207)
      ymax = globmax(207)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(208,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(207))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
!.....ii=208 nskipsf .... cycles skipped between calls to flxvol
      ymin = globmin(208)
      ymax = globmax(208)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(209,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,18._R8,0,0,2,1)
      call crtbcd(nlab(208))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(nsc1,7212)
 7212 format("iqtrubmax and nskipsf vs time")
      call frscj(6)
!
      if(irunaway.gt.0) then
!.....ii=209  recur and apl runaway current
      ymin = min(globmin(209),globmin(210))
      ymax = max(globmax(209),globmax(210))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(210,1,big1,big2)
      call tracec(1hR,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(211,1,big1,big2)
      call tracec(1hP,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(209))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
!.....ii=211 sreav and sumre...runaway source terms
      ymin = min(globmin(211),globmin(212))
      ymax = max(globmax(211),globmax(212))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(212,1,big1,big2)
      call tracec(1hA,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(213,1,big1,big2)
      call tracec(1hT,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,18._R8,0,0,2,1)
      call crtbcd(nlab(211))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(nsc1,7211)
 7211 format("runaway current and source terms")
      call frscj(6)
      endif
!
!
!
!.....ii=129   volt-sec evaluation using Poynting Method
!
      ymin = min(globmin(168),globmin(169))
      ymax = max(globmax(168),globmax(169))
      if(ymax.le.ymin) go to 400
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(169,1,big1,big2)
      call tracec(1hI,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      riendsv = big1(nrecord)
      ristrsv = big1(1)
      call disk(170,1,big1,big2)
      call tracec(1hR,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      rendsv = big1(nrecord)
      rstrsv = big1(1)
      rieendsv = riendsv + rendsv - ristrsv - rstrsv
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(nlab(168))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(s100,2244) rendsv,riendsv,rieendsv
      call setold(tmaxx,ymin+ 0.4_R8*(ymax -ymin ),1,0,1,0)
      call gtextm(s100,80,0,1,4)
 2244 format(/," R    ",1pe11.3,/,                                       &  
     &         " I    ",1pe11.3,/,                                       &  
     &         " D-R+I",1pe11.3   )
!
!.....ii=170  Ejima Coeff
      ymin = globmin(170)
      ymax = globmax(170)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(171,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(170))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
      write(nsc1,4029)
 4029 format(" Poynting Volt-Sec(R&I)  ")
      call frscj(6)
!
      go to 400
!
!.....zmag vs xmag
  981 continue
      xmin = globmin(7)
      xmax = globmax(7)
      ymin = globmin(6)
      ymax = globmax(6)
      diffx = xmax-xmin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffx.le.0) xmax = xmin + max(1._R8,abs(xmin))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymin))
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.800_R8,.290_R8,1.00_R8)
      call disk(8,7,big1,big2)
      call trace(big1,big2,nrecord,-1,-1,0._R8,0._R8)
      icc = 0
      tplot = big3(1)
      dtplot = (big3(nrecord)-big3(1))/10._R8
      do 2991 ic=1,nrecord
      if(big3(ic).lt.tplot) go to 2991
      tplot = tplot + dtplot
      icc = icc+1
      call setold(big1(ic),big2(ic),1,0,1,0)
      call gtext(lab2(icc),1,0)
 2991 continue
      call setld(2._R8,18._R8,0,0,2,1)
      write(s100,1982)
      call gtext(s100,80,0)
 1982 format(" zmag (y) vs xmag (x)")
      write(nsc1,1981) nlab(ii+1),nlab(ii)
 1981 format(3x,a8," vs time and ",a8)
      call setold(xmax,ymax,1,0,1,0)
      write(s100,1991) xmin,xmax,ymin,ymax
      call gtextm(s100,80,0,1,5)
 1991 format(/," xmin",1pe12.4,/,                                        &  
     &         " xmax",1pe12.4,/,                                        &  
     &         " ymin",1pe12.4,/,                                        &  
     &         " ymax",1pe12.4)
      call frscj(6)
      go to 400
!
!.....ii=19    1/q vs nr/b
  781 continue
      xmin = 0._R8
      xmax = max(globmax(19),1.E20_R8)
      ymin = 0._R8
      ymax = max(globmax(74),0.5E-20_R8*xmax)
      diffx = xmax-xmin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffx.le.0) xmax = xmin + max(1._R8,abs(xmin))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymin))
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call disk(20,75,big1,big2)
      call trace(big1,big2,nrecord,-1,-1,0._R8,0._R8)
      icc = 0
      tplot = big3(1)
      dtplot = (big3(nrecord)-big3(1))/10._R8
      do 2791 ic=1,nrecord
      if(big3(ic).lt.tplot) go to 2791
      tplot = tplot + dtplot
      icc = icc+1
      call setold(big1(ic),big2(ic),1,0,1,0)
      call gtext(lab2(icc),1,0)
 2791 continue
      call setld(2._R8,18._R8,0,0,2,1)
      write(s100,1782)
      call gtext(s100,80,0)
 1782 format(" 1/q (y) vs nR/B (x)")
      write(nsc1,1781) nlab(ii-1),nlab(ii)
 1781 format(3x,a8," vs time and ",a8)
      call setold(xmax,ymax,1,0,1,0)
      write(s100,1791) xmin,xmax,ymin,ymax
      call gtextm(s100,80,0,1,5)
 1791 format(/," xmin",1pe12.4,/,                                        &  
     &         " xmax",1pe12.4,/,                                        &  
     &         " ymin",1pe12.4,/,                                        &  
     &         " ymax",1pe12.4)
!
!.....shade unstable part
      xlist(1) = 0.0_R8
      ylist(1) = 0.0_R8
      xlist(2) = xmax
      ylist(2) = (0.5E-20_R8)*xmax
      xlist(3) = xmax
      ylist(3) = 0.0_R8
!     if(xmax.gt.0)
!    1call hatch(xlist,ylist,3,pi/2.,3,0)
!
      call frscj(6)
!
  782 ymin = min(globmin(77),globmin(78),globmin(79))
      if(globmin(215).gt.0) ymin = min(ymin,globmin(215))
      ymax = max(globmax(77),globmax(78),globmax(79),globmax(215))
      if(ymin == ymin .and. ymax == ymax) then
        diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
        if(diffy.le.0) ymax = ymin+max(1._R8,abs(ymin))
        call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
        call setold(tmaxx,ymin,1,0,1,0)
        write(s100,6666)
        call gtext(s100,80,0)
        call disk(78,1,big1,big2)
        call tracec(1hl,big2,big1,nrecord,-1,-1,0._R8,0._R8)
        call disk(79,1,big1,big2)
        call tracec(1hv,big2,big1,nrecord,-1,-1,0._R8,0._R8)
        call disk(80,1,big1,big2)
        call tracec(1h0,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      endif
!
!     if(idata.eq.9 .or. idata.eq.10 .or. idata.eq.11) then
!     call disk(216,1,big1,big2)
!     if(big1(1) .le. 1.E40_R8) then
!     call tracec(1hE,big2,big1,nrecord,-1,-1,0._R8,0._R8)
!     endif
!     endif
!
      call setld(2._R8,18._R8,0,0,2,1)
      write(s100,1882)
      call gtext(s100,80,0)
 1882 format("density(10**20m-3)")
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1891) ymin,ymax
      call gtextm(s100,80,0,1,8)
 1891 format(/," ymin",1pe12.4,/,                                        &  
     &         " ymax",1pe12.4,/,                                        &  
     &      //," L..line av",                                            &  
     &       /," V..vol ave",                                            &  
     &       /," 0..central" )
      ymin = globmin(80)
      ymax = globmax(80)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin+max(1._R8,abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1._R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(81,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      if(irfp.ne.1) then
      if (acoef(504).eq.0) then
                           write(s100,1883)
      else
                           write(s100,1884)
      endif
      call gtext(s100,80,0)
 1883 format("int energy(MJ)")
 1884 format("mag energy(MJ)")
      else
      write(s100,2883)
      call gtext(s100,80,0)
 2883 format("nloop-in hyper")
      endif
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(nsc1,1088) nlab(79),nlab(80)
      call frscj(6)
      go to 400
!
!....ii=23    li vs q
  791 continue
      xmin = 1.0_R8
      xmax = 8.0_R8
      ymin = 0.1_R8
      ymax = 1.1_R8
      diffx = xmax-xmin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffx.le.0) xmax = xmin + max(1._R8,abs(xmin))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymin))
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call disk(62,24,big1,big2)
      call trace(big1,big2,nrecord,-1,-1,0._R8,0._R8)
!
      icc = 0
      tplot = big3(1)
      dtplot = (big3(nrecord)-big3(1))/10._R8
      do 2792 ic=1,nrecord
      if(big3(ic) .lt. tplot) go to 2792
      tplot = tplot + dtplot
      icc = icc + 1
      call setold(big1(ic),big2(ic),1,0,1,0)
      call gtext(lab2(icc),1,0)
 2792 continue
!
      call setld(2._R8,18._R8,0,0,2,1)
      write(s100,1792)
      call gtext(s100,80,0)
 1792 format("li/2(y) vs q (x)")
      write(nsc1,1781) nlab(ii-1),nlab(ii)
      call setold(xmax,ymax,1,0,1,0)
      write(s100,1791) xmin,xmax,ymin,ymax
      call gtextm(s100,80,0,1,5)
!
!.....shade unstable part
      call hatch(xlist2,ylist2,20,3*pi/4._R8,3,0)
      call frscj(6)
!
      xmin = 0.4_R8
      xmax = 1.5_R8
      ymin = 0.0_R8
      ymax = 5.0_R8
      diffx = xmax-xmin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffx.le.0) xmax = xmin + max(1._R8,abs(xmin))
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymin))
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.642_R8,.084_R8,.584_R8)
      call disk(174,175,big1,big2)
      call trace(big1,big2,nrecord,-1,-1,0._R8,0._R8)
!
      icc = 0
      tplot = big3(1)
      dtplot = (big3(nrecord)-big3(1))/10._R8
      do 2711 ic=1,nrecord
      if(big3(ic) .lt. tplot) go to 2711
      tplot = tplot + dtplot
      icc = icc + 1
      call setold(big1(ic),big2(ic),1,0,1,0)
      call gtext(lab2(icc),1,0)
 2711 continue
!
      call setld(2._R8,18._R8,0,0,2,1)
      write(s100,1712)
      call gtext(s100,80,0)
 1712 format("Ct(y) vs li(x)")
      write(nsc1,1781) nlab(174),nlab(173)
!
!.....ii=174  Troyon Coeff
      ymin = globmin(174)
      ymax = globmax(174)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(175,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(174))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      call frscj(6)
      if(acoef(760) .lt. 1._R8) go to 400
!
!.....ii=199  Pellet Radius
      ymin = globmin(204)
      ymax = globmax(204)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(205,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(204))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
!.....ii=201  Major Radius vs time
      ymin = globmin(201)
      ymax = globmax(201)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      if(diffx.le.0) xmax = xmin + max(1._R8,abs(xmax),abs(xmin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call disk(202,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,18._R8,0,0,2,1)
      write(s100,7203)
      call gtext(s100,80,0)
 7203 format(" pellet penetration vs time")
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(nsc1,7199)
 7199 format(" pellet radius & penetration")
      call frscj(6)
!
!.....ii=200  dndt..pellet ablation rate
      ymin = globmin(200)
      ymax = globmax(200)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(201,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(200))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
!.....ii=202  Te-Pellet
      ymin = globmin(202)
      ymax = globmax(202)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(203,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,18._R8,0,0,2,1)
      call crtbcd(nlab(202))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(nsc1,7200)
 7200 format("pellet dndt and Te vs time")
      call frscj(6)
!
!.....ii=205  pellet atoms vs time
      ymin = globmin(205)
      ymax = globmax(205)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.700_R8,1.0_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(206,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,50._R8,0,0,2,1)
      call crtbcd(nlab(205))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
!
!.....ii=206  total atoms
      ymin = globmin(206)
      ymax = globmax(206)
      diffy = ymax-ymin-1.E-12_R8*(abs(ymax)+abs(ymin))
      if(diffy.le.0) ymax = ymin + max(1._R8,abs(ymax),abs(ymin))
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call disk(207,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,18._R8,0,0,2,1)
      call crtbcd(nlab(206))
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(nsc1,7201)
 7201 format("pellet atoms and total atoms")
      call frscj(6)
      go to 400
  811 continue
      write(s100,1811)
      call gtextm(s100,80,0,1,3)
      call disk(76,1,big1,big2)
      call tracec(1hv,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call disk(77,1,big1,big2)
      call tracec(1h0,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call frscj(6)
      go to 400
 1811 format("  *...ave",/,                                              &  
     &       "  v...vac",/,                                              &  
     &       "  0..cent*(1/2)")
  950 continue
      if(ibalsw.le.0) go to 400
      index = ii - pglobs
      if(index.gt.npsi) go to 400
      if(index.gt.1) go to 960
      smin = 0._R8
      smax = max( AREAL(npsit-1), AREAL(npsitmx) )
      call mapg(tmin,tmaxx,smin,smax,.142_R8,.800_R8,.290_R8,1._R8)
      call setold(tmaxx,smin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setld(2._R8,22._R8,0,0,2,1)
      write(s100,1950)
      call gtext(s100,80,0)
 1950 format(" balloon stability")
  960 continue
      aindx = AREAL(index)
      do 970 i=1,nrecord
      if(big1(i).eq.5.0_R8) call pointc(1h-,big2(i),aindx,1,-1,-1,0._R8,  &  
     & 0._R8)
      if(big1(i).eq.6.0_R8) call pointc(1h1,big2(i),aindx,1,-1,-1,0._R8,  &  
     & 0._R8)
      if(big1(i).eq.7.0_R8) call pointc(1h2,big2(i),aindx,1,-1,-1,0._R8,  &  
     & 0._R8)
      if(big1(i).eq.8.0_R8) call pointc(1h3,big2(i),aindx,1,-1,-1,0._R8,  &  
     & 0._R8)
      if(big1(i).eq.3.0_R8) call pointc(1h0,big2(i),aindx,1,-1,-1,0._R8,  &  
     & 0._R8)
      if(big1(i).eq.2.0_R8) call pointc(1h.,big2(i),aindx,1,-1,-1,0._R8,  &  
     & 0._R8)
  970 if(big1(i).eq.1.0_R8) call pointc(1h*,big2(i),aindx,1,-1,-1,0._R8,  &  
     & 0._R8)
      if(index.lt.npsi) go to 400
      write(nsc1,1970)
 1970 format(" balloon stability vs time ")
      call frscj(6)
      go to 400
!
  500 continue
      iimpgls = ii - pgls
      index = mod(iimpgls-1,35) + 1
!..rxw:
      if(noplot(22).gt.0) goto 400
      if(ii.gt.pgls+ncoil) go to 400
      if(ii.eq.pgls+ncoil) go to 521
      if(ii.eq.pgls+1) go to 510
      call tracec(lab2(index),big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)    
      go to 400
  510 continue
      call mapg(tmin,tmaxx,curmin,curmax,.142_R8,.800_R8,.290_R8,1._R8)
      call setold(tmaxx,curmin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call tracec(lab2(index),big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)    
      call setld(2._R8,22._R8,0,0,2,1)
      write(s100,1510)
      call gtext(s100,80,0)
 1510 format(" currents(ka)" )
      go to 400
  521 continue
      call tracec(lab2(index),big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)    
      call setold(tmaxx,ene2max,1,0,1,0)
      write(nsc1,1510)
      call frscj(6)
      go to 400
  400 continue
!
!.....special divertor plots
      if(iplate.eq.0 .or. nplate.eq.0) go to 449
!..rxw:
      if(noplot(14).gt.0) goto 449
!
      ii = pglobs+ppsi
      do 440 n=1,nplate
      ii = ii+2
      ymax = max(globmax(ii-1),globmax(ii))
      call disk(ii+1,1,big1,big2)
      call disk(ii,1,big3,big2)
      ymin = 1.E30_R8
      do 2440 i=1,nrecord
      if(big1(i).ne.0) ymin = min(ymin,big1(i))
      if(big3(i).ne.0) ymin = min(ymin,big3(i))
 2440 continue
      if(ymax.le.ymin) go to 539
      call mapg(tmin,tmaxx,ymin,ymax,                                    &  
     &          .142_R8,.800_R8,.700_R8,1._R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      call setld(2._R8,47._R8,0,0,2,1)
      write(s100,1440) n
      call gtext(s100,80,0)
 1440 format("strike pt plate",i3)
      call tracec(1h1,big2(1),big3(1),nrecord,-1,-1,0._R8,0._R8)
      call tracec(1h2,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
  539 continue
      ii = ii+2
      ymax = max(globmax(ii-1),globmax(ii))
      ymin = min(globmin(ii-1),globmin(ii))
      if(ymax.le.ymin) go to 439
      call mapg(tmin,tmaxx,ymin,ymax,                                    &  
     &          .142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      call setld(2._R8,17._R8,0,0,2,1)
      write(s100,1540) n
      call gtext(s100,80,0)
 1540 format("distance - x-pt to plate",i3)
      call disk(ii,1,big3,big2)
      call tracec(1h1,big2(1),big3(1),nrecord,-1,-1,0._R8,0._R8)
      call disk(ii+1,1,big1,big2)
      call tracec(1h2,big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      write(nsc1,1435) n
 1435 format(" divertor plate",i3,"  strike pt")
      call frscj(6)
  439 continue
      ymin = 0._R8
      ymax = 0._R8
      if(nseg(n).le.0) go to 435
      do 436 i=ii+1,ii+nseg(n)
      ymin = min(ymin,globmin(i))
  436 ymax = max(ymax,globmax(i))
      if(ymax.le.ymin) go to 435
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,0.95_R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      call setld(2._R8,17._R8,0,0,2,1)
      write(s100,1436) n
      call gtext(s100,80,0)
 1436 format(" heat flux..plate",i3)
      xl = xlplate(n)
      xr = xrplate(n)
      zl = zlplate(n)
      zr = zrplate(n)
      dis = sqrt((xl-xr)**2 + (zl-zr)**2)
      disncm = 100._R8*dis/nseg(n)
      ymaxp = ymax + .05_R8*(ymax-ymin)
      call setold(tmin,ymaxp,1,0,1,0)
      write(s100,1437) xl,zl,xr,zr,nseg(n),disncm
      call gtext(s100,80,0)
 1437 format("(",f5.2,",",f5.2,") to (",f5.2,",",f5.2,")  ",             &  
     &   i3," segs of len",f6.1," cm")
      llab = 0
      do 432 i=ii+1,ii+nseg(n)
      call disk(i+1,1,big1,big2)
      llab = llab + 1
      call tracec(lab2(mod(llab-1,35)+1),big2(1),big1(1),nrecord,        &  
     &          -1,-1,0._R8,0._R8)
  432 continue
  435 continue
      write(nsc1,1535) n
 1535 format(" divertor plate",i3,"   heat flux")
      call frscj(6)
      ii = ii + nseg(n)
  440 continue
  449 continue
!
!.....special current group plots
!
      if(ncoil.le.0) go to 592
!..rxw:
      if(noplot(15).gt.0) goto 592
!
      do 591 iig=1,ngroupt
      ig=nogroupt(iig)
      if(gigpmax(ig) .le. gigpmin(ig)) go to 591
      call mapg(tmin,tmaxx,gigpmin(ig),gigpmax(ig),                      &  
     &          .142_R8,.800_R8,.700_R8,1._R8)
      call setold(tmaxx,gigpmin(ig),1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setold(tmaxx,gigpmax(ig),1,0,1,0)
      write(s100,1023) gigpmin(ig),gigpmax(ig)
      call gtextm(s100,80,0,1,3)
      call setld(2._R8,47._R8,0,0,2,1)
      write(s100,1490) ig
      call gtext(s100,80,0)
      llab = 0
      do 490 ii=pgls+1,pgls+ncoil
      il = ii-pgls
      if(iabs(igroupc(il)).ne.ig) go to 490
      call disk(ii+1,1,big1,big2)
      llab = llab + 1
      call tracec(lab2(mod(llab-1,35)+1),big2(1),big1(1),nrecord,        &  
     &          -1,-1,0._R8,0._R8)
  490 continue
 1490 format("group",i3," curr(ka)")
!
!.....special voltage group plots
!
      if(ncoil.le.0) go to 491
      if(gvgpmax(ig) .le. gvgpmin(ig)) go to 491
      call mapg(tmin,tmaxx,gvgpmin(ig),gvgpmax(ig),                      &  
     &          .142_R8,.800_R8,.290_R8,.584_R8)
      call setold(tmaxx,gvgpmin(ig),1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setold(tmaxx,gvgpmax(ig),1,0,1,0)
      write(s100,1023) gvgpmin(ig),gvgpmax(ig)
      call gtextm(s100,80,0,1,3)
      call setld(2._R8,17._R8,0,0,2,1)
      write(s100,1590) ig
      call gtext(s100,80,0)
      llab = 0
      do 590 ii=pgls+ncoil+1,pgls+2*ncoil
      il = ii-(pgls+ncoil)
      if(iabs(igroupc(il)).ne.ig) go to 590
      call disk(ii+1,1,big1,big2)
      llab = llab + 1
      call tracec(lab2(mod(llab-1,35)+1),big2(1),big1(1),nrecord,        &  
     &          -1,-1,0._R8,0._R8)
  590 continue
  491 continue
      write(nsc1,1591) ig
      call frscj(6)
 1590 format("group",i3," volt(kv)")
 1591 format(" group",i3," current and voltage")
!
!.....codeing here to plot gap voltages and currents
      if(iseries(ig).eq.0) go to 591
      izero = pgls + 3*pncoil + 6*pngroup + 6
      i1 = izero + (ig-1)*3 + 1
      i2 = izero + (ig-1)*3 + 2
      i3 = izero + (ig-1)*3 + 3
      call disk(i2+1,1,big1,big2)
      gapmin = 0._R8
      gapmax = 0._R8
      do 691 i=1,nrecord
      gapmin = min(gapmin,big1(i))
      gapmax = max(gapmax,big1(i))
  691 continue
      if(gapmin .ge. gapmax) go to 591
      call mapg(tmin,tmaxx,gapmin,gapmax,.142_R8,.800_R8,.700_R8,1._R8)
      call setold(tmaxx,gapmin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setold(tmaxx,gapmax,1,0,1,0)
      write(s100,1023) gapmin,gapmax
      call gtextm(s100,80,0,1,3)
      call setld(2._R8,45._R8,0,0,2,1)
      write(s100,1691) ig
      call gtext(s100,80,0)
 1691 format("group",i3," gap cur(ka)")
      call trace(big2,big1,nrecord,                                      &  
     &          -1,-1,0._R8,0._R8)
!
      call disk(i3+1,1,big1,big2)
      gapmin = 0._R8
      gapmax = 0._R8
      do 631 i=1,nrecord
      gapmin = min(gapmin,big1(i))
      gapmax = max(gapmax,big1(i))
  631 continue
      if(gapmin .ge. gapmax) gapmax = gapmin + 1._R8
      call mapg(tmin,tmaxx,gapmin,gapmax,.142_R8,.800_R8,.290_R8,        &  
     & .584_R8)
      call setold(tmaxx,gapmin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setold(tmaxx,gapmax,1,0,1,0)
      write(s100,1023) gapmin,gapmax
      call gtextm(s100,80,0,1,3)
      call setld(2._R8,17._R8,0,0,2,1)
      write(s100,1692) ig
      call gtext(s100,80,0)
 1692 format(" group",i3," gap volt(kv)")
      call trace(big2,big1,nrecord,                                      &  
     &           -1,-1,0._R8,0._R8)
!
      write(nsc1,1693) ig
 1693 format(" group",i3," gap cur and volt")
      call frscj(6)
!
  591 continue
  592 continue
!
      izero = pgls + 3*pncoil
      if(gcmax.le.gcmin) go to 792
!..rxw:
      if(noplot(16).gt.0) goto 792
      gcmin = min(gcmin,acoef(703))
      gcmax = max(gcmax,acoef(704))
!
      call mapg(tmin,tmaxx,gcmin,gcmax,.142_R8,.800_R8,.290_R8,1._R8)
      call setold(tmaxx,gcmin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setld(2._R8,20._R8,0,0,2,1)
      write(s100,1793)
      call gtext(s100,80,0)
 1793 format("group curr(kA-T)..actual(solid),prepro(dot)")
!
!.....special printed output requested by G.Sheffield
      write(nout,2690)
 2690 format("1   * * *  summary of group currents(kA)  * * *",//,       &  
     &   "        time    group(1)    group(2)    group(3)    group(4)",  &  
     &                                                                   &  
     &   "    group(5)    group(6)    group(7)    group(8)" )
      iadscr = 0
      do 681 irec=1,nrecord
      do 671 l=1,lenscr
  671 temp(l) = pltsav(iadscr+l)
      iadscr = iadscr + lenscr
!
!.....remultiply group 1 currents by 10
      temp(izero+2) = temp(izero+2)*10._R8
      write(nout,2681) temp(1),(temp(ii),ii=izero+2,izero+16,2)
  681 continue
 2681 format(1p9e12.4)
!


      open(399,file='grp_curr.out',status='unknown',iostat=ios)
      iadscr = 0
      do 3681 irec=1,nrecord
      do 3671 l=1,lenscr
 3671 temp(l) = pltsav(iadscr+l)
      iadscr = iadscr + lenscr
!
!.....remultiply group 1 currents by 10
!      temp(izero+2) = temp(izero+2)*10._R8
      if (l .le.49) &
     & write(399,3399) temp(1),(temp(ii),ii=izero+2,izero+2*pngroup+1,2)
 3681 continue
 3399 format(1p50e12.4)
      close(399)




      llab = 0
      do 690 ii=izero+1,izero+2*pngroup,2
      llab = llab + 1
      call disk(ii+1,ii+2,big1,big3)
      sum = 0._R8
      do 685 i=1,nrecord
  685 sum = sum + big1(i)
      if(sum.eq.0) go to 690
!
!...only plot primary groups
      do 9912 nn=1,ncoil
      l = iabs(igroupc(nn))
      if(l.eq.llab) go to 9913
 9912 continue
      go to 690
 9913 continue
!
      call tracec(lab2(mod(llab-1,35)+1),big2(1),big1(1),nrecord,        &  
     &           -1,-1,0.0_R8,0.0_R8)
  690 continue
      call setold(tmaxx,gcmax,1,0,1,0)
      write(s100,1023) gcmin,gcmax
      call gtextm(s100,80,0,1,3)
      write(s100,2025)
      call gtextm(s100,80,0,1,4)
 2025 format(//," ... preprogrammed",/,                                  &  
     &          " --- actual ")
      write(s100,2026)
      call gtextm(s100,80,0,1,4)
 2026 format(//,"NOTE: Group 1 ",/,                                      &  
     &          " divided by 10 ")
!
      write(nsc1,1690)
 1690 format("current groups,actual & prepro")
      call frscj(6)
!
!     special plots of feedback groups
      do 9940 ii=1,numfb
      ig = nrfb(ii)
      if(ig.eq.0) go to 9940
!
!     first determine min and max
      if = izero + 2*ig - 1
      call disk(1,if+1,big2,big1)
      fbmin = 1.E6_R8
      fbmax =-1.E6_R8
      do 2685 i=1,nrecord
      fbmin = min(fbmin,big1(i))
      fbmax = max(fbmax,big1(i))
 2685 continue
      fbdiff = fbmax - fbmin
!     write(nterm,3333) ig,if,fbmin,fbmax
!3333 format(" ig,if,fbmin,fbmax =",2i5,1p2e12.4)
      if(fbdiff.le.0) go to 9940
      amax = fbmax + .05_R8*fbdiff
      amin = fbmin - .05_R8*fbdiff
      call mapg(tmin,tmaxx,amin,amax,.142_R8,.800_R8,.290_R8,1._R8)
      call setold(tmaxx,amin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setld(2._R8,20._R8,0,0,2,1)
      write(s100,1993) ig
 1993 format(" feedback + group current, group",i3)
      call trace(big2,big1,nrecord,-1,-1,0.0_R8,0.0_R8)
      call setold(tmaxx,gcmax,1,0,1,0)
      write(s100,1023) fbmin,fbmax
      call gtextm(s100,80,0,1,3)
      call setld(2._R8,20._R8,0,0,2,1)
      write(s100,9930) ig
      call gtext(s100,80,0)
      write(nsc1,9930) ig
 9930 format(" feedback currents, group",i3)
      call frscj(6)
 9940 continue
  792 continue
!
      izero = pgls + 3*pncoil + 2*pngroup
      ymax = max(gvmax,gvmax0)
      ymin = min(gvmin,gvmin0)
      if(ymax.le.ymin) go to 842
!..rxw:
      if(noplot(17).gt.0) goto 842
!
      call mapg(tmin,tmaxx,ymin,ymax,.142_R8,.800_R8,.290_R8,1._R8)
      call setold(tmaxx,ymin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setld(2._R8,20._R8,0,0,2,1)
      write(s100,1843)
      call gtext(s100,80,0)
 1843 format("group volt(kV/T)..actual(solid),prepro(dot)")
      llab = 0._R8
!
!.....special printed output requested by G.Sheffield
      write(nout,2790)
 2790 format("1   * * *  summary of group volt/turn(kV)  * * *",//,      &  
     &   "        time    group(1)    group(2)    group(3)    group(4)",  &  
     &                                                                   &  
     &   "    group(5)    group(6)    group(7)    group(8)" )
      iadscr = 0
      do 682 irec=1,nrecord
      do 672 l=1,lenscr
  672 temp(l) = pltsav(iadscr+l)
      iadscr = iadscr + lenscr
      write(nout,2681) temp(1),(temp(ii),ii=izero+2,izero+9)
  682 continue
!
      do 844 ii=izero+1,izero+pngroup
      llab = llab + 1
      call disk(1,ii+1,big2,big1)
      sum = 0._R8
      do 845 i=1,nrecord
  845 sum = sum + big1(i)
      if(sum.eq.0) go to 844
      call tracec(lab2(mod(llab-1,35)+1),big2(1),big1(1),nrecord,        &  
     &           -1,-1,0._R8,0._R8)
  844 continue
      llab = 0
      do 846 ii=izero+pngroup+1,izero+2*pngroup
      llab = llab + 1
      call disk(1,ii+1,big2,big1)
      sum = 0._R8
      do 847 i=1,nrecord
  847 sum = sum + big1(i)
      if(sum.eq.0) go to 846
      call pointc(lab2(mod(llab-1,35)+1),big2,big1,nrecord,              &  
     &           -1,-1,0.0_R8,0.0_R8)
  846 continue
      call setold(tmaxx,ymax,1,0,1,0)
      write(s100,1023) ymin,ymax
      call gtextm(s100,80,0,1,3)
      write(s100,2025)
      call gtextm(s100,80,0,1,4)
!
      write(nsc1,1844)
 1844 format("group voltage..actual & prepro")
      call frscj(6)
  842 continue
!
      izero = pgls + 3*pncoil + 4*pngroup
      if(gpmax.le.gpmin) go to 942
!..rxw:
      if(noplot(18).gt.0) goto 942
!
      call mapg(tmin,tmaxx,gpmin,gpmax,.142_R8,.800_R8,.290_R8,1._R8)
      call setold(tmaxx,gpmin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setld(2._R8,20._R8,0,0,2,1)
      write(s100,1943)
      call gtext(s100,80,0)
 1943 format("group power(Mw)")
      llab = 0._R8
      do 944 ii=izero+1,izero+pngroup
      llab = llab + 1
      call disk(1,ii+1,big2,big1)
      sum = 0._R8
      do 945 i=1,nrecord
  945 sum = sum + big1(i)
      if(sum.eq.0) go to 944
      call tracec(lab2(mod(llab-1,35)+1),big2(1),big1(1),nrecord,        &  
     &           -1,-1,0._R8,0._R8)
  944 continue
      write(nsc1,1944)
 1944 format("group power..actual")
      call frscj(6)
  942 continue
!
      izero = pgls + 3*pncoil + 5*pngroup
      if(gemax.le.gemin) go to 862
!..rxw:
      if(noplot(19).gt.0) goto 862
!
      write(nout,3901)
 3901 format(////," final energy in each coil group (MW) ")
      call mapg(tmin,tmaxx,gemin,gemax,.142_R8,.800_R8,.290_R8,1._R8)
      call setold(tmaxx,gemin,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setld(2._R8,20._R8,0,0,2,1)
      write(s100,1863)
      call gtext(s100,80,0)
 1863 format("group v-sec(W)")
      llab = 0._R8
      sum2 = 0._R8
      do 864 ii=izero+1,izero+pngroup
      llab = llab + 1
      call disk(1,ii+1,big2,big1)
      sum = 0._R8
      do 865 i=1,nrecord
  865 sum = sum + big1(i)
      if(sum.eq.0) go to 864
      call tracec(lab2(mod(llab-1,35)+1),big2(1),big1(1),nrecord,        &  
     &           -1,-1,0._R8,0._R8)
      write(nout,3902) llab,big1(nrecord)
 3902 format(" group ",i3, "  has energy",1pe12.4)
      sum2 = sum2 + big1(nrecord)
  864 continue
      write(nout,3903) sum2
 3903 format(/," total PF energy = ",1pe12.4)
      write(nsc1,1864)
 1864 format("group v-sec...actual")
      call frscj(6)
  862 continue
!
      if(gpmint .ge. gpmaxt) go to 872
!..rxw:
      if(noplot(20).gt.0) goto 872
!
      call mapg(tmin,tmaxx,gpmint,gpmaxt,.142_R8,.800_R8,.700_R8,1.0_R8)    
      call setold(tmaxx,gpmint,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      ii = pgls + 3*pncoil + 6*pngroup + 1
      call disk(ii+1,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      ii = pgls + 3*pncoil + 6*pngroup + 3
      call disk(ii+1,1,big1,big2)
      call tracec(1hr,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      ii = pgls + 3*pncoil + 6*pngroup + 4
      call disk(ii+1,1,big1,big2)
      call tracec(1hi,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,42._R8,0,0,2,1)
      write(s100,1872)
      call gtext(s100,80,0)
 1872 format(" total power (MW)" )
!
      if(gemint .ge. gemaxt) go to 882
      call mapg(tmin,tmaxx,gemint,gemaxt,.142_R8,.800_R8,.290_R8,        &  
     & .584_R8)
      call setold(tmaxx,gemint,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      ii = pgls + 3*pncoil + 6*pngroup + 2
      call disk(ii+1,1,big1,big2)
      call tracec(1h*,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      ii = pgls + 3*pncoil + 6*pngroup + 5
      call disk(ii+1,1,big1,big2)
      call tracec(1hi,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      ii = pgls + 3*pncoil + 6*pngroup + 6
      call disk(ii+1,1,big1,big2)
      call tracec(1hr,big2,big1,nrecord,-1,-1,0._R8,0._R8)
      call setld(2._R8,22._R8,0,0,2,1)
      write(s100,1892)
      call gtext(s100,80,0)
 1892 format(" total energy (Mj)" )
  882 continue
      write(nsc1,2882)
 2882 format(" total power and energy")
      call frscj(6)
  872 continue
!
      if(itemp.le.0) go to 692
!..rxw:
      if(noplot(21).gt.0) goto 692
!
      izero = pgls + 3*pncoil + 9*pngroup + 6
      if(templmx.le.templmn) go to 852
      call mapg(tmin,tmaxx,templmn,templmx,.142_R8,.800_R8,.290_R8,      &  
     & 1._R8)
      call setold(tmaxx,templmn,1,0,1,0)
      write(s100,6666)
      call gtext(s100,80,0)
      call setld(2._R8,20._R8,0,0,2,1)
      write(s100,1853)
      call gtext(s100,80,0)
 1853 format("coil temps (degrees k) ")
      llab = 0
      do 854 ii=izero+1,izero+numpf
      llab = llab + 1
      call disk(1,ii+1,big2,big1)
      sum = 0._R8
      do 855 i=1,nrecord
  855 sum = sum + big1(i)
      if(sum.eq.0) go to 854
      call tracec(lab2(mod(llab-1,35)+1),big2(1),big1(1),nrecord,        &  
     &           -1,-1,0._R8,0._R8)
  854 continue
      write(nsc1,1854)
 1854 format(" coil temps (degrees k) ")
      call frscj(6)
  852 continue
!
  692 continue
!
!
!.....plot observation point time historys
!
!..rxw:
      if(noplot(13).le.0)                                                &  
     & call plotobs(tmin,tmaxx,big1,big2,big3,big4,big5)
!
!
!......special output file added 08/11/2011 for Dakota
      maxdakoda  = 5       ! number of output scalars
      allocate(idakoda(maxdakoda),adakoda(maxdakoda),ndakoda(maxdakoda))
!
      idakoda(1) = 218     ! Te_ped
      idakoda(2) = 219     ! Ti_ped
      idakoda(3) = 213     ! Sum of NI current drive
      idakoda(4) = 8       ! Plasma Current
      idakoda(5) = 140     ! H_98
!
      ndakoda(1) = "  Te_ped  "
      ndakoda(2) = "  Ti_ped  "
      ndakoda(3) = "  I_NI    "
      ndakoda(4) = "  I_P     "
      ndakoda(5) = "  H_98    "
!
      do ii=1,maxdakoda
        adakoda(ii) = pltsav(lenscr*(nrecord-1) + idakoda(ii) + 1)
      enddo
      open(299,file='dakota.out',status='unknown',iostat=ios)
        write(299,2998) (ndakoda(ii),ii=1,maxdakoda)
        write(299,2999) (adakoda(ii),ii=1,maxdakoda)
      close(299)
      return
 2998 format(10(a10,2x))
 2999 format(1p10e12.4)
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
