!#include "f77_dcomplx.h"
      subroutine advimp1d
      USE CLINAM
      USE RADTAB
      USE SAPROP
      USE SCR1
      USE SPECIE
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,jp,imppelp,jsrco,jsrcn,in,jscrn,jdiro,jdir,jmax
      INTEGER jmin,nchrgs,ispcr,nimp,jq,nchrgs1,i,is,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 donor,pelabl,dndt,dmdt,tej,plasvol,timpel,fracpel
      REAL*8 sumane0,zpelo,xpelo,dum1,dpsidxo,dpsidzo,dum4,pspelo
      REAL*8 dum5,dum6,dum7,dpsidxn,dpsidzn,pspeln,ends,zone,at
      REAL*8 eval1,epval,rval1,rpval,eval2,rval2,anej,trmpel,trmr
      REAL*8 radpelo,trmrmx,befo,denom,deltpsio,deltpsin,psi3,psi2
      REAL*8 psi4,psi1,fractot,psiza,psizb,fracj,fac,fcc,facvol
      REAL*8 sraves,srcimps,sumanen,delane,rlcpsio,flximp,dfacp
      REAL*8 vfacp,dfacm,vfacm,vtrm1,vtrm2,vtrm3,vtrm4,facn,tekev
      REAL*8 ainzrjq,radrjq,recrjq,sumnq,pimlt,sumpe,pltcmx,pltcmn
      REAL*8 termplot,pltcmnn,pltcmxx,dpltc,plttmn,plttmx,dpltt
      REAL*8 sumrad
      REAL*8 AREAL, sum
!============
!     dimension snp(ppsi),aimp(ppsi),bimp(ppsi),cimp(ppsi),
!    .   dimp(ppsi),ework(ppsi),fwork(ppsi),nqt(pchrgmx,pimp,ppsi),
!    .   dcstoc(ppsi), vpinch(ppsi), totden(ppsi,pimp),
!    .   srcimp(ppsi), tesave(ppsi), anesave(ppsi),
!    .      ajtary(ppsi),dsrave(ppsi),dsrcimp(ppsi)
      character*1 lab2(40), lab3(7)
      data lab2 / "0","1","2","3","4","5","6","7","8","9"                &  
     &,           "a","b","c","d","e","f","g","h","i","j"                &  
     &,           "k","l","m","n","o","p","q","r","s","t"                &  
     &,           "u","v","w","z","y","z","#","$","%","^"/
      data lab3   / "o","c", "i", "b","n","k","a"/
      data donor  / 0.0_R8/
      data pelabl /0._R8/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: snp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aimp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: bimp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cimp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dimp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ework
      REAL*8, ALLOCATABLE, DIMENSION(:) :: fwork
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nqt
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dcstoc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: vpinch
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: totden
      REAL*8, ALLOCATABLE, DIMENSION(:) :: srcimp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tesave
      REAL*8, ALLOCATABLE, DIMENSION(:) :: anesave
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ajtary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dsrave
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dsrcimp
!============      
      IF(.not.ALLOCATED(snp)) ALLOCATE( snp(ppsi), STAT=istat)
      IF(.not.ALLOCATED(aimp)) ALLOCATE( aimp(ppsi), STAT=istat)
      IF(.not.ALLOCATED(bimp)) ALLOCATE( bimp(ppsi), STAT=istat)
      IF(.not.ALLOCATED(cimp)) ALLOCATE( cimp(ppsi), STAT=istat)
      IF(.not.ALLOCATED(dimp)) ALLOCATE( dimp(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ework)) ALLOCATE( ework(ppsi), STAT=istat)
      IF(.not.ALLOCATED(fwork)) ALLOCATE( fwork(ppsi), STAT=istat)
      IF(.not.ALLOCATED(nqt)) ALLOCATE( nqt(pchrgmx,pimp,ppsi),          &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(dcstoc)) ALLOCATE( dcstoc(ppsi), STAT=istat)
      IF(.not.ALLOCATED(vpinch)) ALLOCATE( vpinch(ppsi), STAT=istat)
      IF(.not.ALLOCATED(totden)) ALLOCATE( totden(ppsi,pimp),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(srcimp)) ALLOCATE( srcimp(ppsi), STAT=istat)
      IF(.not.ALLOCATED(tesave)) ALLOCATE( tesave(ppsi), STAT=istat)
      IF(.not.ALLOCATED(anesave)) ALLOCATE( anesave(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ajtary)) ALLOCATE( ajtary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(dsrave)) ALLOCATE( dsrave(ppsi), STAT=istat)
      IF(.not.ALLOCATED(dsrcimp)) ALLOCATE( dsrcimp(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : advimp1d  ' 
!============      
!.....................................................................
!     PART 1:   Pellet Injection
!.....................................................................
!
!.....initialize source parameters
      dndt = 0.0_R8
      dmdt = 0.0_R8
      tej = tevv
      do 1 j = 1, npsi
!.....check for contribution due to density jet
      if(acoef(816).gt.0._R8.and.                                        &  
          (times.ge.acoef(817) .and. times.le.acoef(818))) then
          srcimp(j) = sravejet(j)*acoef(819)/(1._R8-acoef(819))*udsd/udst
      else
          srcimp(j) = 0.0_R8
      endif
         srave (j) = 0.0_R8
         dsrave(j) = 0.0_R8
         dsrcimp(j) = 0.0_R8
    1 continue
         plasvol = 0._R8
      do 12 jp=2,npsit
   12 plasvol = plasvol + vp(jp)
!.....initialize pellet parameters when pellet is injected
!
      if(acoef(760) .le. 0) go to 101
      if( times      .ge. acoef(766) .and.                               &  
        ((times-dts) .lt. acoef(766) .or.  kcycle.eq.0))  then
!
           vxpel   = acoef(761)
           vzpel   = acoef(762)
            xpel   = acoef(763)
            zpel   = acoef(764)
           radpel  = acoef(765)
           timpel  = acoef(766)
           fracpel = acoef(767)
!
           sumane0 = 0._R8
           do 8888 j=2,npsit
           sumane0 = sumane0 + vp(j)*dpsi*ane(j)
 8888      continue
      endif
!
!.....check for second pellet and initialize if necessary
!
      if(acoef(770) .ge. 1._R8.and.                                      &  
     &   (times.ge.acoef(776) .and.                                      &  
     &  ((times-dts) .lt. acoef(776) .or. kcycle.eq.0))) then
!
           vxpel   = acoef(771)
           vzpel   = acoef(772)
            xpel   = acoef(773)
            zpel   = acoef(774)
          radpel   = acoef(775)
          timpel   = acoef(776)
         fracpel   = acoef(777)
!
!.....keep shooting in pellets every acoef(778) seconds until t=acoef(779)
           if(times .lt. acoef(779) - acoef(778)) then
                          acoef(776) = acoef(776) + acoef(778)
           endif
!
      endif
!
      imppelp = imppel
      if(fracpel.gt.0 .and. fracpel.lt.1.0_R8) imppelp = 0
!.....NOTE:   imppelp=0 indicates pure D pellet
!
! * * (1.1) transport pellet
!
      jsrco = npsit+1
      jsrcn = npsit+1
!
      if(dts.le.0 .or. radpel .le. 1.E-6_R8                              &  
     &            .or. times .lt. timpel) go to 101
      zpelo = zpel
      xpelo = xpel
      call grap(1,zpelo,xpelo,dum1,dpsidxo,dpsidzo,dum4,pspelo           &  
     &,                     dum5,dum6,dum7,0)
      xpel  = xpelo + vxpel*dts
      zpel  = zpelo + vzpel*dts
      call grap(1,zpel, xpel, dum1,dpsidxn,dpsidzn,dum4,pspeln           &  
     &,                     dum5,dum6,dum7,0)
      if(pspeln.gt.psilim .or. xpel.lt.ccon .or. xpel.gt.alx             &  
     &                    .or. zpel.lt.-alz .or. zpel.gt.alz) go to 101
      do 2 j = 2, npsit
        if(pspelo.ge. xsv2(j-1) .and. pspelo .lt. xsv2(j)) then
           jsrco = j
        endif
        if(pspeln.ge. xsv2(j-1) .and. pspeln .lt. xsv2(j)) then
           jsrcn = j
        endif
    2 continue
      if(jsrco .gt. npsit .or. jsrcn .gt. npsit) go to 101
!
!.....NOTE:  Pellet starts in zone jsrco at psi = pspelo
!                     ends in zone jscrn at psi = pspeln
!
!.....evaluate temperature and density at pellet end point
      call eeval(pspeln,2,eval1,epval,imag,jmag)
      call reval(pspeln,idens,1,rval1,rpval,imag,jmag)
      call eeval(pspelo,2,eval2,epval,imag,jmag)
      call reval(pspelo,idens,1,rval2,rpval,imag,jmag)
      anej= .5_R8* (rval1+rval2)*udsd
      tej = .5_R8*((eval1+eval2)/anej)*udsh
!
!.....call pellet ablation routine, calculate density source, and shrink pellet
      call pelkut(radpel, tej, anej, imppelp, trmpel, dndt)
!.....NOTE:   modify ablation rate by factor acoef(809)
      dmdt = dndt * acoef(809)
      trmr    =-dmdt*dts/trmpel
      radpelo = radpel
      trmrmx  =-radpelo**3
      if(abs(trmr).gt.abs(trmrmx)) trmr = trmrmx
      radpel    = (max(1.E-24_R8,radpelo**3 + trmr))**(1._R8/3._R8)
!
!.....distribute source amongst several zones if necessary
      jdiro = jdir
      jdir = sign(1._R8,pspeln-pspelo)
!
      befo = 0._R8
      denom = sqrt(vxpel**2+vzpel**2)
      if(denom.gt.0) befo = radpel/denom
      deltpsio = abs((vxpel*dpsidxo + vzpel*dpsidzo)*befo)
      deltpsin = abs((vxpel*dpsidxn + vzpel*dpsidzn)*befo)
      if(jdir.gt.0) then
!        pspeln > pspelo
         psi3 = pspeln
         psi2 = pspelo
         psi4 = pspeln + deltpsin
         psi1 = pspelo - deltpsio
      else
         psi2 = pspeln
         psi3 = pspelo
         psi4 = pspelo + deltpsio
         psi1 = pspeln - deltpsin
      endif
!
!.....distribute finite sized pellet over trajectory
       fractot = 0
      do j=1,npsit-1
        psiza = xsv2(j)
        psizb = xsv2(j+1)
        call pelsource(psi1,psi2,psi3,psi4,psiza,psizb,fracj)
        srcimp(j) = srcimp(j) + dmdt*fracj*fracpel
        srave(j) = (dmdt/(vp(j)*dpsi))*fracj*(1._R8-fracpel)*usdd/usdt
        fractot = fractot + fracj
      enddo
!
            write(6,1111) jsrco,jsrcn,pspeln,tej,anej,dndt, radpel,      &  
     &                    fractot
 1111 format(" j=",2i4," pspeln= ",1pe11.4," te=",1pe11.4,               &  
     & " ne=",1pe11.4," dndt=", 1pe11.4," radpel=",1pe11.4,              &  
     & " fractot=",1pe11.4)
!
!
!
!.....new coding added 10/28/?? to average density source by backaveraging
      if(acoef(790) .ne. 0) then
        fac = acoef(790)
        fcc = 1._R8- fac
        do 33 j=jsrco, jsrcn, jdir
          jmax = npsit
          jmin = j
          if(jdir.gt.0) jmin = 2
          facvol = 0
        do 13 jp = jmax,jmin,-1
   13     facvol = facvol + vp(jp)
        facvol = facvol + acoef(792)*plasvol
        do 11 jp=jmax,jmin,-1
          dsrave(jp)  =  dsrave(jp) + (vp(j)/facvol)*srave(j)
   11     dsrcimp(jp) = dsrcimp(jp) + (vp(jp)/facvol)*srcimp(j)
   33   continue
        do 14 j=2,npsit
          sraves = srave(j)
          srcimps = srcimp(j)
          srave(j)  = fac*dsrave(j)  + fcc*sraves
   14     srcimp(j) = fac*dsrcimp(j) + fcc*srcimps
      endif
!
!
      if(mod(kcycle,20).eq.0)                                            &  
     &write(6,1112) xpel,dmdt,radpel,tepel,anej,sumpel
 1112 format ("xp=",1pe11.4," dmdt=",1pe11.4," radp=",1pe11.4            &  
     &,   " te=",1pe11.4," ane=",1pe11.4," sum=",1pe11.4)
!
  101 continue
!
!.....save some things for printing and ploting
      tepel = tej
      dndtpel = dmdt
!
! * * * check particle balance
!
      sumanen = 0._R8
      do 8889 j=2,npsit
      sumanen = sumanen + vp(j)*dpsi*ane(j)
 8889 continue
      pelabl = pelabl + dmdt*dts
      delane = sumanen-sumane0
! ....impurity radiant emmission
      totradp = 0.0_R8
      do 4 j = 2,npsit
      totradp = totradp + sradion(j)*vp(j)*dpsi
    4 continue
      if(imppel .gt. 0) then
         sumpel = 0._R8
         do 5 j = 2, npsit
         nchrgs = nchrgsr(imppel)
         do 5 ispcr = 1,nchrgs
         sumpel = sumpel + nq(ispcr,imppel,j)*dpsi
    5    continue
      else
         sumpel = pelabl
      endif
!
!....................................................................
!     PART 2:   Diffuse impurity
!....................................................................
!
!
      rlcpsio = acoef(852)/(AREAL(npsit-1)*dpsi)
      do 10 j = 1,npsit
      dcstoc(j) = acoef(851) * fbchi
      if(times.ge.acoef(95) .and. qprof2(j).le.acoef(96))                &  
     &    dcstoc(j) = dcstoc(j)*acoef(124)
      simpe(j) = 0.0_R8
      vpinch(j) = -acoef(851)*rlcpsio                                    &  
     &          * sqrt(gxmja2(j)/xmja2(j))*tpi*qprof2(j)
      if(idens.eq.3) vpinch(j) = -acoef(852)
   10 continue
      vpinch(1) = 0.0_R8
      dcstoc(1) = 0.0_R8
      if(impbnd.eq.0) then
!
!.....zero flux boundary conditions (dcstoc & vpinch = 0)
!
      vpinch(npsit) = 0.0_R8
      dcstoc(npsit) = 0.0_R8
               else
!
!.....diffusive boundary conditions
!
      do 15 nimp = 1,pimp
      nchrgs = nchrgsr(nimp)
      do 15 ispcr = 1,nchrgs
   15 nq(ispcr,nimp,npsit+1) = 0.0_R8
               endif
!
! * * * impurity source    (#/s, cell volume = vp*dpsi)
!
      flximp = acoef(853)
!
      do 500 nimp = 1,pimp
      if(fraci(nimp).le.0.0_R8) go to 500
      nchrgs = nchrgsr(nimp)
      if(iimp.eq.2) then
         do 100 ispcr = 1,nchrgs
!
         dfacp = dcstoc(1)*vp2(1)*dts*rdpsi**2                              &  
                      *(gxmja2(1)/xmja2(1))*(tpi*qprof2(1))**2
         vfacp = dts*rdpsi*vpinch(1)                                        &  
                    * sqrt(gxmja2(1)/xmja2(1))*(tpi*qprof2(1))
!
         do 110 j = 2,npsit
!
! * * * diffusion
!
         dfacm = dfacp
         dfacp = dcstoc(j)*vp2(j)*dts*rdpsi**2                              &  
                      *(gxmja2(j)/xmja2(j))*(tpi*qprof2(j))**2
         aimp(j-1) = -dfacm/vp(j-1)
         bimp(j-1) =  1._R8+ (dfacp+dfacm)/vp(j)
         cimp(j-1) = -dfacp/vp(j+1)
!
! * * pinch velocity
!
         vfacm = vfacp
         vfacp = dts*rdpsi*vpinch(j)                                        &  
            * sqrt(gxmja2(j)/xmja2(j))*(tpi*qprof2(j))
         vtrm1 = 0.5_R8*vfacp*(1._R8-donor*sign(1._R8,vpinch(j)))
         vtrm2 = 0.5_R8*vfacp*(1._R8+donor*sign(1._R8,vpinch(j)))
         vtrm3 = 0.5_R8*vfacm*(1._R8-donor*sign(1._R8,vpinch(j-1)))
         vtrm4 = 0.5_R8*vfacm*(1._R8+donor*sign(1._R8,vpinch(j-1)))
         cimp(j-1) = cimp(j-1) + vtrm1
         bimp(j-1) = bimp(j-1) + vtrm2 - vtrm3
         aimp(j-1) = aimp(j-1) - vtrm4
!
         dimp(j-1) = nq(ispcr,nimp,j)
  110    continue
         if(impbnd.ne.0) then
            dimp(npsit-1)=dimp(npsit-1)-cimp(npsit-1)*nq(ispcr,nimp,npsit+1)
         endif
!
! * * * impurity source    (#/s, cell volume = vp*dpsi)
!
         if(ispcr.eq.1) then
            dimp(npsit-1) = dimp(npsit-1) + dts*flximp/dpsi
               if(nimp.eq.imppel) then
                  jdir = sign(1._R8,pspeln-pspelo)
                  do 111 j = 2,npsit
                  dimp(j-1) = dimp(j-1)+dts*srcimp(j)/dpsi
  111             continue
               endif
         endif
!
         call tridiag(aimp,bimp,cimp,dimp,snp,ework,fwork,npsit-1)
!
         do 120 j = 2,npsit
         nqt(ispcr,nimp,j) = snp(j-1)
  120    continue
  100    continue
      else   !  on iimp
         do 140 j=2,npsit
         if(ilte.eq.0) then
            nqt(1,nimp,j) = facimp*fraci(nimp)*ane(j)*vp(j)
            do 130 ispcr=2,nchrgs
            nqt(ispcr,nimp,j) = 0.0_R8
  130       continue
         else
            do 135 ispcr=1,nchrgs
            facn = 1._R8
            if(adno(j).gt.0.0_R8) facn = adn(j)/adno(j)
            nqt(ispcr,nimp,j) = nq(ispcr,nimp,j)*facn
  135    continue
         endif  !  on ilte
  140 continue
      endif ! on iimp
!..................................................................
!     PART 3:  Impurity ionization
!..................................................................
!
      do 200 j = 2,npsit
!
      do 210 jq = 1,nchrgs
      tekev = te(j)*1.E-3_R8
      call lookup (nimp,jq,tekev,ane(j),ainzrjq,radrjq,                  &  
     &recrjq)
      ainz(jq,j) = ane(j)*ainzrjq
!     if(ainz(jq,j) .lt. 1.e-40) ainz(jq,j) = 1.e-40
!     if(ainz(jq,j) .lt. 1.e-40) then
!
!     write(nterm,4455) jq,j,ainz(jq,j)
!4455 format(" jq,j,ainz =",2i5,e12.5)
!     ainz(jq,j) = 1.e-40
!     endif
      rad(nimp,jq,j) = ane(j)*radrjq
      rec(jq,j) = ane(j)*recrjq
  210 continue
 
      if(ilte.eq.0) then
!
! * * LTE model
!
      sum = 0.0_R8
      nchrgs1 = nchrgs-1
      sumnq = 0.0_R8
      do 220 i = 1,nchrgs
      sumnq = sumnq + nqt(i,nimp,j)
      pimlt = 1.0_R8
      is = i
      if(i.eq.nchrgs) go to 230
      do 240 k = is,nchrgs1
      pimlt = pimlt*rec(k+1,j)/ainz(k,j)
!     if(pimlt.gt.1.e40) pimlt = 1.e40
!     if(pimlt.gt.1.e40) then
!     write(nterm,4466) pimlt
!4466 format(" pimlt =",e12.5)
!     endif
!     endif
  240 continue
  230 sum = sum+pimlt
  220 continue
 
      snp(nchrgs) = sumnq/sum
      nn = nchrgs1
      do 250 i = 1,nchrgs1
      snp(nn) = rec(nn+1,j)*snp(nn+1)/ainz(nn,j)
      nn = nn-1
  250 continue
      else
!
! * * NLTE
!
      aimp(1) = 0.0_R8
      do 310 i = 1,nchrgs
      if(i.gt.1) aimp(i) = -dts*ainz(i-1,j)
      bimp(i) =  1._R8+ dts*(ainz(i,j)+rec(i,j))
      if(i.lt.nchrgs) cimp(i) = -dts*rec(i+1,j)
      dimp(i) = nqt(i,nimp,j)
  310 continue
      cimp(nchrgs) = 0.0_R8
!
      call tridiag(aimp,bimp,cimp,dimp,snp,ework,fwork,nchrgs)
!
      endif
      do 410 i = 1,nchrgs
      if(dts.gt.0.0_R8) then
      simpe(j) = simpe(j) + AREAL(i-1)*(snp(i)-nq(i,nimp,j))             &  
     &                 /(vp(j)*dts*usdt*udsd)
               endif
      nq(i,nimp,j) = snp(i)
  410 continue
      sumpe = sumpe + simpe(j)*dts*usdt*udsd*vp(j)*dpsi
!
  200 continue
!
      do 600 j = npsit+1,npsi
      simpe(j) = 0.0_R8
      nq(1,nimp,j) = facimp*fraci(nimp)*ane(j)*vp(j)
!
      do 610 i = 2,nchrgs
      nq(i,nimp,j) = 0.0_R8
  610 continue
  600 continue
  500 continue
!
!---->special diagnostic printout
!     write(nterm,9333)
!9333 format(" special diagnostic printout of simpe and srave")
!     write(nterm,9331) (srave(j),j=1,npsit)
!     write(nterm,9332) (simpe(j),j=1,npsit)
!9331 format(" srave",1p5e12.4)
!9332 format(" simpe",1p5e12.4)
!     stop
!
      return
!
!...................................................................
!     PART 4:   Impurity Plots
!...................................................................
!
      entry impplot
!
! .... plot 1, impurity charge state density
      pltcmx = -1.E60_R8
      pltcmn =  1.E60_R8
!
      do 700 nimp = 1,pimp
      nchrgs = nchrgsr(nimp)
      if(fraci(nimp).le.0.0_R8) go to 700
      do 705 j = 1,npsit
      totden(j,nimp) = 0.0_R8
  705 continue
!
      do 710 j = 2,npsit
      do 710 ispcr = 1,nchrgs
      totden(j,nimp) = totden(j,nimp) + nq(ispcr,nimp,j)/vp(j)
      termplot = nq(ispcr,nimp,j)/vp(j)
      if(termplot .gt. 0) pltcmn = min(pltcmn,termplot)
  710 continue
!
      totden(1,nimp) = totden(2,nimp)
      call limits(totden(1,nimp),1,npsit,.10_R8,.05_R8,pltcmnn,pltcmxx)
      if(pltcmx.lt.pltcmxx) pltcmx = pltcmxx
      if(pltcmn.gt.pltcmnn) pltcmn = pltcmnn
  700 continue
      if(pltcmx .le.pltcmn) go to 801
!
      pltcmn = max(pltcmn,1.E-6_R8*pltcmx)
      dpltc = + .15_R8*(pltcmx-pltcmn)
      plttmn = xsv2(1)
      plttmx = xsv2(npsit)
      dpltt = + .175_R8*(plttmx-plttmn)
      call mapssl(plttmn,plttmx,pltcmn,pltcmx,.15_R8,.50_R8,.60_R8,      &  
     & .95_R8)
      do 730 nimp = 1,pimp
      nchrgs = nchrgsr(nimp)
      if(fraci(nimp).le.0.0_R8) go to 730
      do 735 ispcr = 1,nchrgs
      do 736 j = 1,npsit+1
  736 aimp(j) = max(pltcmn,nq(ispcr,nimp,j)/vp(j))
      aimp(1) = aimp(2)
      call tracec(lab2(ispcr),xsv2(1),aimp(1),npsit,-1,-1,0._R8,0._R8)
  735 continue
      call map   (plttmn,plttmx,pltcmn,pltcmx,.15_R8,.50_R8,.60_R8,      &  
     & .95_R8)
      call setold(plttmn+dpltt,pltcmx+.2_R8*dpltc,1,0,1,0)
      write(s100,1899)     times
      call gtext(s100,80,0)
 1899 format(                                                            &  
     &" time(s)= ",1pe11.4)
  730 continue
      call setold(plttmn-1.20_R8*dpltt,pltcmn,1,0,1,1)
      write(s100,2008)
      call gtext(s100,80,0)
 2008 format('charge state density (per m3)')
      call setold(plttmn+dpltt,pltcmn-1.1_R8*dpltc,1,0,1,0)
      write(s100,2009)
      call gtext(s100,80,0)
 2009 format('plasma poloidal flux')
!
! .... plot 2, impurity total densities
      call maps(plttmn,plttmx,pltcmn,pltcmx,.65_R8,1.0_R8,.60_R8,.95_R8)    
      do 740 nimp = 1,pimp
      if(fraci(nimp).le.0.0_R8) go to 740
      call tracec(lab3(nimp),xsv2(1),totden(1,nimp),npsit,-1,-1,0._R8,   &  
     & 0._R8)
  740 continue
      call setold(plttmn-1.3_R8*dpltt,pltcmn,1,0,1,1)
      write(s100,2018)
      call gtext(s100,80,0)
 2018 format('impurity total density (per m3)')
      call setold(plttmn+dpltt,pltcmn-1.1_R8*dpltc,1,0,1,0)
      write(s100,2009)
      call gtext(s100,80,0)
  801 continue
!
! .... plot 3, impurity radiant emmission
      sumrad = 0.0_R8
      do 750 j = 2,npsit
      aimp(j) = sradion(j)*1.E-6_R8
      sumrad = sumrad + sradion(j)*vp(j)*dpsi
  750 continue
      aimp(1) = aimp(2)
      call limits(aimp(1),1,npsit,.10_R8,.05_R8,pltcmn,pltcmx)
      if(pltcmx .le. pltcmn) go to 802
      dpltc = + .15_R8*(pltcmx-pltcmn)
      call maps(plttmn,plttmx,pltcmn,pltcmx,.15_R8,.50_R8,.15_R8,.50_R8)    
      call setold(plttmn+dpltt,pltcmx+.2_R8*dpltc,1,0,1,0)
      write(s100,1898)     sumrad* 1.E-6_R8
      call gtext(s100,80,0)
 1898 format(                                                            &  
     &" total (MW)=",1pe11.4)
      call tracec(1h ,xsv2(1),aimp(1),npsit,-1,-1,0._R8,0._R8)
      call setold(plttmn-dpltt,pltcmn,1,0,1,1)
      write(s100,2028)
      call gtext(s100,80,0)
 2028 format('radiant power density (MW/m3)')
      call setold(plttmn+dpltt,pltcmn-1.1_R8*dpltc,1,0,1,0)
      write(s100,2009)
      call gtext(s100,80,0)
  802 continue
!
! ...... plot 4 , pinch velocity
      call limits(vpinch(1),1,npsit,.10_R8,.05_R8,pltcmn,pltcmx)
      if(pltcmx .le. pltcmn) go to 2039
      dpltc = + .15_R8*(pltcmx-pltcmn)
      if(acoef(851).eq.0.0_R8) go to 2039
      call maps(plttmn,plttmx,pltcmn,pltcmx,.65_R8,1.0_R8,.15_R8,.50_R8)    
      call setold(plttmn+dpltt,pltcmx+dpltc,1,0,1,0)
      call tracec(1h ,xsv2(1),vpinch(1),npsit,-1,-1,0._R8,0._R8)
      call setold(plttmn-dpltt,pltcmn+dpltc,1,0,1,1)
      write(s100,2038)
      call gtext(s100,80,0)
 2038 format('impurity pinch velocity (m/s)')
      call setold(plttmn+dpltt,pltcmn-1.1_R8*dpltc,1,0,1,0)
      write(s100,2009)
      call gtext(s100,80,0)
 2039 call frscj(16)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
