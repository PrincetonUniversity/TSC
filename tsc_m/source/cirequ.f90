      subroutine cirequ
!     -----------------
!
!.....solve circuit equations to define ccoil array
!
!.....description of variables used:
!
!     nreg  : total number of external coil groups
!     jnreg(i), i=1,nreg : contains the number of the ext. coil groups
!                           in increasing order
!
!     inductance matrix of ext. coil groups:
!     ameg   (i,j) : i=1,nreg; j=1,nreg
!     ameginv(i,j) : inverse inductance matrix
!
!     resistances of ext. groups (can change in time):
!     rsceg(i) : i=1,nreg
!
!     currents of ext. groups:
!     cegp1(i) : at time point n+1;  i=1,nreg
!     ceg  (i) :               n
!     cegm1(i) :               n-1
!
!     flux on ext. groups due to internal wires and plasma:
!     flxeg  (i) : at time point n  ;  i=1,nreg
!     flxegm1(i) :               n-1;
!     flxegm2(i) :               n-2;
!
!     group voltages:
!     veg  (i) : at time point n ; i=1,nreg
!     vegm1(i) : at time point n-1
!
!     contribution to voltages by ideal ohmic heating system:
!     voheg  (i) : at time point n ; i=1,nreg
!     vohegm1(i) :               n-1
!
!     time steps: dt, dtold
!
      USE CLINAM
      USE COMWOY
      USE SCR6
      USE TCVCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!.. KML 10/25/88
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ii,jg,iic,ic,iw,i,iswtch,ig,kg,npro,iut1,iut2,iswitch
      INTEGER i1,i2,itmp,njg,j,nkg,jj,nstart,nfb,nd
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cpv,cdv,civ,v,rs,rt,flxplas,sumcx,f1,bp,bd,cp,cd,det
      REAL*8 epsminv,esc,rsc,rtc,rf1,rf2,vmin,vmax,ratio,diff,out1
      REAL*8 out2,out3,out4,out5,out6,out7,dtm2,dtm1,sumdia
      REAL*8 sumoff,ans,amutlrc,alin,sumsym,ans1,ans2
      REAL*8 sum
!============
      dimension cpv(5), cdv(5) ,civ(5) ,v(7)
      data rs,rt/.0075_R8,.0625_R8/
!============      
!
!
      if(icirc.eq.1) go to 500
!
!...............................................................................
!..... 1.0  Only for ICIRC=0  Don't solve circuit equations, set ccoil=ccoil0
!...............................................................................
!
!
      if(ncoil.eq.nwire) go to 101
      do 100 ii=1,ncoil-nwire
      ccoil(ii) = ccoil0(ii)
      ccoils(ii) = ccoil(ii)*udsi
  100 continue
  101 continue
!
      return
!
  500 continue
!
!...............................................................................
!..... 2.0 For ICIRC=1, first calculate some group quantities
!...............................................................................
!                flxeg(jg) is total flux at coil group jg due to plasma and
!                          internal coils
!                veg(jg)   is internal work array of units voltage
!                rsceg(jg) is resistance of coil group jg
!                ceg0(jg)  is desired current in coil group jg
!
!
      if(ncoil.eq.nwire) goto 900
      if(itemp.eq.1) call temprise
!
      do 550 jg=1,nreg
        flxeg(jg) = 0._R8
        veg  (jg) = 0._R8
        rsceg(jg) = 0._R8
        iic = 0
!
        do 520 ic=1,ncoil-nwire
          if(igroupc(ic) .ne. jnreg(jg)) goto 520
          iic = iic + 1
          xtmp1(iic) = xcurf
          ztmp1(iic) = zcurf
          xtmp2(iic) = xcoil(ic)
          ztmp2(iic) = zcoil(ic)
          atrn2(iic) = aturnsc(ic)
          veg(jg)    = veg(jg)   + rscoil(ic)*ccoil0(ic)*aturnsc(ic)
          rsceg(jg)  = abs( rsceg(jg) + rscoil(ic)*aturnsc(ic)**2 )
 520    continue
!
!.........sum up flux on ext. coils due to internal wires:
          do 510 iw=1,nwire
            flxeg(jg) = flxeg(jg)                                        &  
     &                + atans(jg,iw) * ccoil(ncoil-nwire+iw)
 510      continue
!
       call gvect(xtmp1,ztmp1,xtmp2,ztmp2,iic,g1,g2,g3,g4,g5,g6,nmult    &  
     &            ,ineg)
        do 530 i=1,iic
          flxplas = g1(i)*tcurdtp                                        &  
     &            + .5_R8*( g4(i)*(cmom(1,2) + cmom(2,1) )               &  
     &            + g5(i)*cmom(2,2)  +  g6(i)*cmom(1,1))
      if(lrswtch.gt.0) flxplas=0
          flxeg(jg) = flxeg(jg) - atrn2(i)*tpi*flxplas
 530    continue
!
        flxeg(jg) = (isym+1) * flxeg(jg)
!
        veg(jg)   = (isym+1) * veg(jg)
        rsceg(jg) = (isym+1) * rsceg(jg)
        ceg0(jg) = veg(jg)/rsceg(jg)
!
 550  continue
!     write(nterm,2550)
!2550 format(1x,"  jg        flxeg     veg       rsceg     ceg0")
!     do 2551 jg=1,nreg
!2551 write(nterm,2552) jg,flxeg(jg),veg(jg),rsceg(jg),ceg0(jg),
!    1  ceg(jg)
!2552 format(i5,5x,1p5e10.2)
!
!..............................................................................
!..... 3.0  For kcycle.le.1 only    initialize some quantities
!..............................................................................
!
!
!.. KML 10/25/88
!
      if(acoef(296).eq.2._R8) then
      do 558 jg=1,5
 558  cisum(jg)=0._R8
      vtran=rt*usdr*(cegp1(1)+cegp1(2)+cegp1(3)+cegp1(4))
      endif
!
!@@@a
      if(kcycle.le.1) then
        do 556 jg=1,nreg
          if (irst1 .ne. 2 .or. acoef(296) .ne. 4) then
            ceg(jg)   = veg(jg) / rsceg(jg)
          else
            if (jg.eq.1) ceg(jg) = (ivist(jg)) * usdi
            if (jg.eq.2) ceg(jg) = (ivist(jg)) * usdi
            if (jg.eq.3) ceg(jg) = ivist(jg) * usdi
            if (jg.eq.4) ceg(jg) = ivist(jg) * usdi
            if (jg.eq.5) ceg(jg) = ivist(jg) * usdi
            if (jg.eq.6) ceg(jg) = ivist(jg) * usdi
            if (jg.eq.7) ceg(jg) = iohist    * usdi
          endif
          cegp1(jg) = ceg(jg)
          cegm1(jg) = ceg(jg)
          ceg0m1(jg) = ceg0(jg)
          ipowm1(jg)  = ipow(jg)
          ipowm2(jg)  = ipow(jg)
          flxegm1(jg) = flxeg(jg)
          ipowp (jg)  = 0.0_R8
          flxp  (jg)  = 0.0_R8
 556    continue
          ipowm1(8)  = ipow(8)
          ipowm2(8)  = ipow(8)
          ipowm1(9)  = ipow(9)
          ipowm2(9)  = ipow(9)
!@@@e
      energyr = 0._R8
!
        goto 590
      endif
!@@@a
!     initialize ipowm1, ipowm2 if restart run:
      if(iswtch.eq.0) then
        iswtch = 1
        do 10030 ig = 1,nreg
          ipowm1(ig)  = ipow(ig)
          ipowm2(ig)  = ipow(ig)
          ipowp (ig)  = 0.0_R8
10030   continue
          ipowm1(8)  = ipow(8)
          ipowm2(8)  = ipow(8)
          ipowm1(9)  = ipow(9)
          ipowm2(9)  = ipow(9)
      endif
!@@@e
!
!.............................................................................
!..... 4.0 calculate OH voltage and define wrk1 for diagnostic printout only
!.............................................................................
      if(iprnt.eq. nskipr+1 ) write(nout,1000) times,kcycle
      if(acoef(296).eq.4._R8) go to 11110
      do 1560 ig=1,nreg
      cxave(jnreg(ig)) = 0.5_R8*(cx(jnreg(ig)) + cxold2(jnreg(ig)) )
      cxold2(jnreg(ig)) = cxave(jnreg(ig))
 1560 continue
!
      do 560 jg=1,nreg
      sumcx = 0
      if(idata.ne.7) go to 561
      do 559 ig=1,nreg
      sumcx = sumcx + ameg(jg,ig)*cxave(jnreg(ig))
  559 continue
      sumcx = sumcx + rsceg(jg)*((dt+dtold)*cxave(jnreg(jg))+cegm1(jg))
  561 continue
!
        voheg(jg) = (isym+1)*wrk2(jg)*(-tpi*reboun)
!
        wrk1(jg) = vegm1(jg) - rsceg(jg)*cegm1(jg)                       &  
     &           - (flxeg(jg) - flxegm2(jg))/(dt+dtold)                  &  
     &                         + gvolt0(jnreg(jg))*aturneg(jg)*usdv      &  
     &           + sumcx
 
 560  continue
      if(acoef(296).eq. 2._R8) go to 905
!
!..............................................................................
!.....5.0 calculate new external group current (cegp1):
!.............................................................................
!
!..first define matrices needed to solve circuit equation:
!
!  ameg(ig,jg)*(cegp1(ig) - cegm1(ig))
!             + (dt+dtold)*(rsceg(ig)*cegp1(ig)) + (flxeg(ig)-flxegm2(ig))
!   = (dt+dtold)*(            gvolt0(jnreg(ig))*aturneg(ig)*usdv
!                 + ameg(ig,jg)*cxave(jnreg(jg)))
!     +              cp*(dt+dtold)*(cegp1-ceg0)
!     +              cd*((cegp1-cegm1)-(ceg0-ceg0m2))
!     + ameg(ig,jg)*(bp*(dt+dtold)*(cegp1-ceg0)
!                   +bd*((cegp1-cegm1)-(ceg0-ceg0m2)))
!
!
!   or....
!
!     a1(ig,jg)*cegp1(jg) - a2(ig,jg)*cegm1(jg) - a3(ig,jg)*ceg0(jg)
!                         + a4*ceg0m2(jg) = rr(ig)
!
!
!.....turn on system slowly to avoid voltage spikes
      f1 = 1._R8
      if(kcycle.lt.200) f1 = (kcycle/200._R8)**4
!
!
      if (acoef(296) .eq. 6._R8) then
         call volt(v(1),v(2),v(3),v(4),v(5),v(6),v(7))
       end if
      if (acoef(296) .eq. 7._R8) then
         call volt2(v(1),v(2),v(3),v(4),v(5),v(6),v(7))
      end if
!
      do 790 ig=1,nreg
      if(acoef(296).eq.1) go to 775
      if(acoef(296).eq.6) go to 772
      if(acoef(296).eq.7) go to 773
      bp = gainpeg(ig)*udsi/udsv*f1
      bd = gaindeg(ig)*udsi/udst/udsv*f1
      cp = 0._R8
      cd = 0._R8
      go to 779
  775 continue
      bp = 0._R8
      bd = 0._R8
      cp = gainpeg(ig)*udsi/udsv*f1
      cd = gaindeg(ig)*udsi/udst/udsv*f1
      go to 779
  773 bp = 0._R8
      bd = 0._R8
      cp = 0._R8
      cd = 0._R8
      gvolt0(jnreg(ig)) = v(jnreg(ig))
      go to 779
  772 bp = 0._R8
      bd = 0._R8
      cp = 0._R8
      cd = 0._R8
!CC... modify for 1-turn coils - remove aturneg(ig)
      if(acoef(296) .eq. 6._R8) then
       gvolt0(jnreg(ig)) = v(jnreg(ig))
         else
       gvolt0(jnreg(ig)) = v(jnreg(ig))/aturneg(ig)
       endif
!CC
  779 continue
      sumcx = 0
      do 780 jg=1,nreg
      a1(ig,jg) = (1._R8+(dt+dtold)*bp+bd)*ameg(ig,jg)
      a2(ig,jg) = (1._R8+bd)*ameg(ig,jg)
      a3(ig,jg) = ((dt+dtold)*bp + bd)*ameg(ig,jg)
      a4(ig,jg) = bd*ameg(ig,jg)
      if(idata.ne.7) go to 780
      sumcx = sumcx + ameg(ig,jg)*cxave(jnreg(jg))
  780 continue
      if(idata.ne.7) go to 783
      sumcx = sumcx + rsceg(ig)*((dt+dtold)*cxave(jnreg(ig))+cegm1(ig))
  783 continue
      a1(ig,ig) = a1(ig,ig)+(dt+dtold)*(rsceg(ig)+cp)+cd
      a2(ig,ig) = a2(ig,ig)+cd
      a3(ig,ig) = a3(ig,ig)+(dt+dtold)*cp+cd
      a4(ig,ig) = a4(ig,ig)+cd
      rr(ig) = (dt+dtold)*(gvolt0(jnreg(ig))*                            &  
     &    usdv*aturneg(ig) + sumcx) - (flxeg(ig)-flxegm2(ig))
  790 continue
!..invert a1
      call minv(a1,nreg,pngroup,wrk3,det,epsminv,0,1)
!..now calculate current at advanced time
      do 770 ig=1,nreg
      sum = 0
      do 760 jg=1,nreg
      do 750 kg=1,nreg
      sum = sum + a1(ig,jg)*( a2(jg,kg)*cegm1(kg)                        &  
     &                      + a3(jg,kg)*ceg0(kg)                         &  
     &                      - a4(jg,kg)*ceg0m2(kg) )
  750 continue
      sum = sum + a1(ig,jg)*rr(jg)
  760 continue
      cegp1(ig) = sum
  770 continue
!
!     write(nterm,3771)
!     do 1771 ig=1,nreg
!     write(nterm,2771) ig,gainpeg(ig),gaindeg(ig),rr(ig),
!    1    cegp1(ig),cegm1(ig),ceg0(ig)
!2771 format(i3,1p6e10.2)
!1771 continue
!3771 format(" ig   gainpeg   gaindeg        rr     cegp1     cegm1",
!    1       "      ceg0")
!
!
!..............................................................................
!..... 6.0 calculate voltage applied to each coil group
!..............................................................................
!
!.....calculate voltage applied to each coil group
!     write(nterm,2776)
 2776 format("  ig      veg0     cegm1     rsceg     ceg0")
      do 776 ig=1,nreg
      veg0(ig) = 0._R8
      if(times.gt.tauseg(ig)) go to 776
      if(times.lt.teineg(ig)) go to 776
      sum = 0._R8
      do 785 jg=1,nreg
      sum = sum + a3(ig,jg)*(ceg0(jg)-cegp1(jg))                         &  
     &          - a4(ig,jg)*(ceg0m2(jg)-cegm1(jg))                       &  
     &          + (dt+dtold)*ameg(ig,jg)*cxave(jnreg(jg))
  785 continue
      if(idata.ne.7) go to 786
      sum = sum +                                                        &  
     & (dt+dtold)*(rsceg(ig)*( (dt+dtold)*cxave(jnreg(ig))+cegm1(ig)))
  786 continue
      vdt(ig) = sum + (dt+dtold)*( gvolt0(jnreg(ig))*usdv*aturneg(ig))
      if(dt+dtold .le. 0) go to 776
      veg0(ig) = vdt(ig)/(dt+dtold)
!     write(nterm,1776) ig,veg0(ig),cegm1(ig),rsceg(ig),ceg0(ig)
 1776 format(i3,1p4e10.2)
  776 continue
!     do 2778 ig=1,nreg
!     do 2778 jg=1,nreg
!     write(nterm,2777) ig,jg,ameg(ig,jg)
!2778 continue
!2777 format(1x,2i3,1pe12.4)
!     stop
!
      go to 988
!..............................................................................
!..... 7.0 For ZT-H only   [ acoef(296)=2 ]
!..............................................................................
!
!..KML 10/25/88. Acoef(296)=2. specifically uses the ZT-H circuit.
!..... a1,a2,a3,a4 are 5x5 matrices
!
 905  continue
      esc=acoef(869)*1000._R8/udsv
      if(vtran.gt.esc) rs=1.E+120_R8
      if(vtran.le.esc) rs=.0075_R8
      rsc=rs/udsr
      rtc=rt/udsr
      rf1=rtc/(rtc + rsc)
      rf2=1._R8/(1._R8/rsc + 1._R8/rtc)
      do 802 ig=1,5
      do 802 jg=1,5
      a2(ig,jg)=ameg(ig,jg)
      a3(ig,jg)=0._R8
      a4(ig,jg)=0._R8
 802  continue
      do 812 ig=1,5
      cpv(ig) = gainpeg(ig)*udsi/udsv*f1
      cdv(ig) = gaindeg(ig)*udsi/udst/udsv*f1
      civ(ig) = gainieg(ig)*udsi*udst/udsv*f1
      rr(ig) =  - (flxeg(ig)-flxegm2(ig))
      a1(ig,5)=ameg(ig,5)
      a1(5,ig)=ameg(5,ig)
  812 continue
      a1(5,5)=a1(5,5) + rsceg(5)*(dt + dtold)
!
      do 817 ig=1,5
      cisum(ig)=cisum(ig)+0.25_R8*civ(ig)*(dt+dtold)*(ceg0m1(ig)         &  
     & +ceg0m2(ig) - ceg(ig) - cegm1(ig))
 817  continue
!
      do 822 ig=1,4
      do 821 jg=1,4
      a1(ig,jg)=ameg(ig,jg) + (dt + dtold)*rf2
 821  continue
      a1(ig,ig)=a1(ig,ig) + (dt + dtold)*rsceg(ig)
      a3(ig,1)=-rf1*(dt + dtold)*cpv(1) - rf1*cdv(1)
      a4(ig,1)=-rf1*cdv(1)
      rr(ig)=rr(ig) - (dt+dtold)*cisum(1)*rf1
 822  continue
!
      do 831 ig=2,5
      a1(ig,ig) = a1(ig,ig) + (dt+dtold)*cpv(ig) + cdv(ig)
      a2(ig,ig) = a2(ig,ig) + cdv(ig)
      a3(ig,ig)=(dt + dtold)*cpv(ig) + cdv(ig)
      a4(ig,ig)=cdv(ig)
      rr(ig)=rr(ig) + (dt+dtold)*cisum(ig)
 831  continue
      do 833 ig=1,4
      a1(ig,1) = a1(ig,1) - rf1*(dt+dtold)*cpv(1)  - rf1*cdv(1)
      a2(ig,1) = a2(ig,1) - rf1*cdv(1)
 833  continue
!
!..invert a1
      call minv(a1,5,pngroup,wrk3,det,epsminv,0,1)
!..now calculate current at advanced time
      do 880 ig=1,5
      sum = 0
      do 860 jg=1,5
      do 850 kg=1,5
      sum = sum + a1(ig,jg)*( a2(jg,kg)*cegm1(kg)                        &  
     &                      + a3(jg,kg)*ceg0(kg)                         &  
     &                      - a4(jg,kg)*ceg0m2(kg) )
  850 continue
      sum = sum + a1(ig,jg)*rr(jg)
  860 continue
      cegp1(ig) = sum
  880 continue
!
!
!.....calculate voltage applied to each coil group
      do 876 ig=1,5
      vdt(ig)=cdv(ig)*(ceg0(ig)-ceg0m2(ig))-cdv(ig)*(cegp1(ig)-cegm1(ig)  &  
!    &                                                                   &  
     & )+(dt+dtold)*cpv(ig)*(ceg0(ig)-cegp1(ig))+ (dt+dtold)*cisum(ig)
      if(dt+dtold .le. 0) go to 876
      veg0(ig) = vdt(ig)/(dt+dtold)
  876 continue
      vtran=veg0(1)*rf1 + rf2*(cegp1(1)+cegp1(2)+cegp1(3)+cegp1(4))
!      if(mod(kcycle,10).eq.0)write(nterm,3)esc,vtran,(veg0(ig),ig=1,5)
!     x,(veg(ig),ig=1,5)
!      if(mod(kcycle,10).eq.0)write(nout,3)esc,vtran,(veg0(ig),ig=1,5)
!     x,(veg(ig),ig=1,5)
! 3    format(12e8.1)
!
!
!
 988  continue
!
!..............................................................................
!.....8.0  limit group voltages by power supplies
!..............................................................................
      do 789 ig=1,nreg
      vmin = vmineg(ig)*usdv*aturneg(ig)
      vmax = vmaxeg(ig)*usdv*aturneg(ig)
      veg(ig) = max(veg0(ig),vmin)
      veg(ig) = min(veg(ig),vmax)
      vdiff(ig) = 0._R8
      if(veg0(ig).eq.0) go to 788
      ratio = veg(ig)/veg0(ig)
      if(ratio .lt. 0) ratio = 0._R8
      if(ratio .gt. 1._R8) ratio= 1._R8
      vdiff(ig) = (1._R8-ratio)*vdt(ig)
  788 continue
  789 continue
!
      do 881 ig=1,nreg
      do 882 jg=1,nreg
  882 a1(ig,jg) = ameg(ig,jg)
  881 a1(ig,ig) = a1(ig,ig)+(dt+dtold)*rsceg(ig)
      call minv(a1,nreg,pngroup,wrk3,det,epsminv,0,1)
      do 781 ig=1,nreg
      sum = 0._R8
      do 782 jg=1,nreg
  782 sum = sum + a1(ig,jg)*vdiff(jg)
  781 cegp1(ig) = cegp1(ig) - sum
      go to 11100
 
!--------------------------------------------------------
!  explicit method to calculate new external coil currents
!    for acoef(296) = 4.
! -------------------------------------------------------
11110 continue
!      --------------------------------
!      woyke-feedback fuer acoef(296) = 4.
!      --------------------------------
!         uebergabe der
!         stromableitungen (waveforms)
!         der induktion von inneren kreisen auf auessere,
!         integration der abweichung
!         und tracking
!         in regler wenn acoef(296) = 4.
!         ------------------------------------
           npro = 9
            do 20000 jg = 1, npro
              ipowp (jg) = (ipow(jg) - ipowm2(jg))                       &  
     &        /((dt+dtold)*udst)
              flxp  (jg) = (flxeg(jg) - flxegm2(jg))/(dt + dtold)*udsv
              iffint(jg) = iffint (jg) + (ipow (jg) - ivff (jg))         &  
     &                     * dt * udst
20000      continue
         do 11115 ig=1,nreg
!
!        regeldifferenz und uebergabe der ist werte fuer stroeme
          if((ig.ge.1) .and. (ig.le.6)) then
            ivist (ig)  = ceg (ig) * udsi
            diff        = ivsoll (ig) - ivist (ig)
          endif
          if(ig.eq.7) then
            iohist  = ceg (ig) * udsi
            diff    = iohsoll - iohist
          endif
 
!        integration
          vegint(ig) = vegint (ig) + diff * dt * udst
!        regelgesetz
          veg   (ig) = gainpeg (ig) * diff + gainieg (ig) * vegint (ig)
!        beschraenkung der stellgroesse
          if ((veg(ig) .lt. vmineg (ig)) .or. (veg(ig) .gt. vmaxeg(ig)))  &  
!    &                                                                   &  
     &              vegint (ig) = vegint (ig) - diff * dt * udst
          veg (ig) = min (veg (ig), vmaxeg (ig))
          veg (ig) = max (veg (ig), vmineg (ig))
!        normierung der stellgroesse
          veg (ig) = veg (ig) * usdv
 
11115 continue
!      ------------------------------------
!       ende feedback (w. woyke)
!      ------------------------------------
!..... sum up voltages (rhs of equation):
!     voltage from ideal heating system:
!     (-reboun*tpi: voltage contribution from ideal ohmic heating
!      system per coil)
 
      do 11120 jg=1,nreg
        voheg(jg) = (isym+1)*wrk2(jg)*(-tpi*reboun)
        wrk1(jg) = vegm1(jg) - rsceg(jg)*cegm1(jg)                       &  
     &           - (flxeg(jg) - flxegm2(jg))/(dt+dtold)                  &  
     &           + vohegm1(jg)
11120 continue
 
!.....calculate new coil group currents:
      do 11140 ig = 1,nreg
        sum = cegm1(ig)
        do 11130 jg = 1,nreg
          sum = sum + ameginv(ig,jg) * wrk1(jg) *(dt+dtold)
11130   continue
      cegp1(ig) = sum
11140 continue
!--------------------------------------------------------
11100 continue
!.....output:
!
      if(iprnt.le.nskipr) goto 511
!
      do 501 jg=1,nreg
      iut1   = jg
      iut2   = jnreg(jg)
      out1   = ceg0(jg)                     *udsi*1.E-3_R8
      out2   = cegp1(jg)                     *udsi*1.E-3_R8
      out3   = out1 - out2
      out4   = veg(jg)                      *udsv
      out5 = 0._R8
      out6 = 0._R8
      out7 = (flxeg(jg)-flxegm2(jg))/(dt+dtold)*udsv
      do 499 ig=1,nreg
      if(acoef(296).eq.1._R8) go to 975
      bp = gainpeg(ig)*udsi/udsv*f1
      bd = gaindeg(ig)*udsi/udst/udsv*f1
      cp = 0._R8
      cd = 0._R8
      go to 979
  975 continue
      bp = 0._R8
      bd = 0._R8
      cp = gainpeg(ig)*udsi/udsv*f1
      cd = gaindeg(ig)*udsi/udst/udsv*f1
  979 continue
      out5 = out5 + bp * ameg(jg,ig)*(ceg0(ig)-cegp1(ig))*udsv
      out6 = out6 + bd * ameg(jg,ig)*(ceg0(ig)-cegp1(ig)                 &  
     &                               -ceg0m2(ig)+cegm1(ig))*udsv
  499 continue
      out5 =out5+cp*(ceg0(jg)-cegp1(jg))                *udsv
      out6 =out6+cd*((ceg0(jg)-cegp1(jg))-(ceg0m2(jg)-cegm1(jg)))*udsv
!
      write(nout,1001) iut1,iut2,                                        &  
     &                 out1,out2,out3,out4,out5,out6,out7
!
  501 continue
 511  continue
!
!...now calculate inductive and resistive power for output
      poweri = 0._R8
      powerr = 0._R8
      do 580 ig=1,nreg
!
      poweri = poweri + ceg(ig)*wrk1(ig)*udsv*udsi*1.E-6_R8
      powerr = powerr + rsceg(ig)*cegm1(ig)*cegm1(ig)*udsr*udsi**2*      &
     &         1.E-6_R8
 580  continue
!
 590  continue
!
      energyi = 0
      do 591 ig=1,nreg
      do 591 jg=1,nreg
  591 energyi = energyi + .5_R8*ameg(ig,jg)*ceg(ig)*ceg(jg)
      energyi = energyi*udsp*1.E-6_R8
      energyr = energyr + powerr*dts
!
!.....new coil currents:
      do 610 jg=1,nreg
        do 600 ii=1,ncoil-nwire
          if(igroupc(ii) .ne. jnreg(jg)) goto 600
          ccoil (ii) = aturnsc(ii)*cegp1(jg)
          ccoils(ii) = ccoil(ii)*udsi
!.........voltages of external coils to be plottet:
!@@@a
!===      vegplot(ii) = aturnsc(ii)*veg(jg)/wrk2(jg)
          vegplot(ii) = aturnsc(ii)*veg(jg)*(1+isym)                     &  
     &                / aturneg(jg)
!@@@e
 600    continue
 610  continue
!
 
!...save time step, fluxes and currents:
!   time step:
      dtm2 = dtm1
      dtm1 = dt
      do 700 jg=1,nreg
!   fluxes:
        flxegm2(jg) = flxegm1(jg)
        flxegm1(jg) = flxeg  (jg)
!   group currents:
        cegm1(jg) = ceg  (jg)
        ceg  (jg) = cegp1(jg)
        ceg0m2(jg) = ceg0m1(jg)
        ceg0m1(jg) = ceg0 (jg)
!       -------------------------
! ***   a. krause 20.04.89  *****
        ipowm2(jg) = ipowm1(jg)
        ipowm1(jg) = ipow  (jg)
! *******************************
!   group voltages:
        vegm1(jg) = veg(jg)
!   oh-group voltage:
        vohegm1(jg) = voheg(jg)
!
 700  continue
          ipowm2(8)  = ipowm1(8)
          ipowm1(8)  = ipow(8)
          ipowm2(9)  = ipowm1(9)
          ipowm1(9)  = ipow(9)
!
 900  continue
      return
!     -------------
      entry cirequi
!     -------------
!
!.....initialization for circuit equations
!     --------------
!
      if(icirc.eq.0)       goto 300
      if(iswitch.eq.1)     goto 300
      iswitch = 1
      if(ncoil.eq.nwire)   goto 300
!
!.....order group numbers:
      do 139 ii = 1,ncoil-nwire
        jtemp(ii) = igroupc(ii)
 139  continue
      do 141 i1 = 1,ncoil-nwire
        do 140 i2 = 2,i1
          if( jtemp(i1) .gt. jtemp(i2-1) ) goto 140
          itmp        = jtemp(i2-1)
          jtemp(i2-1) = jtemp(i1)
          jtemp(i1)   = itmp
 140    continue
 141  continue
!
!.....determine total number of external groups (nreg):
      nreg = 1
      jnreg(1) = jtemp(1)
      do 142 ii = 1,ncoil-nwire -1
        if (jtemp(ii+1) .gt. jtemp(ii) ) then
          nreg = nreg + 1
          jnreg(nreg) = jtemp(ii+1)
!
        endif
 142  continue
!.....jnreg(nreg) contains now the group numbers of ext. coils
!     in increasing order
!
!.....begin calculation of self and mutual inductances and the
!     resistance of external coil groups:
!
      do 200 jg=1,nreg
!
!.......sum up the self inductances and resistances:
        iic=0
        sumdia = 0._R8
        rsceg(jg)=0._R8
        do 151 ii=1,ncoil-nwire
          if( igroupc(ii) .ne. jnreg(jg) ) goto 151
          iic = iic + 1
          xtmp1(iic) = xcoil(ii)
          ztmp1(iic) = zcoil(ii)
          atrn1(iic) = aturnsc(ii)
!..rxw/14/01/88
          dxtmp1(iic) = dxcoil(ii)
          dztmp1(iic) = dzcoil(ii)
!...rxw
!
!         sum up diagonal elements of coil inductances in a group:
          sumdia = sumdia + abs( aindc(ii)*aturnsc(ii)**2 )
!
!         resistance:
          rsceg(jg) = rsceg(jg) + abs( rscoil(ii)*aturnsc(ii)**2 )
!         -------
!
 151    continue
        njg = iic
        wrk2(jg) = sumdia
        rsceg(jg) = (isym+1) * rsceg(jg)
        wrk2(jg)  = (isym+1) * wrk2(jg)
!
!       sum up off diagonal elements of coil inductances in a group:
        sumoff = 0._R8
        do 155 i=1,iic-1
          do 153 j=i+1,iic
!..rxw/14/01/88
           if( (dxtmp1(i).le.0) .or. (dztmp1(i).le.0) .or.               &  
     &         (dxtmp1(j).le.0) .or. (dztmp1(j).le.0) ) goto 851
           ans = -amutlrc(xtmp1(i),ztmp1(i),dztmp1(i),dxtmp1(i),1._R8,   &  
     &                  xtmp1(j),ztmp1(j),dztmp1(j),dxtmp1(j),1._R8,20)  &  
     &                  /usdi
           goto 852
 851       continue
       call gf(ineg,nmult,xtmp1(i),ztmp1(i),xtmp1(j),ztmp1(j),ans)
 852       continue
!...rxw
           sumoff = sumoff -  atrn1(i)*ans*atrn1(j)
 153      continue
 155    continue
!
        alin = sumdia + 2._R8*sumoff
!
!.....additional diagonal components for symmetric case:
        sumsym = 0._R8
        if(isym.eq.1) then
          do 158 i=1,njg
           do 157 j=1,njg
!..rxw/14/01/88
           if( (dxtmp1(i).le.0) .or. (dztmp1(i).le.0) .or.               &  
     &         (dxtmp1(j).le.0) .or. (dztmp1(j).le.0) ) goto 853
           ans = -amutlrc(xtmp1(i),ztmp1(i),dztmp1(i),dxtmp1(i),1._R8,   &  
     &                   xtmp1(j),-ztmp1(j),dztmp1(j),dxtmp1(j),1._R8,   &  
     & 20)                                                               &  
     &                   /usdi
           goto 854
 853       continue
       call gf(ineg,nmult,xtmp1(i),ztmp1(i),xtmp1(j),-ztmp1(j),ans)
 854       continue
!...rxw
              sumsym = sumsym -  atrn1(i)*ans*atrn1(j)
 157       continue
 158      continue
        endif
!
!.....self inductance of external coil group:
        if(isym.eq.0) ameg(jg,jg) = alin
        if(isym.eq.1) ameg(jg,jg) = 2._R8*(alin + sumsym)
!
!
!.....calculate mutual inductance of external coil groups:
!
        do 180 kg=jg+1,nreg
          iic = 0
          do 162 ii=1,ncoil-nwire
            if( igroupc(ii) .ne. jnreg(kg) ) goto 162
            iic = iic+1
            xtmp2(iic) = xcoil(ii)
            ztmp2(iic) = zcoil(ii)
            atrn2(iic) = aturnsc(ii)
!..rxw/14/01/88
            dxtmp2(iic) = dxcoil(ii)
            dztmp2(iic) = dzcoil(ii)
!...rxw
 162      continue
          nkg = iic
!
          if(isym.eq.1) goto 170
!
!.........calculate mutal inductance matrix for unsymmetric case:
          sum = 0._R8
          do 166 i=1,njg
            do 164 j=1,nkg
!..rxw/14/01/88
           if( (dxtmp1(i).le.0) .or. (dztmp1(i).le.0) .or.               &  
     &         (dxtmp2(j).le.0) .or. (dztmp2(j).le.0) ) goto 855
           ans = -amutlrc(xtmp1(i),ztmp1(i),dztmp1(i),dxtmp1(i),1._R8,   &  
     &                   xtmp2(j),ztmp2(j),dztmp2(j),dxtmp2(j),1._R8,20)  &  
!    &                                                                   &  
     &                   /usdi
             goto 856
 855        continue
      call gf(ineg,nmult,xtmp1(i),ztmp1(i),xtmp2(j),ztmp2(j),ans)
 856        continue
!...rxw
              sum = sum - atrn1(i)*ans*atrn2(j)
 164        continue
 166      continue
          goto 179
!
!.........mutual inductance matrix for symmetric case:
 170      continue
          sum = 0._R8
          do 176 i=1,njg
            do 199 j=1,nkg
!..rxw/14/01/88
           if( (dxtmp1(i).le.0) .or. (dztmp1(i).le.0) .or.               &  
     &         (dxtmp2(j).le.0) .or. (dztmp2(j).le.0) ) goto 857
!
           ans1 = -amutlrc(xtmp1(i),ztmp1(i),dztmp1(i),dxtmp1(i),1._R8,  &  
     &                   xtmp2(j),ztmp2(j),dztmp2(j),dxtmp2(j),1._R8,20)  &  
!    &                                                                   &  
     &                   /usdi
           ans2 = -amutlrc(xtmp1(i),ztmp1(i),dztmp1(i),dxtmp1(i),1._R8,  &  
     &                   xtmp2(j),-ztmp2(j),dztmp2(j),dxtmp2(j),1._R8,   &  
     & 20)                                                               &  
     &                   /usdi
           goto 858
 857       continue
      call gf(ineg,nmult,xtmp1(i),ztmp1(i),xtmp2(j), ztmp2(j),ans1)
      call gf(ineg,nmult,xtmp1(i),ztmp1(i),xtmp2(j),-ztmp2(j),ans2)
 858       continue
!...rxw
              sum = sum - atrn1(i)*(ans1+ans2)*atrn2(j)
 199        continue
 176      continue
          sum = 2._R8*sum
!
!
!.....mutual inductances of external coil groups:
!
 179      continue
          ameg(jg,kg) = sum
          ameg(kg,jg) = ameg(jg,kg)
!         -----------
!
 180    continue
!
 200  continue
!
!.....invert inductance matrix of external coil groups
!
      do 230 i=1,nreg
        do 230 j=1,nreg
          ameginv(i,j) = ameg(i,j)
 230  continue
!
!     use fast cray routine minv to invert matrix
      epsminv=1.E-15_R8
      call minv(ameginv,nreg,pngroup,wrk1,det,epsminv,0,1)
!
!.....inverse of inductance matrix of external coil group is now
!     stored in ameginv
!
      do 290 jg=1,nreg
        wrk2(jg) = 0._R8
!
        do 245 iw=1,nwire
          atans(jg,iw) = 0._R8
 245    continue
!
        do 280 ic=1,ncoil-nwire
          if(igroupc(ic) .ne. jnreg(jg)) goto 280
!         sum up winding numbers of each external group:
          wrk2(jg) = wrk2(jg) + abs(aturnsc(ic))
!
!.........calculate greensfunction between ext. coils and internal wires:
!
          do 250 iw=1,nwire
        call gf(ineg,nmult,xcoil(ic),zcoil(ic),xwire(iw),zwire(iw),ans)
            atans(jg,iw) = atans(jg,iw) - aturnsc(ic)*ans
 250      continue
          if(isym.eq.0) goto 280
!
!.........additional contribution if isym=1:
          do 260 iw=1,nwire
        call gf(ineg,nmult,xcoil(ic),zcoil(ic),xwire(iw),-zwire(iw),ans)    
            if(zwire(iw).eq.0._R8) ans = 0._R8
            atans(jg,iw) = atans(jg,iw) - aturnsc(ic)*ans
 260      continue
!
 280    continue
 
      aturneg(jg) = wrk2(jg)*(1._R8+isym)
 290  continue
!......check if inductance matrix is diagonally dominant
      do 1314 i=1,nreg
      do 1315 j=1,nreg
      jj = j
      if(j.eq.i) go to 1315
      if(abs(ameg(i,j)*ameg(j,i)) .gt. abs(ameg(i,i)*ameg(j,j)))         &  
     &  go to 1316
 1315 continue
 1314 continue
      go to 1317
 1316 continue
      ineg=30
      write(nout,3455) jj
 3455 format(" * * * error, row ",i3," in inductance matrix",            &  
     &" has off-diagonal element larger than diagonal")
 3456 format(" i,j,ameg(i,j)*usdi",2i3,1pe12.4)
!
!    print out groupc, self and mutual inductances:
       do 314 i=1,nreg
        do 315 j=1,nreg
         write(nout,3456) i,j,ameg(i,j)*usdi
  315   continue
  314  continue
!
!@@@a rxw/30/3/89
!== 301  continue
!@@@e
 1317 continue
!     -------------------
!.....initialize feedback on external group currents:
!     -------------------
!     feedback on external coil group currents
!
!     acoef-array:
!       acoef(294) >  0 --> call feedeg
!       acoef(294): nstart: first acoef-coefficient
!       acoef(295): nfb: total number of coil groups for feedback
!       acoef(296) = 1: standard pid-feedback model
!                  = else:  feedback coefficients are
!                           multiplied by inductance matrix
!       acoef(j):
!        nstart       < j < nstart+  nfb-1 :  first acoef-coefficient
!        gains:
!        nstart+nfb   < j < nstart+2*nfb-1 :  gainpeg : proportional
!        nstart+2*nfb < j < nstart+3*nfb-1 :  gaindeg : differential
!        nstart+3*nfb < j < nstart+4*nfb-1 :  gainieg : integral
!
!        nstart+4*nfb < j < nstart+5*nfb-1 :  teineg
!        nstart+5*nfb < j < nstart+6*nfb-1 :  tauseg
!
!       teineg:  time when feedback is switched on
!       tauseg:  time when feedback is switched off
!
!     -------------
      nstart = acoef(294)
      nfb    = acoef(295)
      nn=nstart-1
      nd=nfb
!
!.....assign current feedback parameters:
!
      write(nout,1002)
!
      do 190 jg = 1,nreg
!
        lgroup(jg) = 0
        teineg(jg) = 0._R8
        tauseg(jg) = 10000._R8
!
        do 185 i = 1,nd
          if( int(acoef(nn+i)) .eq. jnreg(jg) ) then
!
!...........convert feedback coefficients
!           proportional:
            gainpeg(jg) =  acoef(nn+  nd+i)
!           differential:
            gaindeg(jg) =  acoef(nn+2*nd+i)
!           integral:
            gainieg(jg) =  acoef(nn+3*nd+i)
!
!...........time when feedback system is switched on and off:
            teineg(jg) =  acoef(nn+4*nd+i)
            tauseg(jg) =  acoef(nn+5*nd+i)
!
!.....min and max voltages
      vmineg(jg) = acoef(nn+6*nd+i)*1000._R8
      vmaxeg(jg) = acoef(nn+7*nd+i)*1000._R8
      if(vmineg(jg) .eq. 0) vmineg(jg) = -1.E6_R8
      if(vmaxeg(jg) .eq. 0) vmaxeg(jg) =  1.E6_R8
!
            if(tauseg(jg).eq.0._R8) tauseg(jg) = 10000._R8
!
            lgroup(jg) = 1
!
!...........print out information about feedback on external coil currents:
!
            out1 =  acoef(nn+  nd+i)
            out2 =  acoef(nn+2*nd+i)
            out3 =  acoef(nn+3*nd+i)
!
            write(nout,1003) jg,jnreg(jg),teineg(jg),tauseg(jg),         &  
     &                        out1,out2,out3,vmineg(jg),vmaxeg(jg)
!
          endif
 185    continue
 190  continue
!
!.....write output:
      write(nout,1006)
      do 240 jg=1,nreg
        out1 = rsceg(jg)   * udsr
        out2 = ameg(jg,jg) * usdi
        if(out1.ne.0) out3 =out2/out1
        write(nout,1005) jnreg(jg),out1,out2,out3,wrk2(jg)*usdi,         &  
     &     aturneg(jg)
 240  continue
!
!.....end of initialization
!     ---------------------
!
 300  continue
!
      return
!
 1000 format(//,' time = ',1pe12.4,4x,'kcycle=',i7/                      &  
     &       ' ncoil  igroupc',4x,'ceg0(ka) ',4x,'ceg(ka)  ',            &  
     &       3x,'diff(ka)',5x,'volts',                                   &  
     &       8x,'vprop(v) ',4x,'vdiff(v)     flxeg    ')
 1001 format(2i7,2x,8(1pe13.4))
 1002 format(//' feedback on external coil currents applied:'/           &  
     &       '  ncoil',3x,'igroup ',2x,'tfbon',8x,'tfboff',7x,           &  
     &       'vgainp',7x,'vgaind',7x,'vgaini',7x,' vmin ',7x,' vmax '/)
 1003 format(2i7,2x,7(1pe13.4))
 1006 format(/' external coils:'/' ngroup',7x,'resistance',5x,           &  
     &      'inductance',5x,'lr time',8x,                                &  
     &      'inductance (nc)  turns' )
 1005 format(i7,4x,1p5e15.7)
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
