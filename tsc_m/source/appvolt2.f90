      subroutine appvolt2
!     ------------------
!     subroutine used to model the asdex upgrade feedback2 system:
!     pid-feedback control on currents:   cwire(ii)-ccoil(i)
!     feedback system  can be switched on and off
! ***** acoef(290) = 4. for woyke model ********************
!
!     acoef-array:
!       acoef(290): = 2,4 --> call appvolt2
!                           ( in sr. initeq(7.00, initialization)
!                             and in sr. advan3(2.40) )
!       acoef(291): nstart (see below)
!       acoef(292): nfb: maximum number of feedback systems
!       acoef(293): nwprnt, to print out coil currents and voltages
!       acoef(j):
!        nstart       < j < nstart+  nfb-1 :  igroupw
!        nstart+nfb   < j < nstart+2*nfb-1 :  vgainp
!        nstart+2*nfb < j < nstart+3*nfb-1 :  vgaind
!        nstart+3*nfb < j < nstart+4*nfb-1 :  vgaini
!        nstart+4*nfb < j < nstart+5*nfb-1 :  tfbein
!        nstart+5*nfb < j < nstart+6*nfb-1 :  tfbaus
!        nstart+6*nfb < j < nstart+7*nfb-1 :  voltmax
!        nstart+7*nfb < j < nstart+8*nfb-1 :  tramp
!
!       tfbein:  time when feedback system is switched on
!       tfbaus:  time when feedback system is switched off
!       voltmax: maximum voltage
!       tramp:   ramp time
!
      USE CLINAM
      USE COMWOY
      USE SCR1
      USE TCVCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nasup,iswtch,nstart,nfb,nd,ii,jj,iabs,iiii,i,n
      INTEGER ig,igg,iw,jw,iut1,iut2
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 coiturn,tprein,tpraus,out1,out2
      REAL*8 out3,out4,out5,pt,skal,din,dio,vgainp,vgaind,vgaini
      REAL*8 vgainpm,dintg,vold,voltsav,dvolt,dvdt,dvdtmx,factor
      REAL*8 vgmax,out6,out7,out8
!============
!     dimension tfbein(pncoil),tfbaus(pncoil)
!     dimension cursum(pngroup), nsum(pngroup)
!============      
      parameter(nasup=96)
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tfbein
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tfbaus
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cursum
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nsum
!============      
      IF(.not.ALLOCATED(tfbein)) ALLOCATE( tfbein(pncoil), STAT=istat)
      IF(.not.ALLOCATED(tfbaus)) ALLOCATE( tfbaus(pncoil), STAT=istat)
      IF(.not.ALLOCATED(cursum)) ALLOCATE( cursum(pngroup), STAT=istat)
      IF(.not.ALLOCATED(nsum)) ALLOCATE( nsum(pngroup), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : appvolt2  ' 
!============      
!
!
      coiturn = 5._R8
!
!===  if(ifrst(10) .eq. 0) go to 300
      if(iswtch .eq. 1) goto 300
!
!.....initialize:
!===  ifrst(10) = 0
      iswtch = 1
      nstart = acoef(291)
      nfb    = acoef(292)
      nwprnt = acoef(293)
      nn=nstart-1
      nd=nfb
!
      tprein=1000._R8
      tpraus=0._R8
!
!.....assign  power supply parameters
!
      do 200 ii=1,nwire
        save1(ii)=0._R8
!
        do 100 jj=1,nd
          if( int(acoef(nn+jj)) .eq. iabs(igroupw(ii)) ) then
            save1(ii)=10._R8
            ilwire(ii)=0._R8
            if(acoef(290).eq.4) go to 100
!
!...........initial conditions on currents:
          if(kcycle.le.0) then
            diold(ii)=0._R8
            diol2(ii)=0._R8
           endif
!
!...........convert feedback coefficients to dimensionless units
!           proportional:
            vgain(ii)   =  acoef(nn+  nd+jj)
!           differential:
            vgain2(ii)  =  acoef(nn+2*nd+jj)
!           integral:
            vgain3(ii)  =  acoef(nn+3*nd+jj)
!
!...........time when feedback system is switched on and off:
!           time when feedback system is switched on:
            tfbein(ii)  =  acoef(nn+4*nd+jj)
!           time when feedback system is switched off:
            tfbaus(ii)  =  acoef(nn+5*nd+jj)
            save2(ii)   =  tfbein(ii)
!
            if(tfbaus(ii).eq.0._R8) tfbaus(ii)=100._R8
            tprein=min(tprein,tfbein(ii))
            tpraus=max(tpraus,tfbaus(ii))
!
!...........maximum voltage: (predefined: 1000 volts)
            voltmax(ii) = acoef(nn+6*nd+jj)
            if(voltmax(ii).eq.0._R8) voltmax(ii) = 1000._R8
!
!...........ramp time: (predefined: 1 ms)
            tramp(ii)   = acoef(nn+7*nd+jj)
            if(tramp(ii).eq.0._R8)   tramp(ii) = 1.E-3_R8
!
!...........print out information about feedback on currents:
!
            out1 =  acoef(nn+  nd+jj)
            out2 =  acoef(nn+2*nd+jj)
            out3 =  acoef(nn+3*nd+jj)
            out4 =  voltmax(ii)
            out5 =  tramp(ii)
!
            iiii = iiii+1
            if(iiii.eq.1) write(nasup,1002)
!
            write(nasup,1003) ii,igroupw(ii),tfbein(ii),tfbaus(ii),      &  
     &                        out1,out2,out3,out4,out5
!
          endif
 100    continue
 200  continue
!
!.....end of initialization
!
!.....calculate new feedback voltages:
!
 300  continue
      if(kcycle .le. 0) return
!
 
      if(acoef(290).ne.4) go to 441
!     4) feedback 2 kontrollspulen:
!     ================================
 
!     regeldifferenz:
 
      do 410 i = 1, 2
         fb2cd (i) = fb2ck1 (i)                                          &  
     &             * ( icsoll (i) - icist (i))
410   continue
 
!     pidt1 - kern:
 
      do 430 i =1, 2
        fb2cx  (i) = fb2ctn (i)
        fb2cei (i) = fb2cei (i) + fb2ctn (i) * fb2cd (i) * dt * udst
 
        ucsoll (i)   = fb2cei (i) + fb2cd (i)
 
430   continue
 
!     anti - windup massnahmen:
 
      do 440 i = 1, 2
 
       if (ucsoll (i) .lt. ucmin (i)) then
         ucsoll (i) = ucmin (i)
         if (fb2cx (i) .lt. 0) then
              fb2cei (i) = fb2cx (i)
         endif
       endif
 
       if (ucsoll (i) .gt. ucmax (i)) then
         ucsoll (i) = ucmax (i)
         if (fb2cx (i) .gt. 0) then
              fb2cei (i) = fb2cx (i)
         endif
       endif
 
440   continue
  441 continue
!
!.....calculate current in each group
      do 80 n=1,ngroup
      ig = nogroup(n)
      cursum(ig) = 0._R8
      nsum(ig) = 0._R8
      do 90 ii=1,nwire
      igg = igroupw(ii)
      if(igg.ne.ig) go to 90
      i = ncoil - nwire + ii
      cursum(ig) = cursum(ig) + ccoil(i)
      nsum(ig) = nsum(ig) + 1
   90 continue
      if(nsum(ig).ne.0) cursum(ig) = cursum(ig)/nsum(ig)
   80 continue
!
!.....loop:
      do 501 ii =1,nwire
!
        i   = ncoil - nwire + ii
        iw  = iwire(ii)
        jw  = jwire(ii)
        pt  = 0._R8
      etay(iw,jw) = rswire(ii)*dxdz/(tpi*xary(iw))
!
        ptsave(iw,jw) = -rswire(ii)*cwire0(ii)/tpi
        rewire(ii)    =  ptsave(iw,jw)
!@@@a
!      uebergabe der spannungen an kontrollspulen
!      w.woyke 20.7.89
!      a.krause 10.4.90
        if (acoef(290).eq.4) then
        if (igroupw(ii).eq.8 .or. igroupw(ii).eq.9) then
           skal           =   22.0E-06_R8/ (rswire (ii) * udsv / udsi)
           ptsave (iw,jw) = - ucsoll (ii) / skal / tpi * usdv / coiturn
           rewire (ii)    =   ptsave (iw, jw)
           resave (ii)    = - ptsave (iw, jw)
           save1 (ii)     = 0.0_R8
         endif
         go to 500
        endif
!@@@e
!
        din    = cwire0(ii) - ccoil(i)
        dio    = diold(ii)
        diold(ii) = din
!
        if(save1(ii) .lt. 1._R8) go to 500
        if( (times.lt.tfbein(ii)) .or. (times.gt.tfbaus(ii)) )           &  
     &    goto 500
!
!.....gain factors:
      vgainp = vgain(ii)*udsi/udsv
      vgaind = vgain2(ii)*udsi/udst/udsv
      if(vgain3(ii).gt.0)                                                &  
     &vgaini = vgain3(ii)*udst*udsi/udsv
!
      vgainpm = etav*tpi*xary(iw)/dxdz
      if(vgainp.gt.vgainpm) vgainp = vgainpm
!
      if(kcycle.eq.1) dio = din
      diol2(ii) = diol2(ii) + din*dt
      dintg     = diol2(ii)
!
!.....asdex-upgrade feedback system2:
!.....new voltages:
!
      vold = voltn(ii)
!
      ig = igroupw(ii)
      voltn(ii) = vgainp*din                                             &  
     &  +aindw(ii)*cx(ig)*aturnsw(ii)*udsi                               &  
     &          + vgaind*(din-dio)/dt                                    &  
     &          + vgaini*dintg
      voltsav = voltn(ii)
!
      dvolt = voltn(ii) - vold
      dvdt  = dvolt/dt
      dvdtmx= abs(voltmax(ii)*usdv/(tramp(ii)*usdt))
      if( abs(dvdt) .gt. dvdtmx ) then
        factor    = dvdtmx/abs(dvdt)
        voltn(ii) = vold + factor*dvdt*dt
      endif
      voltn(ii) = max(voltn(ii),-voltmax(ii)*usdv)
      voltn(ii) = min(voltn(ii), voltmax(ii)*usdv)
!
      factor = 1._R8
      if(voltsav.ne.0) factor = voltn(ii)/voltsav
!
      if(factor.lt.0._R8) factor = 0._R8
      if(factor.gt.1._R8) factor = 1._R8
!
      vgainp = vgainp*factor
      vgaind = vgaind*factor
      vgaini = vgaini*factor
!
      vgmax = etav*tpi*xary(iw)/dxdz
      if(vgainp.gt.vgmax) vgainp = vgmax
!
      etay(iw,jw) = (rswire(ii)+vgainp)*dxdz/(tpi*xary(iw))
      ptsave(iw,jw) = -vgainp*cwire0(ii)/tpi                             &  
     &  -aindw(ii)*cx(ig)*aturnsw(ii)*udsi/tpi                           &  
     &      - vgaind*(din-dio)/dt/tpi                                    &  
     &      - vgaini*dintg/tpi
      rewire(ii) = -vgainp*(cwire0(ii)-ccoil(i))/tpi                     &  
     &  -aindw(ii)*cx(ig)*aturnsw(ii)*udsi/tpi                           &  
     &      - vgaind*(din-dio)/dt/tpi                                    &  
     &      - vgaini*dintg/tpi
!
!
!.....output:
!.....write coil currents and voltages every "nwprnt" cycles
!
        if(kwprnt .ne. nwprnt) go to 499
!
        iut1   = ii
        iut2   = igroupw(ii)
        out2   = ccoil(i)                   *udsi*1.E-3_R8
        out3   = cwire0(ii)                 *udsi*1.E-3_R8
        out4   = out3 - out2
        out5   = -tpi*rewire(ii)            *udsv
        out6   = vgainp * din               *udsv
        out7   = vgaind * (din-dio) /dt     *udsv
        out8   = vgaini * dintg             *udsv
!
        write(nasup,1001) iut1,iut2,out2,out3,out4,out5,out6,out7,out8
  499 continue
      resave(ii) = voltn(ii)/tpi
!
!.....special logic to force series connections for vgain3 < 0
      if(vgain3(ii).lt.0) then
           ig = igroupw(ii)
           vgmax = etav*tpi*xary(iw)/dxdz/ndiv
!
           etay(iw,jw) = (rswire(ii))*dxdz/(tpi*xary(iw))
           ptsave(iw,jw) = -vgmax*(cursum(ig)-ccoil(i))/tpi
           rewire(ii) = -vgmax*(cursum(ig)-ccoil(i))/tpi
           go to 500
      endif
!
      go to 501
!
 500  continue
      resave(ii) = -ptsave(iw,jw)
 501  continue
!
      if(kwprnt.eq.nwprnt) kwprnt=0
!
      return
!
!..............................................................
 1000 format(//,' time = ',1pe12.4,4x,'kcycle=',i7/                      &  
     &       ' nwire  igroupw',4x,'ccoil(ka)',4x,'cwire0(ka)',           &  
     &       3x,'diff(ka)',5x,'volts',                                   &  
     &       8x,'vprop(v) ',4x,'vdiff(v)     vint(v)    ')
 1001 format(2i7,2x,8(1pe13.4))
 1002 format(//' feedback on currents applied:'/                         &  
     &       '   wire',3x,'igroupw',2x,'tfbon',8x,'tfboff',7x,           &  
     &       'vgainp',7x,'vgaind',7x,'vgaini',7x,                        &  
     &       'voltmax',6x,'tramp'/)
 1003 format(2i7,2x,7(1pe13.4))
!..............................................................
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
