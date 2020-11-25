      subroutine sawtooth
!
      USE sawtooth_mod
      USE trtsc
      USE CLINAM
      USE SAPROP
      USE SAWCOM

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER lprintdiag,ii,l,lsaw,j,is,issave,iisave
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 qcrit,f1saw,f2saw,denom,etasaw
      REAL*8 tflux,crashtimem,qave,sumo,psisaw,hyperfracc,gps,coef3
      REAL*8 rmid,emid,pmid,temidd,timidd,rminor1,rmajor1,elong1
      REAL*8 delta1,r1,s1,aave,deltaq,pi0,bp1sq,pleft,pright,prsaw
      REAL*8 vsaw,sumv,sumcp,sumal,rminorav,betap1,gval,gpval
      REAL*8 gppval,betai0,betaalpha,etparmks,valphen,rhoi,pileft
      REAL*8 piright,anleft,anright,pisaw,betai1,rp,rn,tau
      REAL*8 omegastari,omegastara,gammap,dwideal
      REAL*8 sum
!============
!     common /sawcom/  dwcore, dwcorec, dw, dwc, dwc2, ratio,scrit,
!    1dwbussac,dwelong,dwko,dwfast
      data lprintdiag/0/
!
!
!
!
!
!...... sawtooth models presently available:
!
!       isaw = 1    ...  traditional time-averaged TSC sawtooth model that
!                        uses qsaw, acoef(120), and acoef(124)
!       isaw = 2    ...  discrete sawtooth triggered by times on the
!                        type 75 input card
!       isaw = 3    ...  full Porcelli model
!
!
!.....define lsaw, f1saw, and f2saw
!
      if(isaw.le.1) then
        qcrit = qsaw
      else
        qcrit = 1.0_R8
      endif
      do 310 ii=1,npsit-1
      l = npsit-ii
      if(qprof2(l).lt.qcrit.and.qprof2(l+1).ge.qcrit)                    &  
     &go to 320
  310 continue
!.....no q=qcrit surface in plasma
      l = 0
      f1saw = 0
      f2saw = 1
      go to 321
!.....qcrit lies between qprof2(l) and qprof2(l+1)
  320 continue
      denom = qprof2(l+1) - qprof2(l)
      f1saw = (qprof2(l+1) - qcrit) / denom
      f2saw = (qcrit - qprof2(l)  ) / denom
!
!.....special logic for disruption model
      if(times.gt.acoef(95).and.1._R8/adi(npsit).lt.acoef(96)) then
      l = npsit-3
      f1saw = 0._R8
      f2saw = 1._R8
      endif
 321  lsaw = l+1
 
!
!
!
      go to(1001, 1002, 1003), isaw
 1001 continue
      if(lsaw.le.1) return
      etasaw = f1saw*etpara(lsaw-1)+f2saw*etpara(lsaw)
      do 401 j=2,npsit+1
      if(j-1 .le. lsaw)  then
       etpara(j-1) =etasaw*(1._R8-acoef(120))+acoef(120)*etpara(j-1)
      if(acoef(120).gt.1.0_R8) etpara(j-1) = acoef(120)*etasaw
      if(etpara(j-1) .gt. etav) etpara(j-1) = etav
                 endif
      fsaw = 1._R8
      tfluxs = (lsaw -1)*dpsi
      tflux  = (j-1)*dpsi
!     if(tflux .le. tfluxs) fsaw = acoef(124)
      if(tfluxs .gt. 0.0_R8) then
      if(tflux .eq. (lsaw+3)*dpsi) fsaw = 0.025_R8*acoef(124) +          &  
     & 0.975_R8
      if(tflux .eq. (lsaw+2)*dpsi) fsaw = 0.05_R8*acoef(124) + 0.95_R8
      if(tflux .eq. (lsaw+1)*dpsi) fsaw = 0.1_R8*acoef(124) + 0.9_R8
      if(tflux .eq. (lsaw)*dpsi  ) fsaw = 0.3_R8*acoef(124) + 0.7_R8
      if(tflux .eq. (lsaw-1)*dpsi) fsaw = 0.5_R8*acoef(124) + 0.5_R8
      if(tflux .eq. (lsaw-2)*dpsi) fsaw = 0.7_R8*acoef(124) + 0.3_R8
      if(tflux .eq. (lsaw-3)*dpsi) fsaw = 0.9_R8*acoef(124) + 0.1_R8
      if(tflux .lt. (lsaw-3)*dpsi) fsaw = acoef(124)
      dse1(j) = dse1(j)*fsaw
      dse2(j) = dse2(j)*fsaw
      dse3(j) = dse3(j)*fsaw
!
      dsi1(j) = dsi1(j)*fsaw
      dsi2(j) = dsi2(j)*fsaw
      dsi3(j) = dsi3(j)*fsaw
!
      cs1(j)  = cs1(j)*fsaw
      endif
  401 continue
      return
!
!--------------------------------------------------------------------
!       isaw = 2    ...  discrete sawtooth triggered by times on the
!                        type 75 input card
!--------------------------------------------------------------------
 1002 continue
!
!.....check if time is during a sawtooth crash
      hypermult = 0._R8
      if(numsaw.le.0) return
      crashtimem = max(crashtime,20._R8*dts)
      do 501 is = 1,numsaw
      issave = is
      if(times .ge. sawtime(is) .and.                                    &  
     &   times .le. sawtime(is) + crashtimem) go to 502
  501 continue
      if(use_transp .and. from_xplasma .and. saw_on) call tr_saw_off
      saw_on = .false.
      return
!
!.....sawtooth crash
!
  502 continue
      nskipsf = 1
!
      if(lmix(issave) .eq. 0 .and. lsaw.gt.1) then
!
!.......compute mixing region for Kadomtsev reconnection
        sum = 0
        do 330 ii=2,npsit-1
        iisave = ii
        qave = 0.5_R8*(qprof2(ii)+qprof2(ii-1))
        sumo = sum
        sum = sum + (1._R8/qave - 1._R8)
        if(sum .lt. 0 .and. sumo .gt. 0) go to 340
 330    continue
        lmix(issave) = 1
        go to 350
 340    lmix(issave) = max(lsaw,iisave)
 350    continue
        psisaw = xsv2(lmix(issave))
        r1mix = rminora(lmix(issave))
        hyperfrac = (psilim - psisaw)/(psilim - psimin)
        if(lmix(issave).le.1) hyperfrac = 1._R8
        hyperfracc= 1._R8- hyperfrac
!
!
        write(nterm,6402) issave,kcycle,lsaw,lmix(issave),               &  
     &                    times, dts, crashtimem,hyperfracc
        write(nout,6402) issave,kcycle,lsaw, lmix(issave),               &  
     &                    times, dts, crashtimem,hyperfracc
 6402 format(" Sawtooth #",i3," cycle, lsaw, lmix, time",i7,2i4,1pe10.3,  &  
     &       " dt="1pe9.2, " crashtimem=",1pe9.2,                        &  
     &       " hyperfracc",1pe9.2)
                                   endif
!
      hypermult = acoef(64)*gzero**2*r1mix**4/(crashtimem*xmag**2)
      transmult = acoef(124)
      if(lmix(issave).le.1) return
!
!     turn-on sawteeth in TRANSP
      saw_on = .true.
      saw_t0 =  sawtime(issave)
      saw_t1 =  saw_t0+crashtimem
      if(use_transp .and. from_xplasma .and. .NOT. saw_tr)                      &
     &                                          call tr_saw_on

      do 402 j=2,lmix(issave)
      gps = gja2(j)/vp2(j)*(tpi*qprof2(j))
      coef3 = - transmult*r1mix**2/crashtimem*udst*gps
      dse3(j) =  coef3
      dsi2(j) =  coef3
      dsi3(j) = -coef3
!
      rmid = .5_R8*(adn(j+1)/vp(j+1)+adn(j)/vp(j))
      emid = .5_R8*(ade(j+1)/vpg(j+1)+ade(j)/vpg(j))
      pmid = .5_R8*(adp(j+1)/vpg(j+1)+adp(j)/vpg(j))
      temidd = emid/rmid
      timidd = (pmid-emid)/rmid
      if(timidd.le.0) timidd = 0.1_R8*temidd
      dsi1(j) = -coef3*timidd
      dse1(j) = -coef3*temidd
!
      dsi0(j) = 0._R8
      dse0(j) = 0._R8
      cs1(j) = 0._R8
  402 continue
      return
!--------------------------------------------------------------------
!.....isaw=3    Porcelli model
!--------------------------------------------------------------------
!
 1003 continue
      if(lsaw.le.3) go to 1002
!
      rminor1 = f1saw*rminora(lsaw-1)+f2saw*rminora(lsaw)
      rmajor1 = f1saw*rmajora(lsaw-1)+f2saw*rminora(lsaw)
      elong1  = f1saw*elonga (lsaw-1)+f2saw*elonga (lsaw)
      delta1  = f1saw*deltaa (lsaw-1)+f2saw*deltaa (lsaw)
      r1 = rminor1*sqrt(elong1)
      s1 =(f1saw*(qprof2(lsaw)-qprof2(lsaw-2))/                          &  
     &           (rminora(lsaw)-rminora(lsaw-2))                         &  
     &   + f2saw*(qprof2(lsaw+1)-qprof2(lsaw-1))/                        &  
     &           (rminora(lsaw+1)-rminora(lsaw-1)))*r1
      aave = rminora(npsit)*sqrt(elonga(npsit))
      deltaq = 1._R8- qprof2(1)
      pi0 = (adp(2)-ade(2))/vpg(2)
      bp1sq = f1saw*gxmja2(lsaw-1)/(vp2(lsaw-1)*(tpi*qprof2(lsaw-1)))    &  
     &    +   f2saw*gxmja2(lsaw  )/(vp2(lsaw  )*(tpi*qprof2(lsaw  )))
      pleft = .5_R8*(adp(lsaw-1)/vpg(lsaw-1)+adp(lsaw)/vpg(lsaw))
      pright= .5_R8*(adp(lsaw+1)/vpg(lsaw+1)+adp(lsaw)/vpg(lsaw))
      prsaw = f1saw*pleft + f2saw*pright
      vsaw = f2saw*vp(lsaw)*dpsi
      sum = 0.5_R8*(pleft-prsaw)*vsaw
      sumv= vsaw
      sumcp = (.5_R8*(rminor1+rminora(lsaw-1)))**(1.5_R8)                &  
     &        *(adp(lsaw)-ade(lsaw))/(vpg(lsaw)*pi0)                     &  
     &        *(rminor1-rminora(lsaw-1))/rminor1
      sumal = 0._R8
      do j=2,lsaw-1
        rminorav = .5_R8*(rminora(j-1)+rminora(j))/rminor1
        sum   = sum + (adp(j)/vpg(j)-prsaw)*vp(j)*dpsi
        sumv  = sumv+ vp(j)*dpsi
 
        sumcp = sumcp + rminorav**(1.5_R8)                               &  
     &        *((adp(j)-ade(j))/(vpg(j)*pi0))                            &  
     &        *(rminora(j)-rminora(j-1))/rminor1
        sumal = sumal + rminorav**(1.5_R8)                               &  
     &        *.5_R8*(alphapr(j+1)-alphapr(j-1))
      enddo
      betap1 = (2._R8/bp1sq)*(sum/sumv)
      call geval(psimin,2,gval,gpval,gppval,imag,jmag)
      betai0 = 2._R8*pi0*xmag**2/gval**2
      betaalpha = -(2._R8/bp1sq)*sumal
      etparmks = etpara(2)*udsr
      valphen = (gval/xmag)/                                             &  
     &           sqrt(anhy(1)*2.5_R8*1.6726E-27_R8*4._R8*pi*1.E-7_R8)
      rhoi = 1.43E-4_R8*sqrt(ti(lsaw))*xmag/gval
      pileft = .5_R8*((adp(lsaw-1)-ade(lsaw-1))/vpg(lsaw-1)              &  
     &       +     (adp(lsaw)  -ade(lsaw)  )/vpg(lsaw))
      piright= .5_R8*((adp(lsaw+1)-ade(lsaw+1))/vpg(lsaw+1 )             &  
     &       +     (adp(lsaw)  -ade(lsaw)  )/vpg(lsaw))
      anleft = .5_R8*(adn(lsaw-1)/vp(lsaw-1)+adn(lsaw)/vp(lsaw))
      anright= .5_R8*(adn(lsaw+1)/vp(lsaw+1)+adn(lsaw)/vp(lsaw))
      pisaw = f1saw*pileft + f2saw*piright
      betai1 = 2._R8*pisaw*xmag**2/gval**2
      rp = abs((rminora(lsaw)-rminora(lsaw-1))                           &  
     &   * .5_R8*(pright+pleft)/(pright-pleft))
      rn = abs((rminora(lsaw)-rminora(lsaw-1))                           &  
     &   * .5_R8*(anright+anleft)/(anright-anleft))
      tau = te(2)/ti(2)
      omegastari = 0.989_R8*ti(2)*xmag/(gval*r1*rp)
      omegastara = 8.126E5_R8/ (gval * r1)
      gammap =0.8E3_R8*(xmag*r1/12.96_R8)**(-6._R8/7._R8)*(gval/(5.7_R8*  &  
     & xmag))**(2._R8/7._R8)                                             &  
     &                * (ane(2)/1.E20_R8)**(-3._R8/7._R8)*(te(2)/        &  
     & 2.E4_R8)**(1._R8/14._R8)                                          &  
     &                *s1**(6._R8/7._R8)
!
!
!
      call porcelli( r1, elong1, delta1, aave, s1 , deltaq,              &  
     &    betap1, sumcp, betai0, betaalpha, etparmks,                    &  
     &    valphen, xmag, rhoi, betai1, rp, rn, tau,                      &  
     &    omegastari, valphen, omegastara, gammap,                       &  
     &    dwcore, dwcorec, dw, dwc, dwc2, ratio ,                        &  
     &    scrit, dwbussac,dwelong,dwko,dwfast,dwideal)
!
!
!
!.....temporary diagnostic coding
      lprintdiag = lprintdiag + 1
      if(mod(lprintdiag,500) .eq. 1) then
        write(nout,6501) kcycle, times, lsaw, qprof2(1)
 6501   format(6x, "sawdiag....kcycle",i7," time=",1pe10.2,              &  
     &         " lsaw =",i5,"  q_0 = ",1pe12.4 )
        write(nout,6503) dwcore,dwcorec,dw,dwc,dwc2,ratio,               &  
     &                   dwbussac,dwelong,dwko,dwfast,                   &  
     &                   s1,scrit,betap1
 6503   format(9x," porcelli",1p6e12.4,/,9x," porcell2", 1p7e12.4)
      endif
!
!......test if sawtooth crash is in progress.  If not, test if
!      criteria for new crash is satisfied
!
      if(numsaw.le.0) go to 610
      crashtimem = max(crashtime,20._R8*dts)
      do 609 is=1,numsaw
      issave = is
      if(times .ge. sawtime(is) .and.                                    &  
     &   times .le. sawtime(is) + crashtimem) go to 602
  609 continue
      go to 610
  602 continue
!
!.....sawtooth in progress
      go to 1002
!
  610 continue
!
!----------------------------------------------------
!.....check 4 conditions for new sawtooth
!----------------------------------------------------
!
!.....condition #1...q0 < qsaw
!
      if(qprof2(1) .le. qsaw) then
        numsaw = numsaw + 1
        if(numsaw.ge.psaw) then
          ineg=57
          return
        endif
        sawtime(numsaw) = times
        write(nout,1701) numsaw, times, qprof2(1)
 1701   format(" Sawtooth #",i3,": time" , 1pe10.3,                      &  
     &         " trigger: q0=",1pe9.2," < qsaw")
        go to 1002
      endif
!
!.....condition #2...(Porcelli Eq. 13)
      if(-dwcore .gt. dwcorec) then
        numsaw = numsaw + 1
        if(numsaw.ge.psaw) then
          ineg=57
          return
        endif
        sawtime(numsaw) = times
        write(nout,1702) numsaw, times,-dwcore, dwcorec
 1702   format(" Sawtooth #",i3,": time" , 1pe10.3,                      &  
     &  " trigger: -dwcore =", 1pe9.2," > ", 1pe9.2," Eq.(13)" )
        go to 1002
      endif
!
!.....condition #3...(Porcelli Eq. 14)
      if(-dw .gt. dwc) then
        numsaw = numsaw + 1
        if(numsaw.ge.psaw) then
          ineg=57
          return
        endif
        sawtime(numsaw) = times
        write(nout,1703) numsaw, times,-dw, dwc
 1703   format(" Sawtooth #",i3," at time" , 1pe10.3,                    &  
     &  " trigger: -dw =",1pe9.2," > ",1pe9.2," Eq.(14)")
        go to 1002
      endif
!
!.....condition #4...(Porcelli Eq. 15)
      if(-dw .gt. dwc2 .and. ratio.lt.1) then
        numsaw = numsaw + 1
        if(numsaw.ge.psaw) then
          ineg=57
          return
        endif
        sawtime(numsaw) = times
        write(nout,1704) numsaw, times,-dw, dwc2, ratio
 1704   format(" Sawtooth #",i3," at time" , 1pe10.3,                    &  
     &  " trigger: -dw = ",1pe9.2," > ", 1pe9.2, " &",                   &  
     &         " ratio=",1pe9.2," < 1.0  Eq.(15)")
        go to 1002
      endif
      go to 1002
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
