!
!
!                                                                      |
!                                                                      |
!     File: grids.f(or) ends                ---------------------------|
 
 
 
 
 
 
 
 
!     jrf.F begins -----------------------------------------------------
!     .                                                                |
!     .                                                                |
!     Copyright 1991 by E. J. Valeo, C. F. F. Karney and N. J. Fisch,
!     Princeton University Plasma Physics Laboratory.
!     Ref: Current in wave-driven plasmas, by CFFK and NJF, Physics of
!     Fluids <29> pp 180-192 (1986).  See esp. eqn. 21b.
!     In constructing units, the rest of the code has velocities normalized
!     to c, and f contains a density in cm^{-3}.
!     Dql has a c^2 taken out of it.
!     E is in v/m we assume.
!     \Gamma   = n_e q^4     \ln \lambda / (4 \pi \epsilon_0^2 m^2)
!              = n_e q^4 c^4 \ln \lambda \mu_0/(4\pi) \mu_0 / m^2
!                 with units of velocity^3 / time
!     v_r      = -sign(qE) \sqrt{m\Gamma/\abs(qE)}
!                 where q has the sign of the electron charge = - 1.6e-19
!                 and the electrons run away in the negative direction if
!                 v_r is positive
!     v_{norm} = \abs(v_r) / c
!              = sqrt ( \Gamma / \abs(qE/m) ) / c
!              = sqrt ( n/E q^3/m c^2 lnLambda mu_0 (mu_0/4PI) )
!                 dimensionless, in units of c
!     gmrun    = n_e q^4 ln(Lambda) / ( 4 PI eps_0^2 m_e^2 )          /  c^3
!              = n_e q^4 c^1 / m_e^2 * lnLambda * (mu_0/(4PI))*mu_0
!                 with units of time^{-1}
!
!     nuRuna   = gmrun/vnorm^3
!                 with units of time^{-1}
!     NeAry in cm^{-3}
!     gmrun in s^{-1}
!     nuRuna in s^{-1}
!     js in ampere/m^2
!     dfdv in cm^{-3}
!     dlnJdlnE = dJ/dE * (E/J) is dimensionless
!     and can be shown for small E (E is the E_{dc}) to be approximately
!     eE_{dc}/(m Gamma) * v_phase^2 * (1/(2 \mu)) ( (2+Z+3\mu^2)/(3+Z) )
!     =~ eE_{dc} / lnLamba /n_\parallel^2 /
!                  [ e^2/(4\pi\epsilon_0 c^2/\omega_p^2)] * Order(1)
!
!     ------------------------------------------------------------------
!
!     Keeping the signs straight is a problem.  We say that positive
!     phase velocity tends to drive current in the supportive, or correct
!     direction.  Positive current drive.  Thus, presumably, the Edc
!     is also in the positive direction.  Electrons tend to run away
!     and be accelerated in the opposite direction.  Thus for positive
!     electron velocity which goes with positive phase velocity, the
!     normalized   u \equiv v/v_r is negative, and the mu for this positive
!     v is negative.
!
!     v  E      u    mu
!     +  +      -     -
!     -  +      +     +
!     +  -      +     +
!     -  -      -     -
!
!     Code is written with the presumtion of E positive and mu minus
!     for positive velocity and mu plus for negative velocity.  If
!     E is negative, the quantitiy muminus becomes +1; muplus becomes -1
!
!     ------------------------------------------------------------------
!
!     Getting jrf and the dc electric field to be self consistent
!     is a problem.  This fragment from TSC concerns the iteration:
!     .
!      subroutine curdrive ! TSC Code Fragment to write lhcdoua
!c
!c.....define quantities needed for current drive
!c.....(neutral beam and bootstrap current)
!c
!      include 'clinam.i'
!      include 'scr3.i'
!      include 'saprop.i'
!      data ifstrt/0/
!      data ifircd/0/
!c
!      dimension wrat(ppsi),anucoll(ppsi),tauconf(ppsi),disedge(ppsi)
!      dimension floss(ppsi),tauf(ppsi),form(ppsi),ajavlho(ppsi)
!      dimension djdets2(ppsi),djdelhn(ppsi)
!c
!      ifircd=0
!      if(ilhcd.eq.0 .or. ifk.ne.2 .or. ifircd .ne. 0) go to 49
!      iouttsc(1:6) = 'tscout'
!      iouttsc(7:7) = isuffix(1:1)
!      open(nlhcd2,file=iouttsc,status='unknown',iostat=ios2)
!      ilhcdou(1:6) = 'lhcdou'
!      ilhcdou(7:7) = isuffix(1:1)
!      open(nlhcd,file=ilhcdou,status='unknown',iostat=ios11)
!   49 continue
!c
!
!........
!c
!c     set counter keeping track of LSC/TSC iterations
!c
!      ICLH = 0
!c     re-cycle point if current is not well-converged
!c
! 1654 continue
!
!      call lhwrite
!      call lsc(bpowsmw,nlhcd,nlhcd2,nsc1,nout,nlhcdin,nterm,irayt
!     1        ,iplot,ierr)
!      iraytd = irayt
!      if(ierr.gt.0) ineg=45
!      if(ierr.lt.0) write(nout,1991) ierr
! 1991 format("LSC returned ierr=",i5)
!      if(ierr.lt.0) write(nterm,1991) ierr
!      if(ineg.ne.0) return
!      rewind nlhcd2
!      read(nlhcd2,1001) (powtsc(l),l=2,npsit)
!      read(nlhcd2,1001) (currtsc(l),l=2,npsit)
!      read(nlhcd2,1001) (djdetsc(l),l=2,npsit)
!      read(nlhcd2,1001) (djdets2(l),l=2,npsit)
! 1001 format(1p5e16.6)
!  651 continue
!c
!      powsum = 0.
!      cursum = 0.
!      do 653 l=2,npsi
!      powsum = powsum + powtsc(l)*vp(l)*dpsi
!      cursum = cursum + currtsc(l)*1.e4*vp(l)*dpsi/(tpi*xplas)
!  653 continue
!c     write(nout,9662)
!c     write(nout,9661) (powtsc(l),l=2,npsit)
!c     write(nout,9663)
!c     write(nout,9661) (currtsc(l),l=2,npsit)
!c     write(nout,9664)
!c     write(nout,9661) (djdetsc(l),l=2,npsit)
!c     write(nout,9665)
!c     write(nout,9661) (voltlp(l),l=2,npsit)
!c     write(nout,9666)
!c     write(nout,9661) (.5*udsv*(as(l)+as(l-1)),l=2,npsit)
! 9666 format(" udsv*as ")
! 9661 format(1p5e12.4)
! 9662 format(" powtsc")
! 9663 format(" currtsc")
! 9664 format(" djdetsc")
! 9665 format(" voltlp")
!      iplot=0
!      irayt=0
!c
!c
!      do 660 j=2,npsit
!      ajavlho(j) = ajavlh(j)
!      djdelh(j) = 0.
!      ajavlh(j) = currtsc(j)*1.e4*usdi
!c     if( abs(voltlp(j)) .ge. .01)
!c    1djdelh(j) = djdetsc(j)*ajavlh(j)*tpi / (voltlp(j)*usdv)
!      djdelh(j) = djdets2(j)*1.e4*usdi/(usdv)
!  660 continue
!      ajavlh(1) = ajavlh(2)
!      djdelh(1) = djdelh(2)
!      djdetsc(1) = djdetsc(2)
!      if(ifk.eq.1) go to 700
!c
!c     correct current by Newton's iteration
!c
!      AJAVLHM = 0.
!      do 710 j=2,npsit
!      aterm = .5*(etpara(j)+etpara(j-1))*xplas*djdelh(j)
!      if(aterm.lt.0) aterm=0.
!      fac = aterm / (1.+aterm)
!      ajnew = ajavlh(j) + fac*(ajavlho(j) - ajavlh(j))
!      ajavlh(j) = ajnew
!      AJAVLHM = max(AJAVLHM,abs(ajavlh(j)))
!      currtsc(j) = ajavlh(j) / (1.e4*usdi)
!  710 continue
!      ajavlh(1) = ajavlh(2)
!      currtsc(1) = currtsc(2)
!c     write(nout,9663)
!c     write(nout,9661) (currtsc(l),l=2,npsit)
!c
!      cursum2 = 0.
!      do 753 l=2,npsi
!      cursum2 = cursum2 + currtsc(l)*1.e4*vp(l)*dpsi/(tpi*xplas)
!  753 continue
!      write(nout,9660) kcycle,iraytd,bpowsmw,powsum,cursum,cursum2
!      write(nterm,9660) kcycle,iraytd,bpowsmw,powsum,cursum,cursum2
! 9660 format(" LSC:N=",i5," iraytd=",i2," ,powin,out,cur",
!     1     1p4e 9.1)
!c
!c     see if the correction to the current is reasonably big,
!c     and if so, rewrite the equilibrium file and recall LSC;
!c     but not more than 50 times
!      RATHLH = 0.
!      do 1770 j=2,npsit
!      DENOM = ABS(ajavlh(j))
!      IF(DENOM .LE. .01*AJAVLHM) GO TO 1770
!      RATHLH = max(RATHLH,ABS(ajavlh(j)-ajavlho(j))/DENOM)
! 1770 continue
!      write(nterm,4700) kcycle,iclh,rathlh
! 4700 format(i5,i3,1pe12.4)
!      IF(RATHLH .LE. 0.1) GO TO 700
!      ICLH = ICLH+1
!c
!c     under-relax to avoid oscillations in current
!c
!      write(nout,1790) kcycle,iclh
!      do 1780 j=2,npsit
!      DENOM = ABS(ajavlh(j))
!      IF(DENOM .LE. .01*AJAVLHM) GO TO 1780
!      diff = (ajavlh(j)-ajavlho(j))/DENOM
!      if(abs(diff).le.0.1) go to 1780
!      write(nout,1781) j,ajavlh(j),ajavlho(j),voltlp(j),djdelh(j)
!     1     ,diff
! 1780 continue
! 1790 format(" kcycle=",i4,"  iclh=",i2,/,
!     1"  j      ajavlh     ajavlho      voltlp      djdelh diff")
! 1781 format(i3,1p6e12.4)
!      backa = 0.5
!      if(iclh.gt.10) backa = .25
!      if(iclh.gt.15) backa = .10
!      if(iclh.gt.25) backa = .05
!      DO 1700, J=2, NPSIT
!      AJAVLH(J) = AJAVLHO(J) + backa*( AJAVLH(J) - AJAVLHO(J) )
! 1700 CONTINUE
!c
!      IF(ICLH .LT. 50) GO TO 1654
!      INEG = 45
!  700 continue
!      ajavlh(npsit+1) = ajavlh(npsit)
!      djdelh(npsit+1) = djdelh(npsit)
!      djdetsc(npsit+1)= djdetsc(npsit)
!c
!c.....define surface centered arrays
!      do 720 j=1,npsit
!      ajavlh2(j) = .5*(ajavlh(j)+ajavlh(j+1))
!      djdetsc2(j) = .5*(djdetsc(j)+djdetsc(j+1))
!  720 continue
!c
!      return
!c
!  701 continue
!      do 702 j=1,npsit+1
!      powtsc(j) = 0.
!      currtsc(j) = 0.
!      djdetsc(j) = 0.
!      ajavlh(j) = 0.
!      djdelh(j) = 0.
!      ajavlh2(j) = 0.
!      djdetsc2(j) = 0.
!  702 continue
!      return
!c
!      end
!
!     ------------------------------------------------------------------
!
      SUBROUTINE VmaxCalc
      USE FeBins
      USE Jrf
      USE params
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip,ir,iz,iv
      REAL*8    vl
      do  5   ip=1,NPSIDIM
        VparMaxN(ip)= -0.05_R8
        VparMaxP(ip)= +0.05_R8
 5    continue
!
      do 20   ir=1,nrays
        do 10 iz=1,nzones
        iv = ivind(iz,ir)
        ip = izind(iz,ir)
          if (iv .eq. 0 .or. ip .eq. 0) goto 15
        vl = vpar(iv)
 
 
        if(VparMaxP(ip) .lt. vl) VparMaxP(ip) = vl
        if(VparMaxN(ip) .gt. vl) VparMaxN(ip) = vl
 10   continue
 15   continue
 20   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
