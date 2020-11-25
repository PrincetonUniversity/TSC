!#include "f77_dcomplx.h"
!
!     -----------------------------------------------------------------
!
      SUBROUTINE GetKscat( kri, kzi, kphi, Br, Bz, Bph, tau,             &  
     &                     kro, kzo, kpho, thi, tho )
!     Get K scattered by fluctuations.
!     Copyright 1993 by F. W. Perkins, H. Takahashi, D. W. Ignat, & E. J. Valeo
!     Given the initial k in r,z,phi space and the B at that point
!     return the out k
!     by transforming to a space in grad psi, transverse, and parallel
!     referred to as                  rr    ,   tt      ,     ll
!     rotating the k randomly around the ll (field) direction thus preserving
!     k_parallel (kll)  and k_perp (kperp) but not k-phi
!     The parameter tau controls the randomness.
!     If tau is large then the scattered k peaks perpendicular to the wall
!     If tau is small, then the reflection is specular, without randomization
!
!     The angle is measured with respect to the tt direction
 
!      tt ( perp to B and to grad psi )
!      ^
!      |  theta(i)
!      |          .
!      |     .                           (parallel to B directionout of paper)
!      | .
!      +------- >  rr (grad poloidal flux
!
!     Following Perkins the Prob(ability) of scattering from thetai to theta is
!
!          Prob ~ sin(theta) exp{ - (theta-thetai)^2/tau }
!     and one forms the integral
!          ProbIntl(theta,thetai,tau) which is zero at theta=0 and 1 at theta=PI
!     so that if a random number is chosen between 0 and 1, a unique angle
!     is determined.  If tau is big, PI/2 is most likely; if tau is small
!     then thetai is most likely
!
!
!     IMPLICIT NONE
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*40 MyString
      INTEGER i,j,  isimp, ifirst, NumThets
      INTEGER ithetai, itheta
      INTEGER idum
      REAL*8    kri, kzi, kphi
      REAL*8    kro, kzo, kpho, thi, tho
      REAL*8    krr, ktt, kll, kperp, signth
      REAL*8    ktoti, ktoto, ktott
      REAL*8    Br, Bz, Bph, Bp, B, tau
      REAL*8    TAUMAX, TAUMIN, tauold, thetai, dtheta
      REAL*8    cosalph, sinalph, cosbeta, sinbeta, simp(2)
      REAL*8    PI, PIby2, SMALL, ERR
      REAL*8    random
      REAL*8    ran3
      PARAMETER (NumThets=101)
      PARAMETER (TAUMAX=100.0_R8, TAUMIN=1.0E-05_R8)
!     if  exp{ - (th_out_deg-th_in_deg)^2/del_th_deg^2 }
!     and exp{ - (th_out_rad-th_in_rad)^2/tauFWPerkins }
!     then    (del_th_deg,tauFWP) = (0.1deg,2.e-6) (1.0deg,3e-4) (360deg,40.0)
      REAL*8    ProbIntl(NumThets, NumThets),theta(NumThets)
      REAL*8    clevelin(20), tauDEG
      DATA    ifirst, SMALL, ERR / 1, 1.0E-25_R8, 1.0E-5_R8/
      DATA    simp(1), simp(2) / 4.00_R8, 2.00_R8/
 
      REAL*8    RE41
      REAL*8                                                             &  
     &        ZERO, ONE
      DATA    ZERO, ONE/                                                 &  
     &         0.0_R8, 1.0_R8/
      REAL*8 AREAL
 
!     If first call, set up PI and theta array
!     Also if first call set the ProbIntl to the result for large tau
      if (ifirst .eq. 1) then
        ifirst = 0
        PIby2  = asin(1.00_R8)
        PI     = PIby2 + PIby2
        dtheta = PI / AREAL(NumThets - 1)
        do 10 i = 1,NumThets
          theta(i) = dtheta * AREAL(i-1)
            ProbIntl(i,1) = 0.5_R8* ( 1._R8- cos (theta(i)) )
 
          do  9 j = 2,NumThets
            ProbIntl(i,j) = ProbIntl(i,1)
 9        continue
 10     continue
        tauold = TAUMAX
      endif
!
!
!     Fill the array of ProbIntl unless the values from the last call
!     or from the first call set up are ok
!
      if ( tau .ne. tauold .and.                                         &  
     &   (tau .le. tauold .or. tauold .ne. TAUMAX) ) then
 
        do 20 j=1,NumThets
          ProbIntl(1,j) = ZERO
 20     continue
 
        do 25 j=1,NumThets
          thetai = theta(j)
          do 24 i=2,NumThets
            isimp = mod(i,2) + 1
            ProbIntl(i,j) = ProbIntl(i-1,j) +                            &  
     &       simp(isimp)*sin(theta(i))*exp(-(theta(i)-thetai)**2/tau)
 24       continue
 25     continue
 
        do 30 j=1,NumThets
        do 30 i=2,NumThets
          ProbIntl(i,j) =  ProbIntl(i,j)/ProbIntl(NumThets,j)
 30     continue
 
!ccc  Contour plot begins
!     ------------------------------------------------------------------
!     kclev1: >0 --> clevelin contains levels to be used;
!                    dots for index less than kclev2;
!                    solid for index greater or equal kclev2
!     kclev1: =0 --> first contour at clevelin(1) with
!                    next one up by clevelin(2), and so on and so on
!     kclev1: <0 --> rcontr to choose -kclev1 equally spaced values between
!                    clevelin(1) and clevelin(2)
!     clevelin:      array of contour levels; this is output if kclev1<0
!     kclev2:        separates dots from solid lines
!     call EZrcon(ix1(1),ix2(1), jy1(1), jy2(1),
!    ^             xa,    xb,     ya,     yb   ,
!    ^  kclev1,clevelin,kclev2,
!    ^  PsiContr,
!    ^  i1stdim,
!    ^  xAry, ixmin, ixmax, ixstep,
!    ^  yAry, jymin, jymax, jystep )
      do i=1,9
      clevelin(i)=0.1_R8*AREAL(i)
      enddo
      call EZrcon(150,   500,  250,    600,                              &  
     &            0.0_R8,   3.2_R8,  0.0_R8,    3.2_R8,                  &  
     &  9,clevelin,5,                                                    &  
     &  ProbIntl,                                                        &  
     &  NumThets,                                                        &  
     &  theta, 1, NumThets-1, 1,                                         &  
     &  theta, 1, NumThets-1, 1)
      call EZwrit( 150, 150,                                             &  
     &            'ProbIntl; theta_in ordinate$',0,0)
      call EZwrit( 150, 125,                                             &  
     &            'theta_out abcsissa(radians)$',0,0)
      tauDEG = sqrt(tau)*180._R8/3.14_R8
      write(MyString,'(''tauDEG, tauFWP: '',                             &  
     &                   f5.1,1x,1pe9.2,''$'')')tauDEG,tau
      call EZwrit( 150, 100,MyString,0,0)
      call EZfini(0,0)
      call MkGrfLst   (' ProbIntl plot ')
!ccc  Contour plot ends
      endif
!     .
!     .                                 Here is the normal starting point !!
!     .
      Bp  = sqrt(Br*Br + Bz*Bz) + SMALL
      B   = sqrt(Br*Br + Bz*Bz + Bph*Bph) + SMALL
      cosalph = abs(Bph)/B
      sinalph =     Bp  /B
      cosbeta =     Br  /Bp
      sinbeta =     Bz  /Bp
 
      krr =         sinbeta*kri -         cosbeta*kzi +    ZERO*kphi
      ktt = cosalph*cosbeta*kri + cosalph*sinbeta*kzi - sinalph*kphi
      kll = sinalph*cosbeta*kri + sinalph*sinbeta*kzi + cosalph*kphi
 
      kperp = sqrt(krr*krr + ktt*ktt)
      ktott = sqrt(krr*krr + ktt*ktt + kll*kll)
      if ( krr  .lt. ZERO ) then
        signth = - 1.00_R8
        krr    = abs(krr)
        thetai = asin (krr/kperp)
        if (ktt .lt. ZERO ) thetai = PI - thetai
      else
        signth = + 1.00_R8
        thetai = asin (krr/kperp)
        if (ktt .lt. ZERO ) thetai = PI - thetai
      endif
      RE41 = thetai/dtheta
      ithetai = int(RE41) + 1
!     ithetai = ifix(thetai/dtheta) + 1
      thi = theta(ithetai)
!
      if (tau .lt. TAUMIN) then
        krr = -signth*krr
        ktt = +       ktt
        tho = thi
      else
!
        idum=1
        random = ran3(idum)
!
        do 50 i=2,NumThets
          itheta = i
          if (ProbIntl(i,ithetai) .gt. random ) go to 51
 50     continue
 51     continue
 
 
        krr = -signth*kperp* sin(theta(itheta))
        ktt =         kperp* cos(theta(itheta))
        tho = theta(itheta)
 
      endif
 
      kro =  sinbeta*krr + cosalph*cosbeta*ktt + sinalph*cosbeta*kll
      kzo = -cosbeta*krr + cosalph*sinbeta*ktt + sinalph*sinbeta*kll
      kpho=     ZERO*krr - sinalph*        ktt + cosalph*        kll
 
      ktoti = sqrt(kri*kri + kzi*kzi + kphi*kphi)
      ktoto = sqrt(kro*kro + kzo*kzo + kpho*kpho)
!
!     Take this out June 2000, below
!      if ( abs(ktoti-ktoto) .gt. ERR*ktoti .or.
!     ^     abs(ktoti-ktott) .gt. ERR*ktoti .or.
!     ^     abs(ktott-ktoto) .gt. ERR*ktoto  ) then
!c       call LSCpause
!        write(6,'(''GetKscat error; ktoti, ktoto, ktott:'',
!     ^                           3(1x,1pe10.3))') ktoti, ktoto, ktott
!c       call LSCpause
!      endif
!     Take this out June 2000, above
!
      tauold = tau
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
