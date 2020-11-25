!#include "f77_dcomplx.h"
!
!     -----------------------------------------------------------------
!
!     AccesiSB.F                        -------------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE AccesiSB (r1,r2, z1,z2, Rmag)
!     .                                 Makes contour plot of psi surfaces
!     .                                 and accessibility surfaces
!     .                                 and n_\parallel enhancement surfaces
!     .                                 PsiContr AccContr EnhContr
!     AccContr uses the classic accessibility formula
!      n_//_for_accessibility >= sqrt(Epar) + omega_pe/omega_ce
!     EnhContr uses a form recently derived by Kupfer, Moreau, Ono, Takahashi
!     as follows:
!      n_// == n_tor Btor/Btot + n_pol Bpol/Btot;
!      n_tor = n_//_0 R_launch_point/R
!      n_pol <= n_perp = n_// omega_pe/omega/sqrt(Epar)!electrostatic estimate
!
!     n_// <= n_//_0  { (R_launch / R)  /
!                     [1 - (Bpol/Btot) (omega_pe/omega) / sqrt(Epar)] }
!     .                                 EnhContr == {}
      USE dielec
      USE params
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
      REAL*8    r1, r2, z1, z2, Rmag
      REAL*8  AREAL
 
!
!     Here are variables used in the contour plot:
      CHARACTER*30 MyString
      INTEGER kclev1, kclev2, i1stdim,j2nddim,                           &  
     &                                  ixmin, ixmax, ixstep,            &  
     &                                  jymin, jymax, jystep
      INTEGER nCntrDim,nPnts
      PARAMETER (i1stdim=50, j2nddim=50)
      PARAMETER (nPnts=i1stdim*j2nddim)
      PARAMETER (ixmin =1 ,jymin =1,  kclev2=1)
      PARAMETER (ixstep=1 ,jystep=1 )
      PARAMETER (nCntrDim=50)
      REAL*8    PsiContr(i1stdim,j2nddim)
      REAL*8    AccContr(i1stdim,j2nddim)
      REAL*8    EnhContr(i1stdim,j2nddim)
      REAL*8    Rpunct, Npunct
      REAL*8    xAry(i1stdim), xa, xb, xdumy
      REAL*8    yAry(j2nddim), ya, yb
      REAL*8    clevelin(nCntrDim), ContrMin, ContrMax
!     Here are variables for this subroutine
      REAL*8                                                             &  
     &                  r0,z0, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael,  &  
     &                  Bpol, Btot, BpolAry(i1stdim), BpolMax,           &  
     &                  psiNorm, psiEdge, psiCent, AccMax
!     REAL*8               r0,z0, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
!     REAL*8    Bpol, Btot, BpolAry(i1stdim), BpolMax
!     REAL*8    psiNorm, psiEdge, psiCent, AccMax
      INTEGER i, j, k
      INTEGER ix1(4), ix2(4), jy1(4), jy2(4)
!
!     Location reminder:
!     3 = up  left   4 = up  right
!     1 = low left   2 = low right
!     ------------------------------------------------------------------
      DATA    ix1(1),ix2(1),ix1(2),ix2(2),ix1(3),ix2(3),ix1(4),ix2(4) /  &  
     &          150,   400,   650,   900,   150,   400,   650,   900  /
      DATA    jy1(1),jy2(1),jy1(2),jy2(2),jy1(3),jy2(3),jy1(4),jy2(4) /  &  
     &           75,   325,    75,   325,   450,   700,   450,   700  /
 
      REAL*8                                                             &  
     &        ZERO, ONE, RE81, RE82
      DATA ZERO, ONE/                                                    &  
     &      0.0_R8, 1.0_R8/
 
!  round up from max found--> nice value; #major divns; #minor divns
      call EZrnd2 (r2,xb, i,j)
      call EZrnd2 ((xb-r1), xdumy, i,j)
      call EZrnd2 (z2,yb, i,j)
      if (2._R8*yb .gt. xdumy) then
        xa = xb - 2.0_R8*yb
        ya = - yb
      else
        xa = xb - xdumy
        yb = xdumy / 2._R8
        ya =-yb
      endif
 
!     Fill arrays at points where data is given
      do 10 i = 1, i1stdim
        xAry(i) = r1 + AREAL(i-1)*(r2-r1)/(AREAL(i1stdim)-1)
 10   continue
      do 11 j = 1, j2nddim
        yAry(j) = z1 + AREAL(j-1)*(z2-z1)/(AREAL(j2nddim)-1)
 11   continue
!
!     .                                 Get a normalization for psi
         call plasma2d(r2,ZERO,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
!        call plasma2d (r2,0., psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
         psiEdge = 0.98_R8*psi
        call plasma2d(Rmag,ZERO,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
!       call plasma2d(Rmag,0.,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
         psiCent = psi
         psiNorm = (psiEdge-PsiCent)
!
         AccMax = 1.0_R8
         BpolMax= 0.0_R8
      do 20 i = 1, i1stdim
          r0 = xAry(i)
!          radPaoletti(i) = xAry(i)
        do 19 j = 1, j2nddim
          z0 = yary(j)
!          zeePaoletti(j) = yAry(j)
 
          call plasma2d (r0,z0, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
!     .                                 If outside the plasma set nAcc to 1
!     .                                 otherwise find the accesible nAcc
          if (psi .ge. psiEdge) then
            PsiContr(i,j) = 1._R8
            AccContr(i,j) = 1._R8
            EnhContr(i,j) = 1._R8
            Bpol = 0.0_R8
!            rpsPaoletti(i,j) = sqrt(PsiContr(i,j))
!            ompPaoletti(i,j) = 0.00
!            BthPaoletti(i,j) = 0.00
!            BphPaoletti(i,j) = 0.00
          else
            PsiContr(i,j) = (psi-psiCent)/psiNorm
            call eps (r0, z0, 0._R8, 0._R8)
            AccContr(i,j) = sqrt( pe2/(omc*omc) ) + sqrt(Eper)
            Bpol = sqrt ( Br*Br + Bz*Bz )
            Btot = sqrt ( Bpol*Bpol + RBphi*RBphi/r0/r0 )
            EnhContr(i,j) = ( r0/r2 ) * ( 1._R8                          &  
     &                       -sqrt(abs(1._R8-Epar))*                     &  
     &                       Bpol/Btot/sqrt(abs(Eper))  )                  
            if(EnhContr(i,j) .le. 0.1_R8) EnhContr(i,j) = 0.1_R8
            EnhContr(i,j) = 1.00_R8/ EnhContr(i,j)
 
!            rpsPaoletti(i,j) = sqrt(abs(PsiContr(i,j)))
!            ompPaoletti(i,j) = sqrt( pe2/(omc*omc) )
!            BthPaoletti(i,j) = Bpol
!            BphPaoletti(i,j) = abs(RBphi/r0)
          endif
 
          if(  j .eq. (j2nddim+1)/2 ) then
            BpolAry(i) = Bpol
            if (BpolAry(i) .gt. BpolMax ) BpolMax=BpolAry(i)
          endif
          if (AccContr(i,j) .gt. AccMax) AccMax=AccContr(i,j)
!
 19     continue
 20   continue
      do 25 i = 1, i1stdim
        BpolAry(i) = BpolAry(i)/BpolMax
 25   continue
!     .                                 END of the preparations to plot
 
 
!     .                                 Upper Right Corner (index=4)
!     .                                 Potential n_// Enhancement Mid-Plane
!     .                                 Puncture plot of actual enhancement
!     .
      RE81 = 10._R8
      call EZsets(ix1(4),ix2(4), jy1(4), jy2(4),                         &  
     &             xa,    xb,   ZERO,    RE81,   1)
      call EZrnd2 (xb-xa,xdumy,i,j)
      call EZaxes(i,j,2,5)
      j = (j2nddim+1)/2
      call EZcurv(xAry(1),EnhContr(1,j),i1stdim)
      call EZwrit( ix1(4), (jy2(4)+50),                                  &  
     &            'nll enhance slice$',0,0)
      write(MyString,'(''Bpol/'',f6.3,''T ...$'')')BpolMax
      call EZwrit( ix1(4), (jy2(4)+25),                                  &  
     &            MyString            ,0,0)
      call EZcros(xAry(1),BpolAry(1),i1stdim)
!     .                                 Puncture plot begins
!     .                                 Go thru all the points looking for
!     .                                 a crossing of the midplane.  If found,
!     .                                 then put a point at that major radius
!     .                                 with the n// enhancement found there.
      do i = 2, izone
         if ( ZofRay(i)*ZofRay(i-1) .le. 0.00_R8) then
           Rpunct = 0.5_R8*(RofRay(i)+RofRay(i-1))
           Npunct = 0.5_R8*(npar(i,iray)+npar(i-1,iray)) /               &  
     &                   npar(1,iray)
           call EZpnts(Rpunct,Npunct,1)
         endif
      enddo
!     .                                 Puncture plot ends
 
!     .                                 Upper Right Corner ENDS
!     .
!     .
!     HELP on contour plot usages:
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
!
!     call rcontr(kclev1,clevelin,kclev2,  ! describes contour levels
!    ^  AryContr,                          ! the 2d array of values
!    ^  i1stdim,                           ! first dimension of Ary
!    ^  xAry, ixmin, ixmax, ixstep,        ! x and y locations, with ranges
!    ^  yAry, jymin, jymax, jystep )       ! of indexes to be used
!     ------------------------------------------------------------------
 
      ixmax = i1stdim
      jymax = j2nddim
!     .                                 Lower Left Corner (index = 1)
!     .                                 Psi contours
!     kclev2 =  1 ! by parameter
      ContrMin = 0.01_R8
      ContrMax = 0.99_R8
 
      kclev1 = 11
      do 30 k=1,kclev1
        clevelin(k) = ContrMin + AREAL(k-1)*(ContrMax-ContrMin)/         &  
     &                           AREAL(kclev1-1)
 30   continue
      call EZrcon(ix1(1),ix2(1), jy1(1), jy2(1),                         &  
     &             xa,    xb,     ya,     yb   ,                         &  
     &  kclev1,clevelin,kclev2,                                          &  
     &  PsiContr,                                                        &  
     &  i1stdim,                                                         &  
     &  xAry, ixmin, ixmax, ixstep,                                      &  
     &  yAry, jymin, jymax, jystep )
      call EZwrit( ix1(1), (jy2(1)+50),                                  &  
     &            'Flux contours$',0,0)
!     .                                 Lower Left Corner ENDS
!     .
!     .
!     .                                 Lower Right Corner (index = 2)
!     .                                 Accessibility Contours
      ContrMin = 1.250_R8
      kclev1 = 15
      do 40 k=1,kclev1
        clevelin(k) = ContrMin + AREAL(k-1)*0.25_R8
 40   continue
      call EZrcon(ix1(2),ix2(2), jy1(2), jy2(2),                         &  
     &             xa,    xb,     ya,     yb   ,                         &  
     &  kclev1,clevelin,kclev2,                                          &  
     &  AccContr,                                                        &  
     &  i1stdim,                                                         &  
     &  xAry, ixmin, ixmax, ixstep,                                      &  
     &  yAry, jymin, jymax, jystep )
      call EZwrit( ix1(2), (jy2(2)+50),                                  &  
     &            'nll access cntrs$',0,0)
      call EZwrit( ix2(2), (jy2(2)+50),                                  &  
     &            '1.25 strt$',0,0)
      call EZwrit( ix2(2), (jy2(2)+25),                                  &  
     &            '0.25 aprt$',0,0)
!     .                                 Lower Right Corner ENDS
!
!
!     .                                 Upper Left Corner (index = 3)
!     .                                 Potential n_// Enhancement Coutours
      kclev1 = 7
      do 50 k=1,kclev1
        clevelin(k) = ( sqrt(2.00_R8) )**AREAL(k)
 50   continue
 
 
      call EZrcon(ix1(3),ix2(3), jy1(3), jy2(3),                         &  
     &             xa,    xb,     ya,     yb   ,                         &  
     &  kclev1,clevelin,kclev2,                                          &  
     &  EnhContr,                                                        &  
     &  i1stdim,                                                         &  
     &  xAry, ixmin, ixmax, ixstep,                                      &  
     &  yAry, jymin, jymax, jystep )
      call EZwrit( ix1(3), (jy2(3)+50),                                  &  
     &            'nll enhance cntrs$',0,0)
      call EZwrit( ix2(3), (jy2(3)+50),                                  &  
     &            'rt 2 strt$',0,0)
      call EZwrit( ix2(3), (jy2(3)+25),                                  &  
     &            'rt 2 aprt$',0,0)
!     .                                 Upper Left Corner ENDS
!
!
      call EZfini(0,0)
!
      return
      END
!                                                                      |
!                                                                      |
!     Accesibi.F                        -------------------------------|
!
!     .                                                                |
!     .                                                                |
!     Rayio.F  ends                         ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
