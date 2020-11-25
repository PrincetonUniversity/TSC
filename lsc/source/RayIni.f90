!#include "f77_dcomplx.h"
 
 
!
!                                                                      |
!                                                                      |
!     AcDc                                  ---------------------------|
 
 
 
 
 
 
 
 
 
!     Rayini.F  begins                      ---------------------------|
 
!      Ini                                  ---------------------------|
!                                                                      |
!                                                                      |
 
      SUBROUTINE RayIni(RayIniErr)
      USE dielec
      USE Doflags
      USE params
      USE PlPr
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL plasma2d, eps, zplrcnr, ftion, prtout
!
! G77 Info explains:
!    The GNU Fortran language disallows `REAL(EXPR)' and `AIMAG(EXPR)',
! where EXPR is any `COMPLEX' type other than `COMPLEX(KIND=1)', except
! when they are used in the following way:
!
!      REAL(REAL(EXPR))
!      REAL(AIMAG(EXPR))
!
! The above forms explicitly specify that the desired effect is to
! convert the real or imaginary part of EXPR, which might be some `REAL'
! type other than `REAL(KIND=1)', to type `REAL(KIND=1)', and have that
! serve as the value of the expression.
!
!    The GNU Fortran language offers clearly named intrinsics to extract
! the real and imaginary parts of a complex entity without any conversion:
!
!      REALPART(EXPR)
!      IMAGPART(EXPR)
!
!    To express the above using typical extended FORTRAN 77, use the
! following constructs (when EXPR is `COMPLEX(KIND=2)'):
!
!      DBLE(EXPR)
!      DIMAG(EXPR)
!
! Ugly Complex Part Extraction
! ----------------------------
!
!    The `-fugly-complex' option enables use of the `REAL()' and `AIMAG()'
! intrinsics with arguments that are `COMPLEX' types other than
! `COMPLEX(KIND=1)'.
!
!    With `-ff90' in effect, these intrinsics return the unconverted real
! and imaginary parts (respectively) of their argument.
!
!    With `-fno-f90' in effect, these intrinsics convert the real and
! imaginary parts to `REAL(KIND=1)', and return the result of that
! conversion.
!
!    Due to this ambiguity, the GNU Fortran language defines these
! constructs as invalid, except in the specific case where they are
! entirely and solely passed as an argument to an invocation of the
! `REAL()' intrinsic.  For example,
!
!      REAL(REAL(Z))
!
! is permitted even when `Z' is `COMPLEX(KIND=2)' and `-fno-ugly-complex'
! is in effect, because the meaning is clear.
!
!    `g77' enforces this restriction, unless `-fugly-complex' is
! specified, in which case the appropriate interpretation is chosen and
! no diagnostic is issued.
! G77 Info explains, ends.
!
!     Integration is in coords  R     Z     Phi
!     Starting uses             Rad   Pol   Phi (no B in Rad direction)
!
!     Kr       1   **  Bz  Br  **  ** Krad **
!         =  ----  *            *  *        *
!     Kz     Bpol  ** -Br  Bz  **  ** Kpol **
!
!     Kpar = (KtorBtor + KpolBpol)/B
!     Krad2 = Kperp2 - (KpolBtor - KtorBpol)**2/B**2
!     Krad is positive for slow wave, negative for fast.
!     Bpol has the sign of (Bz(R-Rmaj) - BrZ) in order to keep the
!     sign of Kr, Kz, Krad straight.  This should work except in
!     quite pathological shapes.
!
!     NAG routine needs the aNAGco(4), ReZ(3), ImZ(3), Toleran,
!     and integers nNAG, Ifail
!     REAL
!    ^        aNAGco(4), ReZ(3), ImZ(3), toleran
!     INTEGER          nNAG     , ifail
!     REAL
!    ^        x02aaf
      REAL*8                                                             &  
     &        Kpar2, Kpol, Krad, Ktor
      COMPLEX*16 zc(3)
      REAL*8    Azplr(4)
      REAL*8                                                             &  
     &     r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      INTEGER i, DummyEr, RayIniErr
      REAL*8                                                             &  
     &        Bpol,Bpol2,Btot2,                                          &  
     &        DispRela,                                                  &  
     &        det,fast,plas,Qpar,Sig,slow,try,                           &  
     &        zi1,zi2,zi3,zr1,zr2,zr3, AprxSlow, AprxFast, discrimt
      REAL*8                                                             &  
     &        CosT,drad,rad,rstar,SinT,t, SEARCHINCR
      REAL*8 AREAL


      SEARCHINCR = 0.005_R8
      T = THET0 * 6.28318_R8
      RayIniErr = 0
                            Rad = Rmax - Rmaj
      If ( Cos(T) .lt. 0._R8) Rad = Rmaj - Rmin
      CosT = Abs(Cos(T)) + 1.0E-20_R8
      SinT = Abs(Sin(T)) + 1.0E-20_R8
      Rad  = min ( Rad/CosT , Zmax/SinT )
      rstar= -SEARCHINCR
      drad = abs(rstar)
 
5     do 30 I=1,10000
      Rad = Rad - drad
      r    = Rmaj + Rad*Cos(T)
      z    =        Rad*Sin(T)
      Kpol = enpol*2.0_R8*3.1416_R8*(fghz/0.3_R8)
      Ktor = enpar*2.0_R8*3.1416_R8*(fghz/0.3_R8)
      if (enpar .eq. 0.00_R8) then
        RayIniErr = 2
        call LSCwarn(' enpar == 0 in RayIni ')
        return
      endif
      call     plasma2d(r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
 
      If ( pe2 .le. pe2min ) Go To 30
!
      Bpol2= Br*Br + Bz*Bz
      Btot2= Bpol2 + (RBphi/r)**2
      Bpol = sqrt(Bpol2)
      Sig  = (r - Rmaj)*Bz - z*Br
      If ( Sig .lt. 0._R8) Bpol = -Bpol
      Kpar2 = ( Kpol*Bpol + Ktor*RBphi/R )**2/Btot2
      call Eps ( r, z, Kpar2, 0.0000_R8)
 
10    Azplr(4) = -(aio/fghz**4 + ael)
      Azplr(3) =  Eper
      Qpar = Kpar2 - woc2*Eper
      Azplr(2) = Qpar*(Epar+Eper) +  woc2*Exy**2
      Azplr(1) = Epar*( Qpar**2 - woc4*Exy**2 )
!
!     -                                 Estimate the roots far from lmc
      discrimt   =  Azplr(2)**2 - 4._R8* Azplr(3) * Azplr(1)
      if ( discrimt .le. 0._R8) then
        discrimt = 0._R8
        RayIniErr= 1
        call LSCwarn(' no accessibity found for this ray ')
        return
      endif
!
      call ZPLRCnr(3,Azplr,ZC)
!     Finds Zeroes of a Polynomial with Laguerre's method if
!     Real Coefficients
!     using Numerical Recipes code so we dont depend on IMSL
      zr1 = REAL (zc(1))
      zr2 = REAL (zc(2))
      zr3 = REAL (zc(3))
      zi1 = AIMAG(zc(1))
      zi2 = AIMAG(zc(2))
      zi3 = AIMAG(zc(3))
      if (PrFlg(RAYWR) .ge. TRUE) then
       AprxSlow = ( - Azplr(2) + sqrt(discrimt)) / (2._R8* Azplr(3))
       AprxFast = ( - Azplr(2) - sqrt(discrimt)) / (2._R8* Azplr(3))
       write(nLSCcom2,'('' Aprox Fast, Slow: '', 2(1x,1pe10.3),/,        &  
     &                  '' Roots found:      '', 3(1x,1pe10.3)   )')     &  
     & AprxFast, AprxSlow, zr1, zr2, zr3
       if (AprxSlow .gt. 0.0_R8) then
          write(nLSCcom2,'('' '',/)')
       endif
      endif
 
!cccccccccccccccccccccccccc
!     aNAGco(1)= Azplr(4)
!     aNAGco(2)= Azplr(3)
!     aNAGco(3)= Azplr(2)
!     aNAGco(4)= Azplr(1)
!     nNAG     = 4
!     toleran  = x02aaf(0.1)
!     call c02aef ( aNAGco, nNAG, ReZ, ImZ, toleran, ifail )
!     zr1      = ReZ(1)
!     zr2      = ReZ(2)
!     zr3      = ReZ(3)
!     zi1      = ImZ(1)
!     zi2      = ImZ(2)
!     zi3      = ImZ(3)
!cccccccccccccccccccccccccc
 
!     The normal case is that all roots are real, and we see if we can start.
22    continue
      plas = max(zr1,zr2,zr3)
      fast = min(zr1,zr2,zr3)
      slow = -1.0E+30_R8
      if (zr1.lt.plas .and. zr1.gt.fast) slow = zr1
      if (zr2.lt.plas .and. zr2.gt.fast) slow = zr2
      if (zr3.lt.plas .and. zr3.gt.fast) slow = zr3
      if (slow .eq. -1.0E+30_R8) go to 30
      if (lfast .eq. 1) try = fast
      if (lfast .eq. 0) try = slow
24    if (try .lt. 0._R8) go to 30
      try = try - ( Kpol*RBphi/R - Ktor*Bpol )**2/Btot2
      if ( try .lt. 0._R8) go to 30
      try  = +sqrt(try)
 
                       krad = + try
      if(lfast .eq. 1) krad = - try
 
      go to 35
30    continue
      if (PrFlg(RAYWR) .ge. TRUE ) write(nLSCcom2,600)                   &  
     &                                       zr1,zr2,zr3,fast,slow,plas
      RayIniErr =2
      call LSCwarn(' cant find a starting point for this ray')
      return
35    continue
      y(1)  =  r
      y(2)  =  z
      y(3)  =  0._R8
      y(4)  =  ( Krad*Bz + Kpol*Br ) / Bpol
      y(5)  =  (-Krad*Br + Kpol*Bz ) / Bpol
      y(6)  =    Ktor*r
      y(7)  =  0._R8
!     to begin the time as fn of distance
      y(8) = begin
      det = DispRela(y(1),y(2),    y(4),y(5),y(6))
!
      if (PrFlg(RAYWR) .ge. TRUE ) then
        write(nLSCcomm,602)r,z,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
        write(nLSCcomm,601)y(1),y(2),y(3),y(4),y(5),y(6),y(7),           &  
     &        det,d1,d2,d4,fast,slow,plas
      endif
!
70    call ftion(DummyEr)
!     only to initialize dkper & wdddw for damping calculation
      call prtout(0)                                                        
!
600   format (' error in ini',/,' roots are ',3e10.3,/,                  &
     &        ' fast,slow,plas are',3e10.3     )                           
601   format ( ' y is ',/,1x,1p7e10.2,/,' d(y) is ', 1pe10.2,/,          &
     &         ' D elements are:', 1p3e10.2,/,                           &  
     &         ' fast,slow,plas are' ,1p3e10.2)                            
602   format ( ' ',/,' New ray',/,                                       &
     &         '    r    z    psi     Br     Bz  RBphi',                 &  
     &         '    omc    Tee   pe2     pi2   aio     ael',/,           &  
     &              2(f5.2),6(f7.3),f7.0,3(f7.3))
!
      return
      END
!                                                                      |
!                                                                      |
!     ini ends                              ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
