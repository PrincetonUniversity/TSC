 
 
!     FTION  function for integration       ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE ftion(BoundsEr)
      USE dielec
      USE Doflags
      USE MKSetc
      USE params
      USE PIetc
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL plasma2d, eps, epsdr
!     KdB   k \cdot B = ( k_r B_r + k_z B_z + (n/r) B_{\phi} )
!     KdK   k \cdot k = ( k_r^2   + k_z^2   + (n/r)^2 )
!     \p = \partial , the Greek for partial derivative
!     \R = {\bf r}  , the vector r
!     \K = {\bf k}  , the vector k
!     \abs{#1} = \mid #1 \mid , the absolute value
!     s  =  s       , the abs value of \R; ds = path length
!     D  = D(\R,\K,\omega) ; D = 0 is the dispersion relation
!     D  = Kper4 * Eper
!        + Kper2 * [ Qpar*(Epar+Eper) + woc2*Exy**2 ]  ! Kper2 * bb
!        +  Epar * [(Qpar )**2 - (woc2*Exy)**2 ]       ! Epar  * cc
!     bkb(3) = 0.5 \p k_{\parallel}^2 / \p (k_r,k_z,n)
!     per(3) = 0.5 \p k_{\perp}^2     / \p (k_r,k_z,n)
!     dRdwt  = -\frac{\p D/\p\K}{\omega \p D/\p\omega}
!     dKdwt  = +\frac{\p D/\p\R}{\omega \p D/\p\omega}
!     dsdwt  = \abs{dRdwt}
!     dRds   = dRdwt/dsdwt      == f(1,2,3)
!     dKds   = dKdwt/dsdwt      == f(4,5,6)
!     dwtds  = 1./dsdwt
!     wdDdw  = \omega \p D/\p \omega
!     DKpar ==  \p D/\p k_{\parallel}^2
!     DKper ==  \p D/\p k_{\perp}^2
!     Qpar = Kpar2 - woc2*Eper
!     QparE= \p Qpar/\p Kper2
!     QparA= \p Qpar/\p Kpar2
!     QparW= \p Qpar/\p \omega^2  \cdot \omega^2
 
      INTEGER BoundsEr
      REAL*8    DispRela
      REAL*8    KdB,KdK,  Kper2,Kper4,Kpar2
      REAL*8    Btot2, Qpar, QparE, QparA, QparW, bb, cc
      REAL*8    bkb(3), per(3), dRdwt(3), dKdwt(3),  dsdwt
      REAL*8    DR,DZ
      REAL*8    r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
 
 
 
 
!     DATA DR,DZ/-0.0005,0.0005/
      DATA DR,DZ/-0.0001_R8,0.0001_R8/
 
      r = Y(1)
      z = Y(2)
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      BoundsEr = 0
            If ( pe2 .le. pe2min     ) then
               BoundsEr = 1
!
!                 write(nLSCcom2,'(/,/,'' pe2 '', g11.4,1x,
!     ^           ''.le.'',1x, g11.4,1x, ''in ftion '')')pe2, pe2min
!
               return
            endif
 
      Btot2 =  Br*Br + Bz*Bz + (RBphi/r)**2
      KdK = Y(4)**2 + Y(5)**2 + (Y(6)/R)**2
 
      KdB = Y(4)*Br + Y(5)*Bz + Y(6)/R/R *RBphi
      Kpar2 = KdB**2/Btot2
      Kper2 = KdK - Kpar2
      call eps   ( r, z, Kpar2, Kper2 )
      call epsdr ( r, z, Kpar2, Kper2 )
      bkb(1) = Br*KdB/Btot2
      bkb(2) = Bz*KdB/Btot2
      bkb(3) = rBphi*KdB/Btot2 /R/R
 
      per(1) = Y(4) - bkb(1)
      per(2) = Y(5) - bkb(2)
      per(3) = Y(6)/r**2 - bkb(3)
10    continue
      Kper4 = Kper2**2
      Qpar = Kpar2 - woc2*Eper
      QparE= -woc2*D11er
      QparA= 1._R8-woc2*D11ar
      QparW= -woc2*(Eper+D11w0)
 
      bb = Qpar*(Epar+Eper) + woc2*Exy*Exy
      cc = Qpar*Qpar - woc4*Exy*Exy
 
      denom = Kper4*D11w0                                                &  
     &     + Kper2*(QparW*(Epar+Eper) + Qpar*(D33w0+D11w0)               &  
     &     +        woc2*Exy*Exy + 2._R8*woc2*EXY*D12w0      )           &  
     &     + cc*D33w0                                                    &  
     &     + Epar*(2._R8*Qpar*QparW -2._R8*woc4*Exy*Exy                  &  
     &            -2._R8*woc4*Exy*D12w0                 )
      denom = 2._R8*denom
      wdDdw = denom
 
      DKper = 2._R8*Eper*Kper2 + bb                                      &  
     &      +  D11er *Kper4                                              &  
     &      + Kper2*(QparE*(Epar+Eper) + Qpar*(D33er+D11er)              &  
     &             + woc2*2._R8*Exy*D12er                      )         &  
     &       + Epar*(Qpar*QparE - woc4*Exy*D12er)*2._R8                  &  
     &       +D33er*cc
 
      DKpar = Kper4*D11ar                                                &  
     &     + Kper2*(QparA*(Epar+Eper) + Qpar*(D33ar+D11ar)               &  
     &     +       2._R8*woc2*Exy*D12ar   )                              &  
     &     + D33ar*cc                                                    &  
     &     + Epar*2._R8*(Qpar*QparA - woc4*Exy*D12ar )
 
      dRdwt(1) = -2._R8* ( DKper*per(1) + DKpar*bkb(1) )/denom
      dRdwt(2) = -2._R8* ( DKper*per(2) + DKpar*bkb(2) )/denom
      dRdwt(3) = -2._R8* ( DKper*per(3) + DKpar*bkb(3) )/denom
 
      dsdwt    = sqrt ( dRdwt(1)**2 + dRdwt(2)**2 + r*r*dRdwt(3)**2 )
 
      dKdwt(1) =(DispRela(Y(1)+DR,Y(2),     Y(4),Y(5),Y(6) ) -           &  
     &           DispRela(Y(1)-DR,Y(2),     Y(4),Y(5),Y(6) ) )           &  
     &                                                 /(2._R8*DR)/denom
      dKdwt(2) =(DispRela(Y(1),Y(2)+DZ,     Y(4),Y(5),Y(6) ) -           &  
     &           DispRela(Y(1),Y(2)-DZ,     Y(4),Y(5),Y(6) ) )           &  
     &                                                 /(2._R8*DZ)/denom
      dKdwt(3) = 0._R8
      f(1)     = dRdwt(1)/dsdwt
      f(2)     = dRdwt(2)/dsdwt
      f(3)     = dRdwt(3)/dsdwt
      f(4)     = dKdwt(1)/dsdwt
      f(5)     = dKdwt(2)/dsdwt
      f(6)     = 0._R8
      f(7)     = 1._R8/dsdwt
      dDdkABS  = wdDdw  * dsdwt
      return
      END
!                                                                      |
!                                                                      |
!     FTION ends                            ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
