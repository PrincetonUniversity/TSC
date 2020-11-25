 
!     FUNCTION DispRela                     ---------------------------|
!                                                                      |
!                                                                      |
      REAL*8 FUNCTION  DispRela (r, z, Kr, Kz, Kphi)
      USE dielec
      USE MKSetc
      USE params
      USE PIetc
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8    r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      REAL*8    Kr, Kz, Kphi, Kpar2, Kper2, KdK
      REAL*8    Btot2, Qpar
!
      KdK  = Kr**2 + (Kz)**2 + (KPhi/R)**2
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      Btot2 =  Br**2 + Bz**2 + (RBphi/r)**2
      Kpar2 = ( Kr*Br + Kz*Bz + Kphi*rBphi/r/r )**2/Btot2
      Kper2 = KdK - Kpar2
 
      call eps ( r, z, Kpar2, Kper2 )
      Qpar = Kpar2 - woc2*Eper
      D1 = Kper2**2*Eper
      D2 = Kper2*( (Epar+Eper)*Qpar + woc2*Exy**2 )
      D4 = Epar*( Qpar**2 - woc4*Exy**2 )
      DispRela  = D1  +  D2  +  D4
      return
      END
!                                                                      |
!                                                                      |
!     DispRela ends                         ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
