!
      SUBROUTINE WhchRoot (r, z, Kr, Kz, Kphi, NperFs, NperSl)
      USE dielec
      USE MKSetc
      USE params
      USE PIetc
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8    NperFs, NperSl, AAD1, BBD2, CCD4, B2m4AC
      REAL*8    r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      REAL*8    Kr, Kz, Kphi, Kpar2, Kper2, KdK
      REAL*8    Btot2, Qpar
!
      REAL*8    R41
      REAL*8                                                             &  
     &        ZERO, ONE
      DATA    ZERO, ONE/                                                 &  
     &         0.0_R8, 1.0_R8/
!
      KdK  = Kr**2 + (Kz)**2 + (KPhi/R)**2
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      Btot2 =  Br**2 + Bz**2 + (RBphi/r)**2
      Kpar2 = ( Kr*Br + Kz*Bz + Kphi*rBphi/r/r )**2/Btot2
      Kper2 = KdK - Kpar2
 
      call eps ( r, z, Kpar2, Kper2 )
      Qpar = Kpar2 - woc2*Eper
      AAD1 = Eper
      BBD2 = ( (Epar+Eper)*Qpar + woc2*Exy**2 )/woc2
      CCD4 = Epar*( Qpar**2 - woc4*Exy**2 )/woc4
      B2m4AC = BBD2**2 - 4._R8*AAD1*CCD4
      NperFs = - BBD2 / (2._R8*AAD1)
      NperSl = NperFs
      if (B2m4AC .gt. 0._R8) then
        B2m4AC = sqrt(B2m4AC)/(2._R8*AAD1)
        NperFs = NperFs - B2m4AC
        NperSl = NperSl + B2m4AC
      endif
 
        NperFs = max(NperFs, ZERO)
        NperSl = max(NperSl, ZERO)
        NperFs = sqrt( NperFs)
        NperSl = sqrt( NperSl)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
