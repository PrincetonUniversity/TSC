!
!     -----------------------------------------------------------------
!
      SUBROUTINE BounceIt (yok, y, n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,i
      REAL*8 yok(n), y(n)
      REAL*8    r,z,   psi, Br, Bz, RBphi, omc, Tee, pe2, pi2, aio, ael
      REAL*8    B2, crr, crz, czr, czz, KrNew, KzNew
!     Go back to the last ok solution point, and re-arrange kr and kz
!     such that k \dot B is unchanged but k \cross B changes sign.  Then
!     put this new information into the starting condition for y.
!     y(1)  y(2)  y(3)  y(4)  y(5)  y(6)
!      r     z    phi   k_r   k_z    n
!
 
 
      r   = yok(1)
      z   = yok(2)
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
 
      B2  =  Br*Br + Bz*Bz
      crr = (Br*Br - Bz*Bz)/B2
      czz = - crr
      crz = 2._R8*  Br*Bz/B2
      czr = crz
 
      KrNew = crr * yok(4) + crz * yok(5)
      KzNew = czr * yok(4) + czz * yok(5)
 
      yok(4) = KrNew
      yok(5) = KzNew
      do 10 i=1,n
        y(i) = yok(i)
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
