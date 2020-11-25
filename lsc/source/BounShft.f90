!
!     -----------------------------------------------------------------
!
      SUBROUTINE BounShft (yok, y, n, ScatKdeg, thi, tho)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,i
      REAL*8    yok(n), y(n), ScatKdeg, tauFWP
      REAL*8    r,z,   psi, Br, Bz, RBphi, omc, Tee, pe2, pi2, aio, ael
      REAL*8    Bph, kro, kzo, kpho, thi, tho
      r   = yok(1)
      z   = yok(2)
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      Bph = RBphi/r
      tauFWP =  (ScatKdeg*3.1416_R8/180.0_R8)**2
      call GetKscat( yok(4), yok(5), yok(6)/r, Br, Bz, Bph, tauFWP,      &  
     &                  kro, kzo,      kpho, thi, tho )
      yok(4) = kro
      yok(5) = Kzo
      yok(6) = kpho*r
      do 10 i=1,n
        y(i) = yok(i)
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
