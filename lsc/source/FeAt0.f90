!     -----------------------------------------------------------------
      SUBROUTINE FeAt0
!     value of fe at t = 0.
      USE FeBins
      USE params
      USE ProfBody
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL FeMaxw
      INTEGER iipsi
!                                       initial conditions for fe(v, psi)
      do 10 iipsi = 1, npsi
         call FeMaxw(iipsi)
!                                       set fe at iipsi vs velocity
10    continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
