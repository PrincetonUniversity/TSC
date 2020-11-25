!     -----------------------------------------------------------------
      SUBROUTINE FeMkNorm
!     determine normalized value of fe at ivZero vs psi
      USE FeBins
      USE params
      USE ProfBody
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip
      REAL*8    TailFact
      TailFact = (1._R8+ TailNeps*sqrt(TailTeps))/(1._R8+ TailNeps)
      do 10 ip = 1, npsi
         FeNorm(ip) = fe0 * NeAry(ip) / Vtherm(ip) * TailFact
10    continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
