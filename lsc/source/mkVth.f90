!     -----------------------------------------------------------------
      SUBROUTINE mkVth
      USE FeBins
      USE params
      USE ProfBody
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip
      do 10 ip = 1, npsi
         Vtherm(ip) = VthNorm * sqrt(TeAry(ip))
         VperpSq(ip) = Vtherm(ip) * Vtherm(ip)
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
