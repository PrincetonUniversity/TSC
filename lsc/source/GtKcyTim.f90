!
!     ------------------------------------------------------------------
!
      SUBROUTINE GtKcyTim (kcycgot, timegot, namegot)
      USE params
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*16 namegot
      INTEGER kcycgot
      REAL*8    timegot
      kcycgot = kcycle
      timegot = times
      namegot = equhd(2) // equhd(3)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
