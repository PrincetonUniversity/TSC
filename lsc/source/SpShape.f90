!
      REAL*8 FUNCTION SpShape(ParValue, center, width)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 ParValue, center, width, relnpar, exparg
      relnpar = (ParValue - center) / width
      exparg = 0.5_R8* relnpar * relnpar
      SpShape = exp( - exparg)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
