!
!
      SUBROUTINE BrodCast(n, vec, scalar)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8 vec(n), scalar
      do 10 i = 1, n
        vec(i) = scalar
 10   continue
      return
      END
!
!     matrix manipulation routines end here
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
