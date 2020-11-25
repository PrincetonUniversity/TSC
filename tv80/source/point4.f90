      subroutine point4(x4, y4)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*4 x4, y4
      REAL*8 x8, y8
      x8 = x4
      y8 = y4
      call point(x8, y8)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
