      subroutine line4(x14, y14, x24, y24)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*4 x14, y14, x24, y24
      REAL*8 x18, y18, x28, y28
      x18 = x14
      x28 = x24
      y18 = y14
      y28 = y24
      call line(x18, y18, x28, y28)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
