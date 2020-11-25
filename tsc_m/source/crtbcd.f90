      subroutine crtbcd(string)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      character*8 string
      call gtext(string,8,0)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
