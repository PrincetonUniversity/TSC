 
 
!     FUNCTION bsi0, bsi1                   ---------------------------|
!                                                                      |
!                                                                      |
      REAL*8 FUNCTION bsi0(x)
!
!     BsI0(x) = e^{-x} * I_0 (x)
!
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 TINY, HUGE, x, xx
      PARAMETER ( TINY = 1.0E-05_R8, HUGE = 20._R8)
      if ( x .lt. HUGE ) then
        bsi0 = 1._R8
        if ( x .lt. TINY ) return
        xx = x/2._R8
        xx = xx*xx
        bsi0 = exp(-x)*(1._R8+xx+xx**2/4._R8+xx**3/36._R8)
        return
 
      else
        bsi0 = 0.62E-04_R8
        call LSCbigK( ' lambda too large in function bsi0')
      endif
 
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
