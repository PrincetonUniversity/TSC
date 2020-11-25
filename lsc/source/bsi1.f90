!                                                                      |
!                                                                      |
!                                                                      |
      REAL*8 FUNCTION bsi1(x)
!
!     BsI1(x) = e^{-x} * I_1 (x)
!
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 TINY, HUGE, x, xx
      PARAMETER ( TINY = 1.0E-05_R8, HUGE = 20._R8)
 
      if ( x .lt. HUGE ) then
        xx = x/2._R8
        bsi1 = exp(-x)*abs(xx)
        if ( x .lt. TINY ) then
          return
        else
          xx = xx*xx
          bsi1=bsi1*(1._R8+xx/2._R8+xx**2/12._R8+xx**3/144._R8)
          return
        endif
 
      else
        bsi1 = 0.16E-03_R8
        call LSCbigK( ' lambda too large in function bsi1')
      endif
      END
!                                                                      |
!                                                                      |
!     bsi0, bsi1 end                        ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
