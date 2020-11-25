!
!     -----------------------------------------------------------------|
!
      SUBROUTINE EZrndu(r,x,i)
!      x is the nice graph max for plotting
!      i is the #of major divisions
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, iexp
      REAL*8    r, x, a
      REAL*8    SmallNum
      REAL*8    RE41
      DATA    SmallNum / 0.00001_R8/
      x = abs ( r )
      i = 1
      if ( r .eq. 0._R8) return
      a = log10 ( x )
      RE41 = a
      iexp = int (RE41)
!     iexp = ifix ( a )
      iexp = 1 - iexp
      if ( a .lt. 0._R8) iexp = iexp + 1
      x = x*(10._R8**iexp)
      if ( x .gt. 10._R8+ SmallNum ) go to 20
      x = 10._R8
      i = 2
      go to 2000
20    if ( x .gt. 20._R8+ SmallNum ) go to 50
      x = 20._R8
      i = 2
      go to 2000
50    if ( x .gt. 50._R8+ SmallNum ) go to 100
      x = 50._R8
      i = 2
      go to 2000
100   x = 100._R8
      i  = 2
2000  x = x / (10._R8**iexp)
      x = sign (x,r)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
