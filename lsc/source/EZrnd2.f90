!
!----------------------------------------------------------------------
!
      SUBROUTINE EZrnd2( r , x, i, ii )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, ii, iexp
      REAL*8    r, x, a, b
      REAL*8    SmallNum
      DATA    SmallNum / 0.00001_R8/
      REAL*8    RE41
!     r is the given maximum value to be graphed
!     x is the nice graph maximum for plotting
!     i  is the # of major divisions
!     ii is the # of minor divisions
!     Chosen maximum is one of 1 1.5 2 2.5 3 4 5 7.5 10 15 etc
 
      x = abs ( r )
      i = 1
      ii= 2
      if ( r .eq. 0._R8) return
 
      a    = log10 ( x )
      RE41 = a
      iexp = int (RE41)
!     iexp = ifix ( a )
      iexp = 1 - iexp
      if ( a .lt. 0._R8) iexp = iexp + 1
      b = 10.0_R8**  iexp
      x = x*b
 
      if      ( x .le. 10._R8+ SmallNum ) then
        x = 10._R8
        i = 2
        ii= 5
 
      else if ( x .le. 15._R8+ SmallNum ) then
        x = 15._R8
        i = 3
        ii= 5
 
      else if ( x .le. 20._R8+ SmallNum ) then
        x = 20._R8
        i = 2
        ii= 5
 
      else if ( x .le. 25._R8+ SmallNum ) then
        x = 25._R8
        i = 5
        ii= 1
 
      else if ( x .le. 30._R8+ SmallNum ) then
        x = 30._R8
        i = 3
        ii= 2
!!!     i = 2 & ii = 3 bring out a bug in sg such that the scale is screwed up
 
      else if ( x .le. 40._R8+ SmallNum ) then
        x = 40._R8
        i = 2
        ii= 4
 
      else if ( x .le. 50._R8+ SmallNum ) then
        x = 50._R8
        i = 2
        ii= 5
 
      else if ( x .le. 75._R8+ SmallNum ) then
        x = 75._R8
        i = 3
        ii= 5
 
      else
        x = 100._R8
        i  = 2
        ii = 5
 
      endif
 
      x = x / b
      x = sign (x,r)
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
