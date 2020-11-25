!
!----------------------------------------------------------------------
!
      SUBROUTINE GNzero(x,y,n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, n
      REAL*8    x(n), y(n)
      REAL*8    SMALL
 
      SMALL = 1.0E-30_R8
 
      do i = 1, n
         if ( abs (y(i)) .lt. SMALL) y(i) = 0.0_R8
      enddo
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
