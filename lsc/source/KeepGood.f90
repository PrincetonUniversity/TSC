!
!     -----------------------------------------------------------------
!
      SUBROUTINE KeepGood(y, yok, n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8    y(n), yok(n)
      do 10 i=1,n
        yok(i) = y(i)
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
