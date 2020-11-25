!     ------------------------------------------------------------------
      REAL*8 FUNCTION MinOfAry(array, n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8    array(n), tmp
      tmp = array(1)
      do 10 i = 2, n
        if(array(i) .lt. tmp) tmp = array(i)
 10   continue
      MinOfAry = tmp
      return
      END
!                                                                      |
!                                                                      |
!     end Low Of 1 ch; MaxOfAry; MinOfAry -----------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
