!     ------------------------------------------------------------------
      REAL*8 FUNCTION MaxOfAry(array, n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8    array(n), tmp
      tmp = array(1)
      do 10 i = 2, n
        if(array(i) .gt. tmp) tmp = array(i)
 10   continue
      MaxOfAry = tmp
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
