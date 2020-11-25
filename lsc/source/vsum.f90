!     ------------------------------------------------------------------
      SUBROUTINE vsum(n, vec, sum)
!                                       Adds n elements of vec, returns sum
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8                                                             &  
     &     vec(n), sum
      sum = 0._R8
      do 10 i = 1, n
         sum = sum + vec(i)
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
