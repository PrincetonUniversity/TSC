 
      SUBROUTINE mkdvp(n, delv, v)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8                                                             &  
     &     delv(n), v(n)
      do 10 i = 1, n - 1
         delv(i) = v(i + 1) - v(i)
 10   continue
      delv(n) = delv(n - 1)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
