 
      SUBROUTINE mkvinv(n, vto, vfrom)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8                                                             &  
     &     vto(n), vfrom(n)
      do 10 i = 1, n
         vto(i) = 1._R8/ vfrom(i)
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
