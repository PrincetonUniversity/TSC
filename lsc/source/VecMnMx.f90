!
!     -----------------------------------------------------------------
      SUBROUTINE VecMnMx(vec, n, vmin, vmax)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n
      REAL*8    vec(n)
      REAL*8    vmin, vmax
      INTEGER i
      i = 1
      vmin = vec(i)
      vmax = vec(i)
      do 20 i = 2, n
              if(vmin .gt. vec(i))then
                      vmin = vec(i)
              else if(vmax .lt. vec(i))then
                      vmax = vec(i)
              endif
 20   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
