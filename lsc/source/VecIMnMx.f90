!     -----------------------------------------------------------------
      SUBROUTINE VecIMnMx(vec, n, vmin, vmax)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,i
!============
      INTEGER vec, vmin, vmax
      DIMENSION vec(n)
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
 
!
!                                                                      |
!                                                                      |
!     GridGen.f(or)     ends                ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
