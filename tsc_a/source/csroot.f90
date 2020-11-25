      subroutine csroot(xr,xi,yr,yi)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 xr,xi,yr,yi
 
!     (yr,yi) = complex sqrt(xr,xi)
!     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
 
      REAL*8 s,tr,ti,pythag,dlapy3gf
      tr = xr
      ti = xi
      s = sqrt(0.5_R8*(dlapy3gf(tr,ti) + abs(tr)))
      if (tr .ge. 0.000_R8) yr = s
      if (ti .lt. 0.000_R8) s = -s
      if (tr .le. 0.000_R8) yi = s
      if (tr .lt. 0.000_R8) yr = 0.5_R8*(ti/yi)
      if (tr .gt. 0.000_R8) yi = 0.5_R8*(ti/yr)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
