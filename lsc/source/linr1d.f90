!
!     ------------------------------------------------------------------
!
      SUBROUTINE linr1d(nx, xdata, ydata, x, y )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER j, jguess, nx
      REAL*8    xdata(nx), ydata(nx), x, y, coef(2)
      REAL*8    d21, dx
!     LINeaR interpolation in 1D
!
!     Try to look much like grnu1d, which see above.
!     The derivative is guaranteed not continuous.
!     Given to make consistent with grnu1d.
!     Linear interpolation of interest to avoid fake oscillations
!     in the interpolant.
!
!     here: yint = sum(1 - 2) * coef(j) * ((x - x(1)) ** (j - 1))
!
!                                       find which interval we are in by
!                                       Numerical Recipes HUNT routine
      jguess = -1
      call huntnr(xdata, nx, x, jguess)
      j = jguess
 
      if (j .eq. 0 ) then
        y  = ydata(1)
!       yp = 0.                         yp not of interest since its not
!                                       continuous.
        return
      endif
      if ( j .eq. nx) then
        y  = ydata(nx)
!       yp = 0.                         yp not of interest
        return
      endif
!
!     This is the computation of the slope, and delta-x
            d21     = (ydata(j + 1) - ydata(j    )) /                    &  
     &                (xdata(j + 1) - xdata(j    ))
 
            coef(1) = ydata(j)
            coef(2) = d21
            dx = x - xdata(j)
!
      y  = coef(2) * dx + coef(1)
!     yp = coef(2)                      yp not of interest
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
