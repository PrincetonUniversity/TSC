!
!     ------------------------------------------------------------------
!
      SUBROUTINE grnu1d(nx, xdata, ydata, jold, coef, x, y, yp )
      USE params
      USE PlPr
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER j, jold, jguess, nx
      INTEGER FORCE_CALCULATION
      PARAMETER(FORCE_CALCULATION = -1)
      REAL*8    xdata(nx), ydata(nx), x, y, yp, coef(4)
      REAL*8    x21,      x32
      REAL*8    d21, d31, d32, d42, dx
!     GRid, NonUniform, interpolation in 1D
!     cubic interpolation yint(x) for points 1, 2, 3, 4 set by conditions
!     yint(2) = ydata(2)
!     yint(3) = ydata(3)
!     d yint(2) / dx = (ydata(3) - ydata(1)) / (xdata(3) - xdata(1) )
!     d yint(3) / dx = (ydata(4) - ydata(2)) / (xdata(4) - xdata(2) )
!
!     here: yint = sum(1 - 4) * coef(j) * ((x - x(2)) ** (j - 1))
!
!                                       find which interval we are in by
!                                       Numerical Recipes HUNT routine
      jguess = jold
      call huntnr(xdata, nx, x, jguess)
      j = jguess
 
      if (j .eq. 0 ) then
        y  = ydata(1)
        yp = 0._R8
        return
      endif
      if ( j .eq. nx) then
        y  = ydata(nx)
        yp = 0._R8
        return
      endif
 
      if((jold .eq. FORCE_CALCULATION) .or. (j .ne. jold))then
!     compute new coefficients
        jold = j
        if((j .gt. 1) .and. (j .lt. (nx - 1))) then
            x32     = (xdata(j + 1) - xdata(j    ))
            d31     = (ydata(j + 1) - ydata(j - 1)) /                    &  
     &                (xdata(j + 1) - xdata(j - 1))
            d32     = (ydata(j + 1) - ydata(j    )) /                    &  
     &                (xdata(j + 1) - xdata(j    ))
            d42     = (ydata(j + 2) - ydata(j    )) /                    &  
     &                (xdata(j + 2) - xdata(j    ))
 
            coef(1) = ydata(j)
            coef(2) = d31
            coef(3) = (3._R8*(d32 - d31) -    (d42 - d31) ) /  x32
            coef(4) = (   (d42 - d31) - 2._R8*(d32 - d31) ) / (x32*x32)
         else if(j .eq. 1)then
            x21     = (xdata(2    ) - xdata(    1))
            d21     = (ydata(2    ) - ydata(    1)) /                    &  
     &                (xdata(2    ) - xdata(    1))
            d31     = (ydata(3    ) - ydata(    1)) /                    &  
     &                (xdata(3    ) - xdata(    1))
            coef(1) = ydata(1)
            coef(2) = (2._R8*d21 - d31)
            coef(3) = (d31 - d21) / x21
            coef(4) = 0._R8
 
         else if(j .eq. (nx - 1) )then
            x32     = (xdata(nx   ) - xdata(nx- 1))
            d31     = (ydata(nx   ) - ydata(nx- 2)) /                    &  
     &                (xdata(nx   ) - xdata(nx- 2))
            d32     = (ydata(nx   ) - ydata(nx- 1)) /                    &  
     &                (xdata(nx   ) - xdata(nx- 1))
            coef(1) = ydata(nx - 1)
            coef(2) =  d31
            coef(3) = (d32 - d31) / x32
            coef(4) = 0._R8
         endif
      endif
!     compute y, yp
      dx = x - xdata(j)
      y  = ((coef(4) * dx + coef(3)) * dx + coef(2)) * dx                &  
     &     + coef(1)
      yp = ((3._R8* coef(4) * dx + 2._R8* coef(3)) * dx + coef(2))
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
