!
!     ------------------------------------------------------------------
!
      SUBROUTINE grrg1d(nx, xmin, dx, ydata, jold, xx, yy, ypr, coef)
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER FORCE_CALCULATION
      PARAMETER(FORCE_CALCULATION = -1)
      INTEGER nx, jold, j
      REAL*8                                                             &  
     &        xmin, dx, ydata(nx), xx, yy, ypr, coef(*), rj,             &  
     &        left1, left2, xrel, ydpr
!     GRid, ReGulalr, interpolation in 1D
!     cubic interpolation yint(x) for points 1, 2, 3, 4 set by conditions
!     yint(2) = ydata(2)
!     yint(3) = ydata(3)
!     d yint(2) / dx = (ydata(3) - ydata(1)) / (2 * dx)
!     d yint(3) / dx = (ydata(4) - ydata(2)) / (2 * dx)
!
!     here: yint = sum(1 - 4) * coef(j) * ((x - x(2)) ** (j - 1))
      rj = (xx - xmin) / dx + 1._R8
      j = aint(rj)
      if(j .lt. 1) then
        yy = ydata(1)
        ypr= 0._R8
           write(nLSCcom2, '( '' x out of range LOW  in grrg1d: '' )')
        return
      endif
      if(j .gt. (nx - 1))then
        yy = ydata(nx)
        ypr= 0._R8
           write(nLSCcom2, '( '' x out of range HIGH in grrg1d: '' )')
        return
      endif
 
      if((jold .eq. FORCE_CALCULATION) .or. (j .ne. jold))then
!     compute new coefficients
         if((j .gt. 1) .and. (j .lt. (nx - 1)))then
            coef(1) = ydata(j)
            coef(2) = 0.5_R8* (ydata(j + 1) - ydata(j - 1))
            left1 = (ydata(j + 1) - (coef(1) + coef(2)))
            left2 = (0.5_R8* (ydata(j + 2) - ydata(j))- coef(2))
            coef(3) = -left2 + 3._R8* left1
            coef(4) = (left2 - 2._R8* left1)
         else if(j .eq. 1)then
            coef(1) = ydata(1)
            ydpr = 0.5_R8* (ydata(3) - ydata(1))
            coef(3) = ydpr - (ydata(2) - ydata(1))
            coef(2) = ydpr - 2._R8* coef(3)
            coef(4) = 0._R8
         else if(j .eq. (nx - 1) )then
            coef(1) = ydata(nx - 1)
            coef(2) = 0.5_R8* (ydata(nx) - ydata(nx - 2))
            coef(3) = (ydata(nx) - coef(1) - coef(2))
            coef(4) = 0._R8
         endif
      endif
!     compute yy, ypr
      xrel = rj - j
      yy = ((coef(4) * xrel + coef(3)) * xrel + coef(2)) * xrel          &  
     &     + coef(1)
      ypr = ((3._R8* coef(4) * xrel + 2._R8* coef(3)) * xrel +           &  
     &      coef(2)) /                                                   &  
     &     dx
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
