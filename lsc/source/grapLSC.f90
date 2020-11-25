!
!----------------------------------------------------------------------
!
      SUBROUTINE grapLSC(imode, xx, zz, psval, psderiv, isw,             &  
     &     xmin, deex, zmin, deez, isave, jsave, nx, nz, nxdim,          &  
     &     psi, fmat)
!***************************************************************
!     *
!     GRid APproximation in the Lower hybrid Simulation Code
!     local bivariate interpolation*
!     *
!***************************************************************
!
!     imode=0   returns everything
!     =1   returns psval only
!.....
!
!     isw=1 forces recalculation of coefficient matrix
!     =0 may not recalculate coefficient matrix if same
!     i,j as last call
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER imode, isw, jval, ival, m, k, l, ii, i, j, nx, nz,         &  
     &     iplus, iminus, jplus, jminus, isave, jsave, frstcall, nxdim
      INTEGER LLC, LRC, ULC, URC, LEFT,                                  &  
     &     RIGHT, TOP, BOTTOM, INTERIOR, FALSE, TRUE
      PARAMETER(                                                         &  
     &     FALSE = 0, TRUE = 1,  INTERIOR = 1,                           &  
     &     LLC = 5,  LRC = 6,    ULC = 7,   URC = 8,                     &  
     &     LEFT = 9, RIGHT = 10, TOP = 11,  BOTTOM = 12)
      REAL*8                                                             &  
     &     zz, xx, gradsq, dpsidx, dpsidz, psval, psixz,                 &  
     &     psixx, psizz, psi, deez, deex, zmin,                          &  
     &     xmin, rival, rjval, cx, cy, psderiv
      REAL*8                                                             &  
     &     vmat, fmat, poly, xpa, zpa, rhs1, rhs2, rhs3, rhs4, sum1,     &  
     &     sum2, sum3, sum4, sum5, sum6, xp, zp, ff, f1
      DATA frstcall / TRUE /
      DIMENSION vmat(4,4),fmat(0:3,0:*),poly(0:3,0:3)
      DIMENSION xpa(0:3),zpa(0:3)
      DIMENSION psi(nxdim, nz), psderiv(0:2, 0:*)
      DATA xpa(0),zpa(0)/ 1._R8, 1._R8/
!
      rjval = (zz-zmin)/deez+1._R8
      rival = (xx-xmin )/deex+1._R8
      jval = rjval
      ival = rival
!
      if(ival .lt. 2)then
         ival = 2
      else if(ival .gt. nx - 2)then
         ival = nx - 2
      endif
      if(jval .lt. 2)then
         jval = 2
      else if(jval .gt. nz - 2)then
         jval = nz - 2
      endif
      if( (frstcall .eq. FALSE) .and.                                    &  
     &     (ival.eq.isave .and. jval.eq.jsave .and. isw.eq.0) )then
         go to 200
      endif
      frstcall = FALSE
      isave = ival
      jsave = jval
!
!     begin INTERIOR point
      do 10 ii=1,4
         go to(101,102,103,104),ii
 101     i      = ival
         j      = jval
         iplus = i + 1
         jplus = j + 1
         iminus = i - 1
         jminus = j - 1
         cx = 0.50_R8
         cy = 0.50_R8
         go to 105
 102     i      = ival
         j      = jval+1
         iplus = i + 1
         jplus = j + 1
         iminus = i - 1
         jminus = j - 1
         cx = 0.50_R8
         cy = 0.50_R8
         go to 105
 103     i      = ival+1
         j      = jval
         iplus = i + 1
         jplus = j + 1
         iminus = i - 1
         jminus = j - 1
         cx = 0.50_R8
         cy = 0.50_R8
         go to 105
 104     i      = ival+1
         j      = jval+1
         iplus = i + 1
         jplus = j + 1
         iminus = i - 1
         jminus = j - 1
         cx = 0.50_R8
         cy = 0.50_R8
 105     continue
!
!...psi:
!
         vmat(1, ii) =      psi(i, j)     - psi(ival, jval)
         vmat(2, ii) = cx * (psi(iplus, j  ) - psi(iminus, j) )
         vmat(3, ii) = cy * (psi(i, jplus) - psi(i, jminus)  )
         vmat(4, ii) = cx * cy * (psi(iplus, jplus) -                    &  
     &        psi(iminus, jplus) - psi(iplus, jminus) +                  &  
     &        psi(iminus, jminus) )
!
 10   continue
!     end INTERIOR point
!
!
!..........................................................
!
!     calculate function values and function derivatives
!     at four corners of the reference cell
!
!..........................................................
!..........................................................c
!     determine coefficients by evaluating polynomial
!     at four corner points c
!..........................................................c
!...point(i,j):c
      fmat(0,0) = 0._R8
      fmat(1,0) = vmat(2,1)
      fmat(0,1) = vmat(3,1)
      fmat(1,1) = vmat(4,1)
!
!...point(i,j+1):
!
      fmat(0,2) = 3._R8*vmat(1,2) - 2._R8*fmat(0,1)                      &  
     &     -    vmat(3,2)
      fmat(0,3) =    vmat(3,2) +    fmat(0,1)                            &  
     &     - 2._R8*vmat(1,2)
      fmat(1,2) = 3._R8*vmat(2,2) - 3._R8*fmat(1,0)                      &  
     &     -    vmat(4,2) - 2._R8*fmat(1,1)
      fmat(1,3) =    vmat(4,2) + 2._R8*fmat(1,0)                         &  
     &     - 2._R8*vmat(2,2) +    fmat(1,1)
!
!...point(i+1,j):
!
      fmat(2,0) = 3._R8*vmat(1,3) - 2._R8*fmat(1,0)                      &  
     &     -    vmat(2,3)
      fmat(3,0) =    vmat(2,3) +    fmat(1,0)                            &  
     &     - 2._R8*vmat(1,3)
      fmat(2,1) = 3._R8*vmat(3,3) - 3._R8*fmat(0,1)                      &  
     &     -    vmat(4,3) - 2._R8*fmat(1,1)
      fmat(3,1) =    vmat(4,3) + 2._R8*fmat(0,1)                         &  
     &     - 2._R8*vmat(3,3) +    fmat(1,1)
!
!...point(i+1,j+1):
!
      rhs1=   vmat(1,4)-vmat(1,3)-vmat(3,3)                              &  
     &     -  (fmat(0,2)+fmat(1,2)                                       &  
     &     +   fmat(0,3)+fmat(1,3) )
      rhs2=   vmat(2,4)-vmat(2,3)-vmat(4,3)                              &  
     &     -  (fmat(1,2)+fmat(1,3) )
      rhs3=   vmat(3,4)-vmat(3,3)                                        &  
     &     -2*(fmat(0,2)+fmat(1,2) )                                     &  
     &     -3*(fmat(0,3)+fmat(1,3) )
      rhs4=   vmat(4,4)-vmat(4,3)                                        &  
     &     -2*(fmat(1,2) )                                               &  
     &     -3*(fmat(1,3) )
!
      fmat(2,2) =  9*rhs1-3*rhs2-3*rhs3+rhs4
      fmat(3,2) = -6*rhs1+3*rhs2+2*rhs3-rhs4
      fmat(2,3) = -6*rhs1+2*rhs2+3*rhs3-rhs4
      fmat(3,3) =  4*rhs1-2*rhs2-2*rhs3+rhs4
 200  continue
!...........................................................
!
!     evaluate function and derivatives at (xx,zz)
!
!...........................................................
!      if(rival .le. 1.)then
!         rival = 1.
!      else if(rival .gt. float(nx))then
!         rival = float(nx)
!      endif
!      if(rjval .le. 1.)then
!         rjval = 1.
!      else if(rjval .gt. float(nz))then
!         rjval = float(nz)
!      endif
      xp     = rival - ival
      zp     = rjval - jval
!
!     begin evaluate function
      do 32 m=1,3
         xpa(m) = xp*xpa(m-1)
 32      zpa(m) = zp*zpa(m-1)
         do 33 l=0,3
            do 33 k=0,3
 33            poly(k,l) = xpa(k)*zpa(l)
!
!
               sum1   = 0._R8
               do 40 k = 0,3
                  do 45 l = 0,3
                     f1     = fmat(k,l)
                     sum1   = sum1 + f1*poly(k,l)
 45               continue
 40            continue
               psval  = sum1 + psi(ival,jval)
!     end evaluate function
!
               if(imode.eq.1) then
                  do 47 k = 0, 2
                     do 48 l = 0, 2
                        psderiv(k, l) = 0._R8
 48                  continue
 47               continue
!
                  return
!
               endif
!
!     begin evaluate derivatives
               sum2   = 0._R8
               do 55 l = 0,3
                  do 50 k = 0,2
                     ff     = (k+1)*fmat(k+1,l)
                     sum2   = sum2 + ff*poly(k,l)
 50               continue
 55            continue
               dpsidx = sum2/deex
!
               sum3   = 0._R8
               do 75 k = 0,3
                  do 70 l = 0,2
                     ff     = (l+1)*fmat(k,l+1)
                     sum3   = sum3 + ff*poly(k,l)
 70               continue
 75            continue
               dpsidz = sum3/deez
               gradsq = dpsidx**2 + dpsidz**2
!
               sum5   = 0._R8
               do 65 l = 0,3
                  do 60 k = 0,1
                     ff     = (k+1)*(k+2)*fmat(k+2,l)
                     sum5   = sum5 + ff*poly(k,l)
 60               continue
 65            continue
               psixx  = sum5/deex/deex
!
!
!
               sum6   = 0._R8
               do 85 k = 0,3
                  do 80 l = 0,1
                     ff     = (l+1)*(l+2)*fmat(k,l+2)
                     sum6   = sum6 + ff*poly(k,l)
 80               continue
 85            continue
               psizz  = sum6/deez/deez
!
!
               sum4   = 0._R8
               do 90 k = 0,2
                  do 95 l = 0,2
                     ff     = (k+1)*(l+1)*fmat(k+1,l+1)
                     sum4   = sum4 + ff*poly(k,l)
 95               continue
 90            continue
               psixz  = sum4/deex/deez
!     end evaluate derivatives
!
!
!
               psderiv(0, 0) = gradsq
               psderiv(1, 0) = dpsidx
               psderiv(2, 0) = psixx
               psderiv(0, 1) = dpsidz
               psderiv(0, 2) = psizz
               psderiv(1, 1) = psixz
               return
               END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
