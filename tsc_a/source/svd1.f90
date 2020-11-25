!#include "f77_dcomplx.h"
      subroutine svd1( m, n, matu, matv,                                 &  
     & ierr, relerr, tau, wmax, wmin, wrat,                              &  
     & rsq, rmsnorm, null, ifstsvd)
!......3.42 svd1
      USE PARAM
      USE SVDCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!----------------------------------------------------------------------
!
!     this routine does least squares fitting to data using the
!     singular value decomposition method. it is substantially
!     a copy of the routine given in the book:
!     computer methods for mathematical computations
!     g.e. forsythe et al
!     prentice hall series in automatic computation (1977)
!
!     the problem fit by least squares is of the form:
!     sum j = 1 to n of coef(j)*f(j,x(i)) = datav(i)
!     where x(i) is the i-th independent coordinate point
!     and datav(i) is the i-th dependent coordinate point
!     the functions f(j,x) are those whose best fit coefficients
!     are sought while the coef(j) are those coefficients.
!     am(i,j) is the design matrix and = f(j,x(i))
!
!     this routine was modified from the original by m. f. reusch
!                                             and by s. c. Jardin
!
!----------------------------------------------------------------------
!
!     this routine determines the singular value decomposition
!                  tr
!     a = u * s * v      of a real m by n rectangular matrix.
!     householder bidiagonalization and a variant of the qr
!     algorithm are used.
!
!     on input:
!
!     nm - must be set to the row dimension of the two dimensional
!     arrays as declared in the calling program dimension statement
!     note that nm must be at least as large as the maximum of
!     m and n.
!
!     m - is the number of rows of a and u
!
!     n - is the number of columns of a and u and the order of v
!
!     a - contains the rectangular input matrix to be decomposed
!
!     matu - should be set to 1 if the u matrix in the decomposition
!     is desired, and to 0 otherwise.
!
!     matv - should be set to 1 if the v matrix in the decomposition
!     is desired, and to 0 otherwise.
!
!     data - are the m values of the dependent variable at the
!     n independent variable coordinate points and is used
!     only if both matu and matv are 1.
!
!     relerr - is the relative error tolerance of the data.
!
!     on output:
!
!     a - is unaltered (unless overwritten by u or v (not done))
!
!     w - contains the n (non-negative) singular values of a
!     (the diagonal elements of s). they are unordered.
!     if an error exit is made, the singular values should be
!     correct for indices ierr + 1, ierr + 2, ... , n.
!
!     u - contains the matrix u (orthogonal column vectors) of the
!     decomposition if matu has been set to 1 otherwise,
!     u is used as a temporary array. u may coincide with a.
!     if an error exit is made, the columns of u corresponding
!     to indicies of correct singular values shoul be correct.
!
!     v - contains the matrix v (orthogonal) of the decomposition if
!     matv has been set to 1 otherwise v is not referenced.
!     v may also coincide with a if u is not needed. if an error exit
!     is made, the columns of v corresponding to indicies of
!     correct singular values should be correct.
!
!     ierr is set to:
!                     zero for a normal return
!                     k    if the k-th singular value has not been
!                          determined after 30 iterations.
!
!     the following are generated only if matu and matv are 1
!
!     wmax and wmin - are the maximum and minimum (.ne.0) singular
!     values respectively, while wrat = wmax/wmin
!
!     rsq - is the rms normed residual error.
!
!     rmsnorm - is the rms norm of the data.
!
!     tau - is the glb to the singular values used to get the coef.
!
!     coef - are the coefficients.
!
!     null - is a count of those singular values less than tau
!
!
!
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,matu,matv,ierr,null,ifstsvd,m,i,j,l,k,ii,mn,kk,k1
      INTEGER its,ll,l1,i1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 relerr,tau,wmax,wmin,wrat,rsq,rmsnorm,g,scale,anorm,s
      REAL*8 f,h,c,y,z,x
      REAL*8 AREAL
!============
      ierr = 0
      if(ifstsvd.eq.0) go to 2000
!
      do 100 i = 1,m
      do 100 j = 1,n
      uu(i,j) = am(i,j)
  100 continue
!
      g = 0.0_R8
      scale = 0.0_R8
      anorm = 0.0_R8
!
      do 300 i = 1,n
      l = i + 1
      rv1(i) = scale * g
      g = 0.0_R8
      scale = 0.0_R8
      s = 0.0_R8
      if(i.gt.m)go to 210
!
      do 120 k = i,m
  120 scale = scale + abs(uu(k,i))
!
      if(scale.eq.0.0_R8)go to 210
!
      do 130 k = i,m
      uu(k,i) = uu(k,i) / scale
      s = s + uu(k,i) * uu(k,i)
  130 continue
!
      f = uu(i,i)
      g =  - sign(sqrt(s),f)
      h = f * g - s
      uu(i,i) = f - g
      if(i.eq.n)go to 190
!
      do 150 j = l,n
      s = 0.0_R8
!
      do 140 k = i,m
  140 s = s + uu(k,i) * uu(k,j)
!
      f = s / h
!
      do 150 k = i,m
      uu(k,j) = uu(k,j) + f * uu(k,i)
  150 continue
!
  190 do 200 k = i,m
  200 uu(k,i) = scale * uu(k,i)
!
  210 sigma(i) = scale * g
      g = 0.0_R8
      s = 0.0_R8
      scale = 0.0_R8
      if(i.gt.m.or.i.eq.n)go to 290
!
      do 220 k = l,n
  220 scale = scale + abs(uu(i,k))
!
      if(scale.eq.0.0_R8)go to 290
!
      do 230 k = l,n
      uu(i,k) = uu(i,k) / scale
      s = s + uu(i,k) * uu(i,k)
  230 continue
!
      f = uu(i,l)
      g =  - sign(sqrt(s),f)
      h = f * g - s
      uu(i,l) = f - g
!
      do 240 k = l,n
  240 rv1(k) = uu(i,k) / h
!
      if(i.eq.m)go to 270
!
      do 260 j = l,m
      s = 0.0_R8
!
      do 250 k = l,n
  250 s = s + uu(j,k) * uu(i,k)
!
      do 260 k = l,n
      uu(j,k) = uu(j,k) + s * rv1(k)
  260 continue
!
  270 do 280 k = l,n
  280 uu(i,k) = scale * uu(i,k)
!
  290 anorm = max(anorm,abs(sigma(i)) + abs(rv1(i)))
  300 continue
!
!     accumulation of right-hand transformations
!
      if(matv.eq.0)go to 410
!
!     for i = n step -1 until 1 do
!
      do 400 ii = 1,n
      i = n + 1 - ii
      if(i.eq.n)go to 390
      if(g.eq.0.0_R8)go to 360
!
      do 320 j = l,n
!
!     double division avoids possible underflow
!
  320 vv(j,i) = (uu(i,j) / uu(i,l)) / g
!
      do 350 j = l,n
      s = 0.0_R8
!
      do 340 k = l,n
  340 s = s + uu(i,k) * vv(k,j)
!
      do 350 k = l,n
      vv(k,j) = vv(k,j) + s * vv(k,i)
  350 continue
!
  360 do 380 j = l,n
      vv(i,j) = 0.0_R8
      vv(j,i) = 0.0_R8
  380 continue
!
  390 vv(i,i) = 1.0_R8
      g = rv1(i)
      l = i
  400 continue
!
!     accumulation of left-hand transformations
!
  410 if(matu.eq.0)go to 510
!
!     for i = min(m,n) step -1 until 1 do
!
      mn = n
      if(m.lt.n)mn = m
!
      do 500 ii = 1,mn
      i = mn + 1 - ii
      l = i + 1
      g = sigma(i)
      if(i.eq.n)go to 430
!
      do 420 j = l,n
  420 uu(i,j) = 0.0_R8
!
  430 if(g.eq.0.0_R8)go to 475
      if(i.eq.mn)go to 460
!
      do 450 j = l,n
      s = 0.0_R8
!
      do 440 k = l,m
  440 s = s + uu(k,i) * uu(k,j)
!
!     double division avoids possible underflow
!
      f = (s / uu(i,i)) / g
!
      do 450 k = i,m
      uu(k,j) = uu(k,j) + f * uu(k,i)
  450 continue
!
  460 do 470 j = i,m
  470 uu(j,i) = uu(j,i) / g
!
      go to 490
!
  475 do 480 j = i,m
  480 uu(j,i) = 0.0_R8
!
  490 uu(i,i) = uu(i,i) + 1.0_R8
  500 continue
!
!     diagonalization of the bidiagonal form
!
!     for k = n step -1 until 1 do
!
  510 do 700 kk = 1,n
      k1 = n - kk
      k = k1 + 1
      its = 0
!
!     test for splitting
!
!     for l = k step -1 until 1 do
!
  520 do 530 ll = 1,k
      l1 = k - ll
      l = l1 + 1
      if(abs(rv1(l)) + anorm.eq.anorm)go to 565
!
!     rv1(1) is always zero, so there is no exit
!     through the bottom of the loop
!
      if(abs(sigma(l1)) + anorm.eq.anorm)go to 540
!
  530 continue
!
!     cancellation of rv1(l) if l greater than 1
!
  540 c = 0.0_R8
      s = 1.0_R8
!
      do 560 i = l,k
      f = s * rv1(i)
      rv1(i) = c * rv1(i)
      if(abs(f) + anorm.eq.anorm)go to 565
      g = sigma(i)
      h = sqrt(f * f + g * g)
      sigma(i) = h
      c = g / h
      s =  - f / h
      if(matu.eq.0)go to 560
!
      do 550 j = 1,m
      y = uu(j,l1)
      z = uu(j,i)
      uu(j,l1) = y * c + z * s
      uu(j,i) =  - y * s + z * c
  550 continue
!
  560 continue
!
!     test for convergence
!
  565 z = sigma(k)
      if(l.eq.k)go to 650
!
!     shift from bottom 2 by 2 minor
!
      if(its.eq.30)go to 1000
!
      its = its + 1
      x = sigma(l)
      y = sigma(k1)
      g = rv1(k1)
      h = rv1(k)
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0_R8* h * y)
      g = sqrt(f * f + 1.0_R8)
      f = ((x - z) * (x + z) + h * (y / (f + sign(g,f)) - h)) / x
!
!     next qr transformation
!
      c = 1.0_R8
      s = 1.0_R8
!
      do 600 i1 = l,k1
      i = i1 + 1
      g = rv1(i)
      y = sigma(i)
      h = s * g
      g = c * g
      z = sqrt(f * f + h * h)
      rv1(i1) = z
      c = f / z
      s = h / z
      f = x * c + g * s
      g =  - x * s + g * c
      h = y * s
      y = y * c
!
      if(matv.eq.0)go to 575
!
      do 570 j = 1,n
      x = vv(j,i1)
      z = vv(j,i)
      vv(j,i1) = x * c + z * s
      vv(j,i) =  - x * s + z * c
  570 continue
!
  575 z = sqrt(f * f + h * h)
      sigma(i1) = z
!
!     rotation can be arbitrary if z is zero
!
      if(z.eq.0.0_R8)go to 580
      c = f / z
      s = h / z
  580 f = c * g + s * y
      x =  - s * g + c * y
!
      if(matu.eq.0)go to 600
!
      do 590 j = 1,m
      y = uu(j,i1)
      z = uu(j,i)
      uu(j,i1) = y * c + z * s
      uu(j,i) =  - y * s + z * c
  590 continue
!
  600 continue
!
      rv1(l) = 0.0_R8
      rv1(k) = f
      sigma(k) = x
      go to 520
!
!     convergence
!
  650 if(z.ge.0.0_R8)go to 700
!
!     sigma(k) is made non-negative
!
      sigma(k) =  - z
!
      if(matv.eq.0)go to 700
!
      do 690 j = 1,n
  690 vv(j,k) =  - vv(j,k)
!
  700 continue
!
      go to 2000
!
!     set error if no convergence to a singular value
!     after 30 iterations
!
 1000 ierr = k
      return
!
 2000 continue
!
!     find maximum and minimum singular values
!
      wmax = -1.0E+30_R8
      wmin =  1.0E+30_R8
!
      do 710 i = 1,n
      wmax=max(wmax,sigma(i))
      if(sigma(i).gt.0.0_R8)wmin=min(wmin,sigma(i))
      coef(i)=0.0_R8
  710 continue
!
      tau = relerr * wmax
!
      if(matu*matv.eq.0)return
!
!     transform data with u and then use v
!     to get the coefficients (coef)
!
      null=0
!
      do 740 j = 1,n
      if(sigma(j).le.tau)go to 735
      s = 0.0_R8
!
      do 720 i = 1,m
      s = s + uu(i,j) * datav(i)
  720 continue
!
      s = s / sigma(j)
!
      do 730 i = 1,n
      coef(i) = coef(i) + s * vv(i,j)
  730 continue
!
      go to 740
  735 null = null + 1
  740 continue
!
!     form square root of sum of residuals and norm
!
      rsq = 0.0_R8
      rmsnorm = 0.0_R8
!
      do 760 i = 1,m
      s = 0.0_R8
!
      do 750 j = 1,n
      s = s + coef(j) * am(i,j)
  750 continue
!
      rsq = rsq + (s - datav(i))**2
      rmsnorm = rmsnorm + datav(i)*datav(i)
  760 continue
!
      rsq = sqrt(rsq / AREAL(m))
      rmsnorm = sqrt(rmsnorm / AREAL(m))
      rsq = rsq / rmsnorm
      wrat = wmax/wmin
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
