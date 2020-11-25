!#include "f77_dcomplx.h"
      subroutine  sumwfp (n1, n2, frdum, fzdum)
!*************************************************************************
!
!        sum wire forces due to poloidal currents
!                 & toroidal field - units: kNt/rad
!                 ROS   19 June 1992
      USE CLINAM
      USE SCR21
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n2,n1,ncall,i,j,m,iabs,ii,jj
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 frdum,fzdum,twopi,currdum,ftscpolr,ftscpolz,delx,delz
      REAL*8 g5,g6,g8,g9,fvv
      REAL*8 AREAL
!============
        data  twopi / 6.2831853_R8/
        data  ncall / 0/
!
      frdum = 0.0_R8
      FZdum = 0.0_R8
        if (mmax.le.0)   return
      if (n1.le.0 .or. n1.gt.nwire)   return
      if (n2.lt.n1.or. n2.gt.nwire)   return
!
        if (gzero.lt.0.000001_R8)   return
        ncall = ncall + 1
        if (ncall.eq.1) call pathfind
!
!          Forces (Nt/rad) for TSC poloidal currents
!                     The current and forces for wire i include contributions
!                                     from all paths that touch wire i.
!                     Paths are assumned to be either horizontal or vertical.
!
!      Sign convention that of pubf8.03   23 mar 90.   ROS.
!
!           delta Phi > 0   --->   cvvpol < 0   --->   compressive forces !!
!
        frdum = 0.0_R8
        fzdum = 0.0_R8
      currdum = 0.0_R8
        do 80  i=n1,n2
        polwir(i) = 0.0_R8
        ftscpolr = 0.0_R8
        ftscpolz = 0.0_R8
        do 70  j=1,8
        if (ipath(i,j) .eq.0)   go to 75
!                                               Have a hit
        m = iabs(ipath(i,j))
        polwir(i) = polwir(i) + polsegc(m) * AREAL(m/ipath(i,j))
        delx = xpc2(m) - xpc1(m)
        delz = zpc2(m) - zpc1(m)
        ii = iwire(i)
        jj = jwire(i)
        g5 = g(ii,jj)*xsqoj(ii)
        g6 = g(ii+1,jj)*xsqoj(ii+1)
        g8 = g(ii,jj+1)*xsqoj(ii)
        g9 = g(ii+1,jj+1)*xsqoj(ii+1)
        fvv   = polsegc(m) * 0.25_R8* (g5+g6 + g8+g9) / (xwire(i) *      &  
     & twopi)
        if (delx.eq.0.0_R8)   ftscpolr = ftscpolr - delz * fvv * 0.5_R8
        if (delz.eq.0.0_R8)   ftscpolz = ftscpolz + delx * fvv * 0.5_R8
   70   continue
   75   frdum = frdum + ftscpolr
        fzdum = fzdum + ftscpolz
        currdum = currdum + polwir(i)
   80   continue
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
