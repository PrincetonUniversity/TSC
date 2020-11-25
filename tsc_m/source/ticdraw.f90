      subroutine ticdraw (xmin,xmax, ymin,ymax, kydo)
!**********************************************************************
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER kydo,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xmax,ymin,ymax,xmin,tscale,delx,xv,tickl
!============
      data  tscale /0.025_R8/
      delx = 0.5_R8
      xv = -5.5_R8
      if (xmax-xmin.gt. 5.0_R8)  then
               delx = 1.0_R8
               xv = -11.0_R8
               endif
      if (xmax-xmin.gt.10.0_R8)  then
               delx = 2.0_R8
               xv = -22.0_R8
               endif
      tickl = tscale * (xmax - xmin)
      if (kydo.eq.-1)  tickl = tickl * 1.0_R8* (ymax-ymin)/(xmax-xmin)
      do 40  k=1,30
      xv = xv + delx
      if (xv.gt.xmax)      go to 50
      call setcrt (xv, ymin)
      call vector (xv, ymin+tickl)
      call setcrt (xv, ymax)
      call vector (xv, ymax-tickl)
   40    continue
   50    if (kydo.le.0)    return
      xv = -5.5_R8
      if (ymin.lt.xv)   xv = -11.0_R8
      if (ymin.lt.xv)   xv = -22.0_R8
      tickl = tscale * (xmax - xmin)
      do 80  k=1,20
      xv = xv + delx
      if (xv.gt.ymax)      return
      call setcrt (xmin,       xv)
      call vector (xmin+tickl, xv)
      call setcrt (xmax,       xv)
      call vector (xmax-tickl, xv)
   80    continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
