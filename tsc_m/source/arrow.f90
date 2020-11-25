      subroutine arrow(ax,ay,vx,vy)
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 ay,vx,vy,ax,pi,r,s,bx,by,theta,psi,ab,rab,sab,delx
      REAL*8 dely,xbar,ybar,cx,cy,dx,dy
!============
      data pi,r,s/ 3.14159265_R8,.2_R8,.1_R8/
!
!.....velocities are assumed scaled to coordinates
      bx = ax + vx
      by = ay + vy
      if(vx.eq.0 .and. vy.eq.0) return
      theta = atan2(vy,vx)
      psi = theta + pi*.5_R8
      ab = sqrt(vx**2 + vy**2)
      rab = r*ab
      sab = s*ab
      delx = sab*cos(psi)
      dely = sab*sin(psi)
      xbar = ax + (1._R8-r)*vx
      ybar = ay + (1._R8-r)*vy
      cx = xbar - delx
      cy = ybar - dely
      dx = xbar + delx
      dy = ybar + dely
      call setcrt(ax,ay)
      call vector(bx,by)
      call vector(cx,cy)
      call setcrt(dx,dy)
      call vector(bx,by)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
