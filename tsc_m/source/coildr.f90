      subroutine coildr(i,j,l)
!
!......draw a box for a coil
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,l,i,iindx,n
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sym,x1,x2,x3,x4,x5,x6,x7,x8,x9,xa,xb,xc,xd,z1,z2,z3,z4
      REAL*8 z5,z6,z7,z8,z9,za,zb,zc,zd,xcent,zcent,rad,x0,z0,thte
      REAL*8 dthte,xn,zn
!============
      if(l.eq.2) go to 20
      do 15 iindx=1,1+isym
      sym = 1._R8
      if(iindx.eq.2) sym = -1._R8
      if(imovie.ge.3 .and. imovie.lt.10) call colora("yellow")
      x1 = xary(i-1)
      x2 = xary(i)
      x3 = xary(i+1)
      x4 = xary(i-1)
      x5 = xary(i)
      x6 = xary(i+1)
      x7 = xary(i-1)
      x8 = xary(i)
      x9 = xary(i+1)
      xa = .25_R8*(x1+x2+x4+x5)
      xb = .25_R8*(x2+x3+x5+x6)
      xc = .25_R8*(x5+x6+x8+x9)
      xd = .25_R8*(x4+x5+x7+x8)
      z1 = sym*zary(j-1)
      z2 = sym*zary(j-1)
      z3 = sym*zary(j-1)
      z4 = sym*zary(j)
      z5 = sym*zary(j)
      z6 = sym*zary(j)
      z7 = sym*zary(j+1)
      z8 = sym*zary(j+1)
      z9 = sym*zary(j+1)
      za = .25_R8*(z1+z2+z4+z5)
      zb = .25_R8*(z2+z3+z5+z6)
      zc = .25_R8*(z5+z6+z8+z9)
      zd = .25_R8*(z4+z5+z7+z8)
      call setcrt(xa,za)
      call vector(xb,zb)
      call vector(xc,zc)
      call vector(xd,zd)
      call vector(xa,za)
      if(l.ne.0) go to 15
      if(imovie.ge.3 .and. imovie.lt.10) call colora("blue")
      call vector(xc,zc)
      call setcrt(xb,zb)
      call vector(xd,zd)
   15 continue
      return
!
!.....active feedback coil
   20 continue
      do 40 iindx=1,1+isym
      sym = 1._R8
      if(iindx.eq.2) sym = -1._R8
      if(imovie.ge.3 .and. imovie.lt.10) call colora("magenta")
      xcent = xary(i)
      zcent = sym*zary(j)
      rad = .25_R8*(xary(i+1)-xary(i-1))
      x0 = xcent+rad
      z0 = zcent
      thte = 0._R8
      dthte = 6.283185307_R8/20._R8
      call setcrt(x0,z0)
      do 30 n=1,20
      thte = thte + dthte
      xn = xcent + rad*cos(thte)
      zn = zcent + rad*sin(thte)
   30 call vector(xn,zn)
   40 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
