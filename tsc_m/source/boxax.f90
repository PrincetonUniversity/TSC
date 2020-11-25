      subroutine boxax (xmin,xmax, ymin,ymax, ixax, jyax)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ixax,jyax
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xmax,ymin,ymax,xmin,div
!============
      character*4  ixaxfmt
      dimension div(2)
!============      
      call setcrt(xmin,ymin)
      call vector(xmax,ymin)
      call vector(xmax,ymax)
      call vector(xmin,ymax)
      call vector(xmin,ymin)
      if (ixax.le.0)   go to 20
      div(1) = xmin
      div(2) = xmax - 1.E-6_R8*(xmax-xmin)
!cccc    print *, ' boxax: xmin,xmax,ymin,ymax, div = ', xmin,xmax,
!cccc     &        ymin,ymax, div(1), div(2)
      ixaxfmt = 'f4.1'
      if (xmax .lt. 10.0_R8)   ixaxfmt = 'f3.1'
      call gaxisf (xmin,ymin, xmax,ymin, 0,1,0, ixaxfmt, 0, div)
   20    if (jyax.le.0)   return
      div(1) = ymin + 1.E-6_R8*(ymax-ymin)
      div(2) = ymax
      call gaxisf (xmin,ymin, xmin,ymax, 0,1,1, 'f4.1', 0, div)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
