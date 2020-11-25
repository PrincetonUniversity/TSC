      subroutine sawdraw(s100,tmin,tmax,ymin,ymax)
!
      USE CLINAM, except_this => s100
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER is
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tmax,ymin,ymax,tmin,x1,y1,vx1,vy1
      character*80 s100(30)
!============
      call setold(tmax,ymin,1,0,1,0)
      write(s100,6666)
 6666 format(" time(sec)")
      call gtext(s100,80,0)
!
      if(isaw.ne.3 .or. numsaw .le.0) return
      do 200 is=1,numsaw
!
      if(sawtime(is) .lt. tmin .or. sawtime(is).gt.tmax)                 &  
     &   go to 100
      x1 = sawtime(is)
      y1 = ymin + 0.1_R8*(ymax-ymin)
      vx1 = 0
      vy1 = -0.1_R8*(ymax-ymin)
      call arrow(x1,y1,vx1,vy1)
  100 continue
  200 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
