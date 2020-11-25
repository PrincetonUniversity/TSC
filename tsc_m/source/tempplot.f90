      subroutine tempplot
!
!
!
!
!
!.....plots surface source and sink terms
!
      USE CLINAM
      USE SAPROP
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tmin, tmax, dmin, dmax, xmin, xmax, xlabset, cmin, cmax
!  
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psit
!
      IF(.not.ALLOCATED(psit)) ALLOCATE( psit(ppsi), STAT=istat)
 
!============      
      if (istat .ne. 0) stop 'Allocation Error : tempplot  '
!============      
!
      if(dpsi.le.0) return
      npts = nx/2
      if(npsit.gt.1) npts = npsit-1
      
!
      do 100 n=2,npts+1

  100 psit(n) = sqrt((float(n-1))/(float(npsit-1)))

      tmin = 0._R8
      tmax = 0._R8
      dmin = 0._R8
      dmax = 0._R8
      cmin = 0._R8
      cmax = 0._R8
      xmin = 0._R8
      xmax = psit(npts+1)
      do 200 n=2,npts
      dmax = max(dmax,d2pi(n),d2pe(n))
      dmin = min(dmin,d2pi(n),d2pe(n))
      tmax = max(tmax,te(n),ti(n))
      tmin = min(tmin,te(n),ti(n))
      cmin = min(cmin,chiesec(n),chiisec(n))
      cmax = max(cmax,chiesec(n),chiisec(n))
  200 continue
      if(dmax.le.dmin) dmax = dmin+1._R8
      if(tmax.le.tmin) tmax = tmin+1._R8
!
      call mapg(xmin,xmax,dmin,dmax,.1_R8,.47_R8,.7_R8,1._R8)
      call tracec(1hi,psit(2),d2pi(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1he,psit(2),d2pe(2),npts,-1,-1,0._R8,0._R8)
      xlabset = xmin-(xmax-xmin)*.20_R8
      call setold(xlabset,dmin,1,0,1,1)
      write(s100,1001)
      call gtext(s100,80,0)
 1001 format("second derivatives")
!
      call mapg(xmin, xmax, tmin,tmax,.59_R8,.96_R8,.7_R8,1._R8)
      call tracec(1hi,psit(2),ti(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1he,psit(2),te(2),npts,-1,-1,0._R8,0._R8)
      call setold(xmax+(xmax-xmin)*.05_R8,tmin,1,0,1,1)
      write(s100,1002)
      call gtext(s100,80,0)
 1002 format("temperatures")
!
      call mapg(xmin, xmax, cmin,cmax,.1_R8,.47_R8,.1_R8,.6_R8)
      call tracec(1hi,psit(2),chiisec(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1he,psit(2),chiesec(2),npts,-1,-1,0._R8,0._R8)
      call setold(xmax+(xmax-xmin)*.05_R8,cmin,1,0,1,1)
      write(s100,1003)
      call gtext(s100,80,0)
 1003 format("chi values from hyper term")
!
 
      call frscj(6)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
