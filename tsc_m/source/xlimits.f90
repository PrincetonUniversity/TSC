      subroutine xlimits(imin,imax)
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!...........................................................
!
!.....find index of xprof
!
!...........................................................
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER imax,imin,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xprof,xx,xprofm
!============
      xprof = acoef(42)
      if(xprof.le.0) xprof = ccon
      do 40 i=3,nx
      imin = i
      xx = xary(i)
      if(xx.gt.xprof) go to 45
   40 continue
      imin = 3
   45 continue
      xprofm = acoef(43)
      if(xprofm.le.0) xprofm = alx
      do 48 i=3,nx
      imax = i
      xx = xary(i)
      if(xx.gt.xprofm) go to 47
   48 continue
      imax = nx
   47 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
