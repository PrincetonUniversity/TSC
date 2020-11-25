      subroutine limits(yary,imin,imax,epsl,epsu,ymin,ymax)
!.....6.98 limits
!.....6.98 limits
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER imin,imax,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 epsl,epsu,ymin,ymax,yary,yval,diff
!============
      dimension yary(1)
!============      
!
      ymin   =  1.E20_R8
      ymax   = -1.E20_R8
      do 10 i = imin,imax
      yval   = yary(i)
      if(yval .lt. ymin) ymin = yval
      if(yval .gt. ymax) ymax = yval
   10 continue
!
      diff   = ymax - ymin
      ymin   = ymin - epsl*diff
      ymax   = ymax + epsu*diff
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
