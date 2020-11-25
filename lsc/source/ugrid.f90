!     matrix manipulation routines
!
      SUBROUTINE ugrid(vector,npts,vmin,vmax)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER npts, ipts
      REAL*8 vector, vmin, vmax, ripts1, rnpts1, dv
      DIMENSION vector(npts)
!
!     generate a uniform grid vector(j) from vmin to vmax
!
      if(npts.le.1)return
      rnpts1=npts-1
      dv=(vmax-vmin)/rnpts1
      do 1 ipts=1,npts
      ripts1=ipts-1
1     vector(ipts)=vmin+ripts1*dv
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
