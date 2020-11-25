      subroutine gtextm(string,ic,ioff,imin,imax)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ioff,imin,imax,ic,imaxx,i
!============
      character*(*) string(*)
      imaxx = imax
      if(imaxx.gt.30) imaxx=30
      do 10 i=imin,imaxx
      call gtext(string(i),ic,ioff)
   10 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
