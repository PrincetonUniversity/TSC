      subroutine setold(x,y,intensity,icase,isize,iorient)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER intensity,icase,isize,iorient
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 y,x
!============
      call setlch(x,y,icase,isize,iorient,-1)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
