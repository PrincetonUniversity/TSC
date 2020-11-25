      subroutine setusr(xc1,xcm,yc1,ycn)
!
! * * set viewport coordinates
!
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL   xcm,yc1,ycn,xc1
!============
      call cpsetr('XC1', xc1)
      call cpsetr('XCM', xcm)
      call cpsetr('YC1', yc1)
      call cpsetr('YCN', ycn)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
