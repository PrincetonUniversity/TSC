      subroutine setvpt(vpl,vpr,vpb,vpt)
!
!  * * set viewport limits
!
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL   vpr,vpb,vpt,vpl
!============
      call cpsetr('VPL', vpl)
      call cpsetr('VPR', vpr)
      call cpsetr('VPB', vpb)
      call cpsetr('VPT', vpt)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
