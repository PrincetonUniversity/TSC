      subroutine second(tcpu)
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 tcpu,tarray
!============
      REAL*4 etime
      external etime
      dimension tarray(2)
      tcpu = etime(tarray)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
