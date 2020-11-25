!
!     -----------------------------------------------------------------
!
      SUBROUTINE NoPlots
      USE params
      USE PlPr
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i
      do 10 i=1, nPlFlg
        PlFlg(i) = 0
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
