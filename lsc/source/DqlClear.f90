!     -----------------------------------------------------------------
      SUBROUTINE DqlClear
      USE DqlBins
      USE params
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      call BrodCast(NVELDIM * NPSIDIM, Dql, 0._R8)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
