 
!     power.F ---------------------------------------------------------+
!                                                                      |
      SUBROUTINE PdepCalc
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      call RfDamp
      call RfHeat
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
