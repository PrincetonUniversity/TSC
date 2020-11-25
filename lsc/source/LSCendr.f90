!     ------------------------------------------------------------------
      SUBROUTINE LSCendr ( ErrMsg )
!     Trys to end ray tracing asap because an error is encountered
      USE params
      USE RayWrk
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
!============
      CHARACTER*(*) ErrMsg
      Lstop = 1
      iEndRy = iEndRy + 1
      call LSCwarn ( ErrMsg )
      return
!     ------------------------------------------------------------------
      ENTRY LSCbigK ( ErrMsg )
      Lstop = 1
      call LSCwarn ( ErrMsg )
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
