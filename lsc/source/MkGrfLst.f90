!     ------------------------------------------------------------------
      SUBROUTINE MkGrfLst ( GrfTitle )
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*(*) GrfTitle
!     write(ns30c1,'(a30 )') GrfTitle ! Abandon list of graph titles
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
