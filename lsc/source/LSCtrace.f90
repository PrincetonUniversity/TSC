!     ------------------------------------------------------------------
      SUBROUTINE LSCtrace ( ErrMsg )
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*(*) ErrMsg
      write(nTSCscrn,'('' LSCtrace called on exiting: '',a40 )')         &  
     &                     ErrMsg
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
