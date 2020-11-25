!
!     -----------------------------------------------------------------
!
      SUBROUTINE FeArrays
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL mkvth, FeMkNorm, Fecvecs
      call mkvth
!                                       generate thermal velocity array
!                                       (vthermal AND Vperpsq vs psi)
      call FeMkNorm
!                                       array of normalization for Fe
!                                       (Fe(v = 0, psi))
      call Fecvecs
!                                       initialize collisional vectors
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
