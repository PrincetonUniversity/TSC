 
!     File: grids.f(or)                     ---------------------------|
!                                                                      |
      SUBROUTINE mkgrids
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL psigrids, vgrids
      call psigrids
!                                       generate psi grid
      call vgrids
!                                       generate parallel velocity grid
      call mkqlsm
!                                       generate smoothing function
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
