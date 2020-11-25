!     -----------------------------------------------------------------
      SUBROUTINE FeInit
!                                       initialize fe solver, fe
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ifirstcall
      EXTERNAL FeConst, FeArrays
      DATA    ifirstcall / 1/
 
      if (ifirstcall .eq. 1 ) then
        ifirstcall = 0
        call FeConst
      endif
 
      call FeArrays
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
