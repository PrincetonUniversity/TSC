!#include "f77_dcomplx.h"
!     ------------------------------------------------------------------
      SUBROUTINE rspwr(pwrlev)
      USE params
      USE RayBins
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ir, itor, ipol
      REAL*8 pwrlev
      REAL*8 AREAL

!     initalize power on rays
      ir = 1
      do 10 itor = 1, ntors
      do 10 ipol = 1, npols
         power(1, ir) = TotPwr/AREAL(npols) * Spec(itor) * pwrlev
         ir = ir + 1
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
