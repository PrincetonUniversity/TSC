!#include "f77_dcomplx.h"
!     ------------------------------------------------------------------
      SUBROUTINE mkqlsm
      USE DqlBins
      USE params
      USE PlPr
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iv
      REAL*8                                                             &  
     &        earg, denomql, qnorm
      REAL*8 AREAL
!     initialize quasilinear smoothing function
      nsmsym = (nsmoo + 1) / 2
      denomql = 1._R8/ AREAL(nsmw)
      do 10 iv = 1, nsmoo
         earg = AREAL(iv - nsmsym) * denomql
         earg = earg * earg / 2._R8
         qlsm(iv) = exp(-earg)
 10   continue
      call vecnorm(nsmoo, qlsm)
!                                       vecnorm: makes nsmoo values of qlsm
!                                       add to 1.
      call vsum(nsmoo, qlsm, qnorm)
!                                       vsum: adds nsmoo values of qlsm,
!                                       returns sum as qnorm.  Should be 1.
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
