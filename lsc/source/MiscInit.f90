!#include "f77_dcomplx.h"
!     ------------------------------------------------------------------
      SUBROUTINE MiscInit
      USE FeBins
      USE params
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8    ONE, FLOAToNV
      REAL*8 AREAL
      DATA    ONE/                                                       &  
     &        1.0_R8/
      FLOAToNV = AREAL(nv)
      call ugrid(ilist, nv,ONE, FLOAToNV )
!     call ugrid(ilist, nv, 1., float(nv))
!     generate an integer list (stored in a floating array)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
