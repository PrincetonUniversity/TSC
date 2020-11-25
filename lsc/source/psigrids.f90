!#include "f77_dcomplx.h"
!     ------------------------------------------------------------------
      SUBROUTINE psigrids
!                                       generate psi grid
      USE params
      USE ProfBody
      USE RayBins
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ugrid
      REAL*8                                                             &  
     &     TTINY, dpp
      INTEGER l
!, delpsi
      PARAMETER(TTINY = 1.E-5_R8)
      REAL*8 AREAL

      dpp = (psimaxx - psiminx) * TTINY / (AREAL(npsi - 1))
      call ugrid(PsiAry, npsi, psiminx + dpp, psimaxx - dpp)
!                                       ENFORCES INTERPOLATION !!
      delpsi = abs(PsiAry(2) - PsiAry(1))
      do 10 l=1,Npsi-1
        MidAry(l) = 0.5_R8*(PsiAry(l+1)+PsiAry(l))
 10   continue
        MidAry(Npsi) = PsiAry(Npsi)
      call Volcalc
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
