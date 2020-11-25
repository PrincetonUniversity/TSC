!#include "f77_dcomplx.h"
!
!     ------------------------------------------------------------------
!
      SUBROUTINE ddtNrnJrn (ip,nDot, jDot, vRunAwayIndex)
      USE DqlBins
      USE FeBins
      USE Jrf
      USE MKSetc
      USE params
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     Added Aug 93 to help quantify the runaway situation for LSC report
!     PPPL 2929
      INTEGER ip,iSMOi
      REAL*8    nDot, jDot, vRunAwayIndex
      REAL*8 AREAL
!
      iSMOi = mod(iSMO,2) + 1
 
        nDot= 0.00_R8
        JDot= 0.00_R8
        vRunAwayIndex= 0.00_R8
!     if (ivrun .gt. 1 .or.  ivrun .lt. nv ) then   BLUNDER repaired AUG 95
      if (ivrun .gt. 1 .and. ivrun .lt. nv ) then
        nDot =             Dql(ivrun, ip, iSMO) * dfdv(ivrun,ip,iSMOi)
        jDot = vpar(ivrun)*Dql(ivrun, ip, iSMO) * dfdv(ivrun,ip,iSMOi)
!
        nDot = abs ( nDot )
        jDot = abs ( JDot )* ECOULB * CLIGHT * 1.0E-5_R8
        vRunAwayIndex = AREAL(ivrun)
      endif
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
