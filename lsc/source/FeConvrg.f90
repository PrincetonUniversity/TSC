!#include "f77_dcomplx.h"
!     -----------------------------------------------------------------
      SUBROUTINE FeConvrg(iramp)
      USE FeBins
      USE params
      USE Ramppwr
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip,iv, iramp, inew, iold
      REAL*8    dum, epss
      REAL*8 AREAL
      DATA    epss / 1.E-32_R8/
      inew = iITR
      iold = mod(inew,2) + 1
      dum = 0._R8
      do 20   ip=1,Npsi
        do 10 iv=2,nv-1
          dum = dum + abs(fe(iv,ip,inew)-fe(iv,ip,iold)) /               &  
     &                   (fe(iv,ip,inew)+fe(iv,ip,iold) + epss)
 10     continue
 20   continue
      FeCvgAry(iramp) = dum/AREAL(nv*npsi)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
