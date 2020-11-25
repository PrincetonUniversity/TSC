!     -----------------------------------------------------------------
      SUBROUTINE FeWeight
      USE FeBins
      USE params
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iv,ip,new,old
      REAL*8    CtrWeght
        CtrWeght = 1._R8- WeghtItr
        new = iITR
        old = mod(iITR,2) + 1
      do 10 iv = 1, nv
      do 10 ip = 1, npsi
        fe(iv,ip,new) = fe(iv,ip,new)*WeghtItr +                         &  
     &                  fe(iv,ip,old)*CtrWeght
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
