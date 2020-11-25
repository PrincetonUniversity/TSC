!
      SUBROUTINE FePrU
!                                       Fe Prime Unsmoothed
      USE DqlBins
      USE FeBins
      USE params
      USE ProfBody
      USE RayBins
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ips, iv
      EXTERNAL BrodCast
      do 20 ips = 1, npsi
        do 10 iv = ivZero +1, nv - 1
          dfdv(iv, ips,1) =+(fe(iv+1, ips, iITR) - fe(iv, ips, iITR)) *  &  
     &         dvip(iv)
 10     continue
!                                       v grid symmetery question:
!                                       make dfdv by looking out to high abs(v)
        do 11 iv = 2, ivZero -1
          dfdv(iv, ips,1) =-(fe(iv-1, ips, iITR) - fe(iv, ips, iITR)) *  &  
     &         dvip(iv)
 11      continue
         dfdv(ivZero, ips,1) = 0._R8
         dfdv(nv    , ips,1) = dfdv(nv - 1, ips,1)
         dfdv( 1    , ips,1) = dfdv(     2, ips,1)
!     compute unsmoothed derivative
!
 20   continue
      return
      END
!     fe    ends                            ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
