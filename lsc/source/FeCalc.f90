!     ------------------------------------------------------------------
!     TABLE OF CONTENTS:
!                       FeCalc(ipsiL, ipsiU) FeMkNorm
!                       FePlus(fe, nuCollx,Dcollx,Dqlx,Vpar,lgFe, ivZero,nv)
!                       FeMinus(...)
!                       FeInit wrFeCA FeArrays Fecvecs mkvth FeAt0
!                       FeMaxw(iipsi) FeConst FePrime
 
!     File:fe.f(or)                         ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE FeCalc(ipsiL, ipsiU)
!                                       compute electron distribution
!                                       function, given quasilinear
!                                       diffusion coefficient and
!                                       plamsa profiles from psi index
!                                       ipsiL to ipsiU
      USE DqlBins
      USE FeBins
      USE params
      USE ProfBody
      USE RayBins
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ipsiL, ipsiU, ip, iNotSmoo, iYesSmoo
      DATA iNotSmoo, iYesSmoo / 1 , 2 /
      EXTERNAL BrodCast, blockcpy,                                       &  
     &         FePlus, FeMinus
!
      call BrodCast(NPSIDIM * NVELDIM, fe(1,1,iITR), 0._R8)
!                                       set to zero to start
      call blockcpy(ipsiU - ipsiL + 1, fe(ivZero, ipsiL, iITR),          &  
     &        NVELDIM, Fenorm(ipsiL), 1)
!                                       copy normalization into
!                                       fe(v = vpar(ivZero), psi)
      do 10 ip = ipsiL, ipsiU
        call FePlus (fe(1, ip, iITR), nuColl(1, ip),                     &  
     &       Dcoll(1, ip), Dql(1,ip, iYesSmoo), Vpar, wkv, ivZero, nv)
!                                       solve for positive velocity
        call FeMinus(fe(1, ip, iITR), nuColl(1, ip),                     &  
     &       Dcoll(1, ip), Dql(1,ip, iYesSmoo), Vpar, wkv, ivZero, nv)
!                                       solve for negative velocity
10    continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
