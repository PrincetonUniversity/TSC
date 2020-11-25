!#include "f77_dcomplx.h"
!
!     plsgrd.f(or)  begins                  ---------------------------|
!                                                                      |
!                                                                      |
 
      SUBROUTINE plsgrd
!     set up grids PsiGrd, xary, zary, xsv
!     isym: symmetry option,  0-no symmetry    1- up/down symmetry
      USE MKSetc
      USE params
      USE plcmx
      USE ProfBody
      USE RayWrk
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, j
      REAL*8                                                             &  
     &     TwoH, TwoK
      REAL*8 AREAL
 
      if (isym .eq. 1) then
            do  10 j=1,NumZh-1
                do  10 i=1,NumR
                   PsiGrd(i, j) = PsiGrd(i, NumZ-(j-1))
 10         continue
      endif
!
 
 
20    TwoH = (xary(NumR)-xary(1))/AREAL(NumR-1)
 
      TwoK = (zary(NumZ)-zary(1))/AREAL(NumZ-1)
      dr = TwoH
      dz = Twok
 
 
      Do 25 i = 2, NumR
      xary(i)  =  xary(1) + AREAL(i-1) * TwoH
25    continue
      Do 30 j = 2, NumZ
      zary(j)  =  zary(1) + AREAL(j-1) * TwoK
30    continue
 
      return
      END
 
 
!                                                                      |
!                                                                      |
!     plsgrd.f(or)  ends                    ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
