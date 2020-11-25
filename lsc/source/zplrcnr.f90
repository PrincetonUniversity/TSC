!#include "f77_dcomplx.h"
 
!     zplrcnr begins    ------------------------------------------------
!                                                                      |
!                                                                      |
      SUBROUTINE zplrcnr (degree, coefr, zeroc)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL zrootsnr
!     ZPLRCNR:
!     finds Zeros of Polynomials with Laguerre's method assuming
!     Real Coefficients by calling a routine from
!     Numerical Recipes.
!     The calling convention is set up to look like IMSL call ZPLRC.
!     f(x) = 0 = coefr(1) + coefr(2)*x + coefr(3)*x**2 ...etc
      INTEGER degree, DEGMAX, i, polish
      PARAMETER (DEGMAX=10)
      REAL*8    coefr(degree+1)
      COMPLEX*16 coefc(DEGMAX+1), zeroc(degree)
      REAL*8 AREAL
!
!     if (degree .gt. DEGMAX ) HALT ' out of bounds in ZPLRCNR '
!
!                                       zroots needs complex coef's
      do 10 i = 1, degree+1
      coefc(i) = CMPLX ( coefr(i) , 0._R8, R8)
  10  continue
!
      polish = 1
      call zrootsnr ( coefc, degree, zeroc, polish )
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
