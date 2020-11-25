!#include "f77_dcomplx.h"
!
!     ------------------------------------------------------------------
!
      SUBROUTINE zrootsnr ( coef, degree, roots, polish )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL laguernr
!     Given the DEGREE and the  DEGREE+1  complex COEF's of the polynomial
!     Sum COEF(i) x**(i-1) this routine successively calls LAGUER and
!     finds all DEGREE complex ROOTS.  The integer variable POLISH
!     should be input as 1 (true) if polishing is desired, of 0 (false)
!     if the roots will be subsequently polished by other means.
!     Ref: Numerical Recipes in Fortran, page 265.
!
!
 
      INTEGER degree, MAXDEGRE, i, j, jj, polish
      REAL*8                                                             &  
     &        EPSILON
      PARAMETER (EPSILON=1.E-06_R8, MAXDEGRE=101)
      COMPLEX*16                                                         &  
     &        coef(degree+1), x, b, c, roots(degree), defl(MAXDEGRE)
      REAL*8 AREAL
!
!                                       Copy coef's for successive deflation
      do 11 j = 1, degree+1
        defl(j) = coef(j)
  11  continue
!
      do 13 j = degree, 1, -1
!                                       Start at 0 to favor smallest
!                                       remaining root
        x = CMPLX ( 0._R8, 0._R8, R8)
        call laguernr ( defl, j, x, EPSILON, 0 )
        if ( abs( AIMAG(x) ) .le. 2._R8*EPSILON**2*abs( REAL(x) ) )     &  
     &      x = CMPLX ( REAL(x) , 0._R8, R8)
        roots( j ) = x
        b = defl( j+1 )
!                                       Forward deflation
        do 12 jj = j, 1, -1
            c = defl(jj)
            defl(jj) = b
            b = x*b + c
  12    continue
  13  continue
!
      if ( polish .eq. 1 ) then
!                                       Polish roots using undeflated coefs
        do 14 j = 1, degree
            call laguernr ( coef, degree, roots(j), EPSILON, 1 )
  14    continue
!
      endif
!
!                                       Sort roots by real part by
!                                       straight insertion
      do 16 j = 2, degree
        x = roots (j)
        do 15 i = j-1, 1, -1
            if ( REAL( roots(i) ) .le. REAL( x ) ) go to 10
            roots ( i+1 ) = roots ( i)
  15    continue
        i=0
  10    roots ( i+1 ) = x
  16  continue
      return
      END
!                                                                      |
!                                                                      |
!     zplrcnr ends      ------------------------------------------------
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
