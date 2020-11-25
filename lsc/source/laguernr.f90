!
!     ------------------------------------------------------------------
!
      SUBROUTINE laguernr ( coef, degree, x, epsilon, polish )
!     Given the  DEGREE  and  DEGREE+1  complex COEF's of the polynomial
!     Sum COEF(i) X**(i-1) and given  epsilon the desired fractional
!     accuracy, and given a complex value  X , this routine improves  X
!     by Laguerre's method until it converges to a root of the given
!     polynomial.  For normal use  POLISH  should be input as 0 (false).
!     When  POLISH  is 1 (true) the routine ignores  EPSILON  and
!     instead attempts to improve  X  (assumed to be a good initial
!     guess) to the achievable roundoff limit.
!     Ref: Numerical Recipes in Fortran, page 264.
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER degree, iter, MAXIT, j, polish
      REAL*8                                                             &  
     &        epsilon, EPSS, err, abx, cdx
!     REAL*8    dxold                     ! Given in book; not used
      COMPLEX*16                                                         &  
     &        coef(degree+1), x, dx, x1,b,d,f,g,h,sq,gp,gm,g2, ZERO
      PARAMETER (ZERO=(0._R8,0._R8), EPSS=6.E-07_R8, MAXIT=100)
!
!     dxold = cabs ( x )                ! Gven in book; not used
!                                       Loop over iterations up to
!                                       allowed maximum
      do 12 iter = 1, MAXIT
        b   = coef ( degree+1 )
        err = abs ( b )
        d   = ZERO
        f   = ZERO
        abx = abs ( x )
!                                       Efficient computation of the
!                                       polynomial and its first 2
!                                       derivatives
        do 11 j = degree, 1, -1
            f   = x*f + d
            d   = x*d + b
            b   = x*b + coef ( j )
            err = abs ( b ) + abx * err
  11    continue
!                                       Estimate of roundoff error in
!                                       evaluating polynomial
        err = EPSS * err
        if ( abs ( b ) .le. err ) then
            return
!                                       The generic case: use Laguerre's
!                                       formula
        else
            g  = d / b
            g2 = g * g
            h  = g2 - 2._R8* f / b
            sq = sqrt ( (degree-1) * (degree*h - g2) )
            gp = g + sq
            gm = g - sq
            if ( abs (gp) .lt. abs (gm) ) gp = gm
            dx = degree / gp
        endif
        x1    = x - dx
        if ( x .eq. x1 ) return
        x     = x1
        cdx   = abs ( dx )
!       dxold = cdx                     ! Given in book; not used
        if ( polish .eq. 0 ) then
            if ( cdx .le. epsilon*abs ( x ) ) return
        endif
 
  12  continue
      call LSCstop( ' too many iterations in LAGUER root finder ')
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
