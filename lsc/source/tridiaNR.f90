!
!     -----------------------------------------------------------------|
!
      SUBROUTINE tridiaNR (a,b,c, r,u,n,                                 &  
     &                              gam,ierror)
!     References:
!     Numerical Recipies in Fortran by William H. Press,
!     et al, Ch 2; p 40.
!
!     Solves for a vector  u  of length  n  the tridiagonal linear set:
!
!     b1   c1   0    0    0    0    0        u1          r1
!     a2   b2   c2   0    0    0    0        u2          r2
!     0    a3   b3   c3   0    0    0        u3          r3
!     0    0    a4   b4   c4   0    0        u4      =   r4
!     ..   ..   ..   ..   ..   ..   ..       ..          ..
!     ..   ..   ..   ..   aN-1 bN-1 cN-1     uN-1        rN-1
!     .              ..   0    aN   bN       uN          rN
!
!     a  b  c  r  are input vectors and are not modified.
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, j, ierror
      REAL*8    a(n), b(n), c(n), r(n), u(n)
      REAL*8    bet
!     PARAMETER (NMAX = 200)  ! as given in the book
!     REAL*8    gam(NMAX)       ! as given in the book
      REAL*8    gam(n)
!
!     Trap error on ill-posed problem:
!     if (b(1) .eq. 0.) stop ' b1 is zero in tridiag' ! as given in the book
      ierror = 0
      if (b(1) .eq. 0._R8) then
         ierror=1
         goto 1313
      endif
!
!     Normal beginning point:
      bet = b(1)
      u(1)= r(1)/bet
!
      do 11 j=2,n
        gam(j) = c(j-1)/bet
        bet    = b(j) - a(j)*gam(j)
!                if (bet .eq. 0.) stop ' tridiag fails' ! as given in the book
                 if (bet .eq. 0._R8) then
                    ierror=2
                    goto 1313
                 endif
        u(j)   = ( r(j) - a(j)*u(j-1))/bet
 11   continue
!
      do 12 j=n-1,1,-1
        u(j) = u(j)-gam(j+1)*u(j+1)
 12   continue
!
      return
!
!     Error condition...
 1313 continue
      return
!
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
