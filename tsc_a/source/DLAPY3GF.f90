      REAL*8 FUNCTION DLAPY3GF( X, Y )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   X, Y, Z
!     ..
!
!  Purpose
!  =======
!
!  DLAPY3GF returns sqrt(x**2+y**2+z**2), taking care not to cause
!  unnecessary overflow.
!
!  Arguments
!  =========
!
!  X       (input) DOUBLE PRECISION
!  Y       (input) DOUBLE PRECISION
!  Z       (input) DOUBLE PRECISION
!          X, Y and Z specify the values x, y and z.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL*8   ZERO
      PARAMETER          ( ZERO = 0.0_R8)
!     ..
!     .. Local Scalars ..
      REAL*8   W, XABS, YABS, ZABS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      Z = 0
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
         DLAPY3GF = ZERO
      ELSE
         DLAPY3GF = W*SQRT( ( XABS / W )**2+( YABS / W )**2+             &  
     &            ( ZABS / W )**2 )
      END IF
      RETURN
!
!     End of DLAPY3GF
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
