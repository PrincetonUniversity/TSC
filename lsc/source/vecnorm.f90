!     ------------------------------------------------------------------
      SUBROUTINE vecnorm(n, vec)
!                                       Makes n elements of vec add to 1,
!                                       returns the necessary divisor as norm
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n
      REAL*8    vec(n), norm, cnorm
      call vsum(n, vec, norm)
!                                       Adds n elements of vec, returns sum as
!                                       norm
      cnorm = 1._R8/ norm
      call svmult(n, vec, cnorm, vec)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
