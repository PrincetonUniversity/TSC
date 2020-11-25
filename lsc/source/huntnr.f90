!
!     ------------------------------------------------------------------
!
      SUBROUTINE huntnr( xx, n, x, jlo )
!     Given an array XX  of length  N , and given a value  X , returns
!     a value  JLO  such that  x  is between  XX(jlo)  and XX(jlo+1).
!     XX  must be monotonic, either increasing or decreasing.
!     JLO=0 or =N is returned to indicate that  X  is out of range.
!     JLO on input is taken as a the initial guess for  JLO on output.
!     Ref:  Numerical Recipes in Fortran, p 91.
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inc, jlo, jhi, jm, n
      REAL*8    xx(n), x
      LOGICAL ascnd
!
!                                       True if table in ascending order.
      ascnd =  xx(n) .gt. xx(1)
!                                       If input guess is not useful, then
!                                       go immediately to bisection.
      if ( jlo .le. 0 .or. jlo .gt. n ) then
        jlo = 0
        jhi = n + 1
        go to 3
      endif
 
!                                       Set the hunting increment.
      inc = 1
!                                       Hunt up.
      if ( x .ge. xx(jlo)  .eqv.  ascnd ) then
 1      jhi = jlo + inc
!                                       Done hunting; off the table.
        if ( jhi .gt. n ) then
            jhi = n + 1
!                                       Not done hunting.
        else if ( x .ge. xx(jhi)  .eqv.  ascnd ) then
            jlo = jhi
!                                       Double increment, try again.
            inc = inc + inc
            go to 1
!                                       Done hunting, value bracketed.
        endif
!                                       Hunt down.
      else
        jhi = jlo
 2      jlo = jhi - inc
!                                       Done hunting; off the table.
        if ( jlo .lt. 1 ) then
            jlo = 0
!                                       Not done hunting.
        else if ( x .lt. xx(jlo)  .eqv.  ascnd ) then
            jhi = jlo
!                                       So double increment, try again.
            inc = inc + inc
            go to 2
!                                       Done hunting, value bracketed.
        endif
      endif
!                                       Begin bisection phase.
 3    if ( jhi - jlo .eq. 1 ) return
      jm = (jhi + jlo) / 2
      if ( x .gt. xx(jm) .eqv. ascnd ) then
        jlo = jm
      else
        jhi = jm
      endif
      go to 3
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
