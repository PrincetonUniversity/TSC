!
!----------------------------------------------------------------------
!
      SUBROUTINE EZnorm ( fin , fout , fmin , fmax , n , iCross0)
!     Takes a set of  n  values in fin and supplies a scaled array
!     in fout such that the maximum is unity, and supplies fmax as
!     that maximum value for information.
!
!     MAX ( fout ) =  1.
!     MAX ( fin  ) =  fmax  =  fmax * MAX ( fout )
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i, iCross0
      REAL*8    fin(n) , fout(n) , fmin , fmax , temp
 
      fmax = fin(1)
      fmin = fin(1)
      do 10 i=2,n
        if (fmax .lt. fin(i) ) fmax = fin(i)
        if (fmin .gt. fin(i) ) fmin = fin(i)
 10   continue
 
      if (abs(fmin) .gt. abs(fmax) ) then
        temp = fmax
        fmax = fmin
        fmin = temp
      endif
 
      if (fmax .eq. 0._R8) then
        fmax = 1._R8
        fmin = 0._R8
      endif
 
      do 20 i=1,n
        fout(i) = fin(i)/fmax
 20   continue
 
      if (fmax * fmin .lt. 0.0_R8) then
        iCross0 = 2
        else
        iCross0 = 1
      endif
 
      return
      END
 
!                                                                      |
!                                                                      |
!     pbxio.f(or)   ends        ---------------------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
