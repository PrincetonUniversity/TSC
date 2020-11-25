 
      SUBROUTINE convolve(n, nsm, smvec, vout, vin)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, nsm, i, j, ilow, iup, jlow, jup
      REAL*8                                                             &  
     &     smvec(-nsm : nsm), vout(n), vin(n)
!
!     convolute source vin with smvec to yield vout
!
!     interior points
      ilow = nsm + 1
      iup = n - nsm
      jlow = - nsm
      jup = nsm
      do 20 i = ilow, iup
         vout(i) = 0._R8
         do 10 j = jlow, jup
            vout(i) = vout(i) + smvec(j) * vin(i - j)
 10      continue
 20   continue
!
!     lower boundary
      ilow = 1
      iup = nsm
      jlow = - nsm
      do 40 i = ilow, iup
         vout(i) = 0._R8
         jup = i - 1
         do 30 j = jlow, jup
            vout(i) = vout(i) + smvec(j) * vin(i - j)
 30      continue
 40   continue
!
!     upper boundary
      ilow = n - nsm + 1
      iup = n
      jup = nsm
      do 60 i = ilow, iup
         vout(i) = 0._R8
         jlow = n - i
         do 50 j = jlow, jup
            vout(i) = vout(i) + smvec(j) * vin(i - j)
 50      continue
 60   continue
!
      return
      END
!                                                                      |
!                                                                      |
!     ql.f(or)  ends                        ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
