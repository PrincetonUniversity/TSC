!     ------------------------------------------------------------------
      SUBROUTINE blockcpy(n, tovec, nto, frvec, nfrom)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i, nto, nfrom
      REAL*8                                                             &  
     &        tovec(nto, *), frvec(nfrom, *)
      do 10 i = 1, n
         tovec(1, i) = frvec(1, i)
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
