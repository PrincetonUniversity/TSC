!     -----------------------------------------------------------------
      INTEGER FUNCTION ivtabl(np)
      USE FeBins
      USE params
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iv
      REAL*8                                                             &  
     &     np, vv
      vv = 1._R8/ np
      iv = 2
 5    continue
      if( (iv .le. nv) .and. (vpar(iv) .lt. vv) ) then
         iv = iv + 1
         go to 5
      endif
!
!                                       v grid symmetery question:
!                                       use vpar nearby of smaller energy
      if( iv .gt. ivZero) then
        ivtabl = iv - 1
      else
        ivtabl = iv
      endif
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
