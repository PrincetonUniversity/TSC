!n
!     ------------------------------------------------------------------
!
      SUBROUTINE jnorm(ip)
      USE FeBins
      USE Jrf
      USE MKSetc
      USE params
      USE PIetc
      USE ProfBody
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip, iv
!     set up normalization, tabulate u lookup values
      gmrun = NeAry(ip) * LnlAry(ip) * ECOULB**4 / ELECMS**2 * CLIGHT *  &  
     &        4._R8*PI * 1.0E-14_R8
      if(EdcAry(ip) .ne. 0._R8)then
         vnorm = sqrt(abs(NeAry(ip) / EdcAry(ip)))   *                   &  
     &           sqrt(ECOULB**3 / ELECMS * CLIGHT **2 * LnlAry(ip) *     &  
     &           4._R8* PI ) * 1.0E-9_R8
         if( EdcAry(ip) .gt. 0._R8)then
            muplus  = +1._R8
            muminus = -1._R8
         else
            muplus  = -1._R8
            muminus = +1._R8
         endif
      else
         vnorm = vnmax
         muplus = 1._R8
         muminus = -1._R8
      endif
      nuRuna = gmrun / (vnorm * vnorm * vnorm)
      do 10 iv = 1, nv
        ugr(iv) = Vpar(iv) / vnorm * muminus
 10   continue
!
!
          ivrun = 0
      if (muminus .eq. -1) then
        do 20 iv = (nv+1)/2, nv-1
          if ( ugr(iv)  .le. -1._R8.and. ugr(iv-1)  .gt. -1._R8) then
            ivrun = iv
            goto 40
          endif
 20     continue
      else
        do 30 iv = 2,(nv+1)/2
          if ( ugr(iv)  .le. -1._R8.and. ugr(iv-1)  .gt. -1._R8) then
            ivrun = iv
            goto 40
          endif
 30     continue
      endif
!
 40   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
