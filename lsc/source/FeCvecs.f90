!     -----------------------------------------------------------------
      SUBROUTINE FeCvecs
!                                       initialize collisional
!                                       diffusion and drag terms
!                                       Dcoll, nuColl with velocity bins
!                                       and psi bins.
      USE DqlBins
      USE FeBins
      USE params
      USE ProfBody
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8                                                             &  
     &        vthsq, vth3, vth5, harg, hvpar, vpnorm, v12
      INTEGER ip, iv
      REAL*8 TransVel, LogThing, TailT32, TailT12
      TailT12 = sqrt (TailTeps)
      TailT32 = sqrt (TailTeps) * TailTeps
      if (TailNeps * TailTeps .gt. 0.00_R8) then
        LogThing = (-log(TailNeps) -log(TailTeps)/2._R8) /               &  
     &             (1._R8-TailTeps)
      else
        LogThing = 1.0E6_R8
      endif
      TailVtrn = sqrt ( 2._R8* LogThing )
      do 20 ip = 1, npsi
         vthsq = vtherm(ip) * vtherm(ip)
         TransVel = vtherm(ip) * TailVtrn
         vth3 = vthsq * vtherm(ip)
         vth5 = vth3 * vthsq
         nu0psi(ip) = NeAry(ip) / vth3
         do 10 iv = 1, nv
!           v12 = 0.5 * (Vpar(iv) + Vpar(iv + 1))!! removed Dec93; replaced by
!           v12 = vpar(iv)                       !! this.  Why was it here??
            v12 = vpar(iv)
            vpnorm = (v12 / vtherm(ip))
            harg = (1._R8+  vpnorm * vpnorm)
            hvpar = 1._R8/ (harg * sqrt(harg))
!           Dcoll(iv, ip) = DcollNorm * nu0psi(ip) * vthsq * hvpar Enright did
!           nuColl(iv, ip) = nuNorm * nu0psi(ip) * hvpar * v12     replace this
 
            Dcoll(iv, ip) = DcollNorm * nu0psi(ip) * vthsq * hvpar *     &  
     &           LnlAry(ip) * BetZAry(ip)
            nuColl(iv, ip) = nuNorm * nu0psi(ip) * hvpar * v12 *         &  
     &           LnlAry(ip) * BetZAry(ip)
!
            if ( TailPeps .gt. 0.00_R8.and.                              &  
     &           abs(vpar(iv)) .gt. TransVel) then
                Dcoll (iv,ip) = Dcoll(iv,ip) * TailT12
                nuColl(iv,ip) = nuColl(iv,ip)* TailT32
            endif
!
10       continue
20    continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
