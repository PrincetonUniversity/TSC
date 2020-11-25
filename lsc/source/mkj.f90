!
!     ------------------------------------------------------------------
!
      SUBROUTINE mkj(ip, jd, iGotRuna, iFillJray)
      USE DqlBins
      USE FeBins
      USE Jrf
      USE MKSetc
      USE params
      USE power_mod
      USE ProfBody
      USE WkAry
      USE PIetc
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip, iv, iGotRuna, iWhichWay, iSMOi, iFillJray
      REAL*8    jd, jdtemp, constfac, xr, eps, w
      REAL*8    FeCutFac, FeCutOff, TrapFac
      DATA    FeCutFac / 1.0E-5_R8/
      REAL*8    dWsdu, dWsduou
!     compute jrf driven, jd, at flux surface ip
      jd = 0._R8
      constfac =   ECOULB * CLIGHT * 1.0E-5_R8/ nuRuna
      iGotRuna = 0
      iSMOi = mod(iSMO,2) + 1
      FeCutOff = FeCutFac*fe(IvZero,ip,iITR)
      eps = sqrt(iVlVec(ip)/(2.*PI**2*xmag))/xmag
      xr = 3.5
!
!     .
!     .                                 first negative velocities
!     .
      do 10 iv = 1, IvZero - 1
!     .                                 Ignore jd where f_e is small
        if ( fe(iv,ip,iITR) .ge. FeCutOff ) then
!     .                                 Muplus, because mu is plus for
!     .                                 negative velocity if the Edc is >0
!     .                                 See above for reversal if E<0
          call WsloPrm(ugr(iv), muplus , ZbrAry(ip),                     &  
     &                 dWsduou, iWhichWay)
          dWsdu = dWsduou  * ugr(iv)
!     .                                 Note that dWsdu*ugr is always positive
          if (iWhichWay .ne. 0 ) iGotRuna = iGotRuna + 1
!     .
!     .                                 trapping effect from Ehst-Karney
!     .                                 for LH limit of large vpar/vth
      w = abs(Vpar(iv))/Vtherm(ip)
      TrapFac = 1.-((eps**0.77*sqrt(xr**2+w**2))/(eps**0.77*xr+w))
!     .                                 For negative velocity the current
!     .                                 driven is negative because df/dv >0
!
          jdtemp      = - Dql(iv, ip, iSMO) * dfdv(iv,ip,iSMOi) *        &  
!cj  &                                         dvplus(iv) *  dWsdu *     &  
     &                          TrapFac *      dvplus(iv) *  dWsdu *     &  
     &                    constfac
          jd = jd + jdtemp
          if(iFillJray .eq. 1) Jray(iv,ip) = jdtemp
        endif
 10   continue
!     .
!     .                                 then positive velocities
!     .
      do 20 iv = IvZero + 1, nv - 1
!     .                                 Ignore jd where f_e is small
!     .                                 Muminus, because mu is neg for
!     .                                 positive veloctiy if the Edc is >0
!     .                                 See above for reversal if E<0
        if ( fe(iv,ip,iITR) .ge. FeCutOff ) then
          call WsloPrm(ugr(iv), muminus , ZbrAry(ip),                    &  
     &                 dWsduou, iWhichWay)
          dWsdu = dWsduou  * ugr(iv)
!     .                                 Note that dWsdu*ugr is always positive
          if (iWhichWay .ne. 0 ) iGotRuna = iGotRuna + 1
!     .
!     .                                 trapping effect from Ehst-Karney
!     .                                 for LH limit of large vpar/vth
      w = abs(Vpar(iv))/Vtherm(ip)
      TrapFac = 1.-((eps**0.77*sqrt(xr**2+w**2))/(eps**0.77*xr+w))
!     .                                 For positive velocity the current
!     .                                 driven is positive because df/dv<0
!
          jdtemp      = - Dql(iv, ip, iSMO) * dfdv(iv,ip,iSMOi) *        &  
!cj  &                                         dvplus(iv) *  dWsdu *     &  
     &                          TrapFac *      dvplus(iv) *  dWsdu *     &  
     &                    constfac
          jd = jd + jdtemp
          if(iFillJray .eq. 1) Jray(iv,ip) = jdtemp
        endif
 20   continue
!      use of 'constfac' removes the need for this Oct 94
!      jd =  jd / nuRuna * ECOULB * CLIGHT * 1.0e-5
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
