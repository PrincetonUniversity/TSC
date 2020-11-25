!     -----------------------------------------------------------------
      SUBROUTINE FeMaxw(iipsi)
!                                       contruct Maxwell distribution
!                                       at Psi(iipsi)
!                                       set fe at iipsi vs velocity
      USE FeBins
      USE MKSetc
      USE params
      USE PIetc
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iipsi, iv, DoCut6
      REAL*8    VthsInv, argument, TOOBIG
      REAL*8    Factor, FactorM1, TransVel, RsltMin, REXPMIN
      REAL*8    LoVfac, HiVfac, LoVterm, HiVterm
      DATA TOOBIG / 200._R8/
      DATA DoCut6/0/
!     DATA REXPMIN / 100. /           Alpha and Ted like 85 better than 100
      DATA REXPMIN /  85._R8/
!
!  dmc 13 Jan 1995 -- replaced "100." with REXPMIN in a couple of
!  places; this avoids compiler problem under osf/1 v3.0, f77 v3.6:
!  evaluation of exp( -100. ) gave an arithmetic error.  In the TRANSP
!  environment, however, it should be noted that exp( -REXPMIN )
!  evaluates by silent underflow to zero.
!
!     produce Maxwellian fe(ipsi) normalized to unity (MKS units)
!     TailPeps  Regarding a ficticious fast electron tail in both
!     TailNeps  directions, the fraction eps(ilon) for Pressure, Ne and Te
!     TailTeps  such that
!     .                  TailNeps = TailPeps * TailTeps
!     .                  TailTeps = T_thermal/T_fast
!     .                  TailPeps = (n_fast T_fast)/(n_thermal T_thermal)
!     .         and
!     .              f(v) = (2 pi v_t^2)^-.5 1. n_e exp[- v^2/(2 v_t^2)]; all v
!     .         f_fast(v) = (2 pi v_f^2)^-.5 1. n_f exp[- v^2/(2 v_f^2)]; all v
!     .
!     .         where n_f/n_e = TailNeps
!     .               v_t/v_f = TailTeps^.5
!     .               n_f v_f^2 / [ n_t v_t^2 ] = TailPeps
!     .
!     TailVtrn  is the transition velocity relative to v_t at which the fast
!               electron tail becomes more important
!                        TailVtrn = v/v_t | transition
!     .                           = sqrt[2 ln(1/Neps 1/Teps^.50)/(1 - Teps^2)]
!
      VthsInv = 1._R8/ (2._R8* Vtherm(iipsi) * Vtherm(iipsi))
      Factor  = 1._R8+ sqrt(TailTeps)*TailNeps
      FactorM1= sqrt(TailTeps)*TailNeps
!     fmxnorm = 1. / (sqrt(2. * PI) * vtherm(iipsi) * CLIGHT * 1.e8)
!
      if (FactorM1 .eq. 0.00_R8) then
      do 10 iv = 3, nv-2
!       fe(iv,iipsi,iITR) = fmxnorm*exp(-Vpar(iv)*Vpar(iv)*VthsInv)
        argument = Vpar(iv)*Vpar(iv)*VthsInv
        if (argument .lt. TOOBIG) then
          fe(iv,iipsi,iITR) = FeNorm(iipsi)* exp( -argument )
        else
          fe(iv,iipsi,iITR) = 0.00_R8
        endif
!
10    continue
!
        fe(1   ,iipsi,iITR) = 0.00_R8
        fe(2   ,iipsi,iITR) = 0.00_R8
        fe(nv-1,iipsi,iITR) = 0.00_R8
        fe(nv  ,iipsi,iITR) = 0.00_R8
      endif
!
      if (FactorM1 .ne. 0.00_R8) then
      RsltMin = exp ( - REXPMIN )
      TransVel = vtherm(iipsi) *  sqrt ( 2._R8*                          &  
     & (-log(TailNeps) -log(TailTeps)/2._R8) / (1._R8-TailTeps) )
!...
      fe(ivZero,iipsi,iITR) = 0.00_R8
      LoVfac = VthSInv
      HiVfac = LoVfac * TailTeps
      do  iv = ivZero - 1, 3, - 1
        if(abs(vpar(iv)) .lt. TransVel) then
          HiVterm = vpar(iv)*LoVfac
        else
          HiVterm = vpar(iv)*HiVfac
        endif
        if(abs(vpar(iv+1)) .lt. TransVel) then
          LoVterm = vpar(iv+1)*LoVfac
        else
          LoVterm = vpar(iv+1)*HiVfac
        endif
        fe(iv,iipsi,iITR) = fe(iv+1,iipsi,iITR) +                        &  
     &    (HiVterm + LoVterm)*(Vpar(iv)-Vpar(iv+1))
      enddo
!
      do  iv = ivZero + 1, nv-2,  1
        if(abs(vpar(iv)) .lt. TransVel) then
          HiVterm = vpar(iv)*LoVfac
        else
          HiVterm = vpar(iv)*HiVfac
        endif
        if(abs(vpar(iv-1)) .lt. TransVel) then
          LoVterm = vpar(iv-1)*LoVfac
        else
          LoVterm = vpar(iv-1)*HiVfac
        endif
        fe(iv,iipsi,iITR) = fe(iv-1,iipsi,iITR) +                        &  
     &    (HiVterm + LoVterm)*(Vpar(iv)-Vpar(iv-1))
      enddo
      do iv = 3, nv-2
        if( fe(iv,iipsi,iITR) .lt. REXPMIN ) then
          fe(iv,iipsi,iITR) = FeNorm(iipsi)*exp(-fe(iv,iipsi,iITR))
        else
          fe(iv,iipsi,iITR) = FeNorm(iipsi) * RsltMin
        endif
      enddo
      fe(1,iipsi,iITR) = 0.00_R8
      fe(2,iipsi,iITR) = 0.00_R8
      fe(nv-1,iipsi,iITR) = 0.00_R8
      fe(nv  ,iipsi,iITR) = 0.00_R8
      endif
 
      if (DoCut6 .eq. 1) then
         do iv=1,nv
            if (abs(vpar(iv)) .gt. 0.6_R8) fe(iv,iipsi,iITR) = 0.00_R8
         enddo
      endif
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
