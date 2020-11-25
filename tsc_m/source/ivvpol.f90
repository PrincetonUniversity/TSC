      subroutine ivvpol
!
!       Calculate DIAMAG.
!
!               Calculate poloidal current cvvpol in vacuum vessel due to
!               change in toroidal flux in plasma, i.e. change in
!               paramagnetism.
!
!               VVL, VVR = Inductance, Resistance of vacuum vessel for
!                          current in poloidal direction (henries, ohms)
!               CVVPOL   = poloidal current in amps at current time point
!               DIAMAG   = toroidal flux in webers relative to vacuum flux
!
!               L di/dt + R i = - d(Phi)/dt           circuit eqn.
!
!
      USE CLINAM
      USE FVV1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nivvpol,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sum2,fac,psimid,diamag,cvvpolhi,diamagol,cvvpold
      REAL*8 timesold
!============
        data  nivvpol /0/
!
        if (irst1.eq.1 .and. nivvpol.eq.0)   go to 95
!
      sum2 = 0._R8
      fac = 1.0_R8
      if(isym.eq.1 ) fac=2.0_R8
      do 90 i=iminn,imaxx
      do 90 j=jminn+1,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 90
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      if(psimid.ge.psilim)   go to 89
      sum2 = sum2 + fac*(g(i,j)-gzero/xsqoj(i))
   89 continue
   90 continue
!
        global(12) = sum2
        diamag     = sum2
!
   95   nivvpol = nivvpol + 1
        if (vvr.lt.1.E-9_R8)   return
        if (nivvpol.lt.1)   return
!
        if (nivvpol.eq.1)   then
                            if(irst1.ne.1)   cvvpol   = 0.0_R8
                            cvvpolhi = 0.0_R8
        else
        cvvpol = (diamagol - diamag + vvl*cvvpold)                       &  
     &              / (vvl + vvr * (times - timesold) )
                            endif
        if (abs(cvvpol).gt.abs(cvvpolhi))   cvvpolhi = cvvpol
        timesold = times
        cvvpold  = cvvpol
        diamagol = diamag
!***        if (nivvpol.lt.101)   then
!***                              write(nout,100)  times, diamag, cvvpol
!***        else
!***        if (mod(nivvpol,10).eq.0) write(nout,100) times, diamag, cvvpol
!***  100   format (' t, diamag, cvvpol =' f10.6, f10.4, e12.4)
!***                               endif
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
