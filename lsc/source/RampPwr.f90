!
 
 
 
 
!     ------------------------------------------------------------------
!     cycle.F begins                                                   |
!     -                                                                |
      SUBROUTINE Ramp_Pwr
      USE DqlBins
      USE FeBins
      USE params
      USE PlPr
      USE power_mod
      USE Ramppwr
      USE RayBins
      USE RayWrk
      USE tscunits
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iramp, ilowest, iPsiL, iPsiU, retcode
      INTEGER i1stcall, iCurious
      INTEGER iITRsave
      REAL*8    PraySumL(NRAMPDIM), PqlSumL(NRAMPDIM)
!     Local  values for debugging: PraySumL PqlSumL  ! L==Local
      EXTERNAL FeCalc
      DATA     i1stcall, iPsiL, iCurious                                 &  
     &        /       1,     1,    0 /
 
      call setpwrl(pwrlevel, nrampup, nflat)
      iPsiU = npsi
        i1stcall = 0
        iITR = 1
        call FeAt0
        call FePrime
        ilowest = 1
 
      do 10 iramp = ilowest, nrampup
!     .                                 rescale power level
         call rspwr(pwrlevel(iramp))
         call setiofl(iramp)
         call RayDamp
         call DqlGen
         iITR = mod(iITR,2) + 1
         call FeCalc(iPsiL, iPsiU)
         call FeWeight
         call FePrime
!        call FeDqlout ! take out June 2000
         if (iCurious .eq. 1) then
            call RfDamp
            call RfHeat
            PraySumL(iramp) = PraySum
            PqlSumL (iramp) = PqlSum
         else
            call PwrDiagn
         endif
!
         call FeConvrg(iramp)
!
 10   continue
      if (iCurious .eq. 1) then
!       call LSCPause
        write(nTSCscrn,'('' iSMO = '',i3,                                &  
     &  '' ; iramp, convergence figure, Pray, Pql: '' )') iSMO
        do 20 iramp = ilowest, nrampup
          write(nTSCscrn,'(i4, 3( 1x,1pe11.4) )')                        &  
     &        iramp, FeCvgAry(iramp),PraySumL(iramp),PqlSumL(iramp)
 20     continue
!       call LSCPause
      endif
!
!     The issue here is to put the original Maxwellian distribution
!     into the other iteration bin of fe so that a graph can be
!     constructed later showing the tail drawn out in fe.  See for
!     example this line in pbixio.F, SUBROUTINE Giruzzi
!       call Dql6Norm(Vpar(iv),fe(iv,ip,iOrigMaxwl),nv
!
        iITRsave = iITR
        iITR = mod(iITR,2) + 1
        call FeAt0
        iITR = iITRsave
!
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
