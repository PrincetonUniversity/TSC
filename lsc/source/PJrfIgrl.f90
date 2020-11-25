 
!
!
!                                                                      |
!                                                                      |
!     Raytrace.F  ends                      ---------------------------|
!
!     OLD code from BounShft saved Apr 93 temporarily
!     REAL*8    B2, Bp, crr, crz, crp, czr, czz, czp, cpr, cpz, cpp
!     REAL*8    KrNew, KzNew, KpNew
!     Go back to the last ok solution point, and re-arrange kr and kz
!     such that k \dot B is enhanced by "NparEnhc" (of order 1.1)
!     and a new k \cross B changes sign, and increases similarly;
!     put this new information into the starting condition for y.
!     This is BOUNcing the ray with a SHiFT of n_parallel
!
!     Vector algebra shows that, in order to preserve k_\parallel
!     and reverse the sign of k_\perp, an old vector k becomes a new
!     vector k' are related as follows:
!
!         \vec{k'} = - \vec{k} + 2 \vec{B} (k \dot B)/B^2
!
!     The integration vector is given by
!       y(1)  y(2)  y(3)  y(4)  y(5)  y(6)
!        r     z    phi   k_r   k_z    n
!
!      Bp  = RBphi/r
!      B2  =  Br*Br + Bz*Bz + Bp*Bp
!
!      crr = (Br*Br - Bz*Bz - Bp*Bp)/B2
!      czz = (Bz*Bz - Br*Br - Bp*Bp)/B2
!      cpp = (Bp*Bp - Br*Br - Bz*Bz)/B2
!      crz = 2. * Br*Bz/B2
!      crp = 2. * Br*Bp/B2
!      czp = 2. * Bz*Bp/B2
!      czr = crz
!      cpr = crp
!      cpz = czp
!
!      KrNew = crr * yok(4) + crz * yok(5) + crp * yok(6)/r
!      KzNew = czr * yok(4) + czz * yok(5) + czp * yok(6)/r
!      KpNew = cpr * yok(4) + cpz * yok(5) + cpp * yok(6)/r
!
!      KrNew = NparEnhc * KrNew
!      KzNew = NparEnhc * KzNew
!      KpNew = NparEnhc * KpNew
!     end temp saved code
 
 
 
 
 
 
!
!     Wr2TSC    begins                      ---------------------------|
!     -                                                                |
!     -                                                                |
      SUBROUTINE PJrfIgrl
!     Computes the integrals of Prf, Pql, Jrf, as function of psi
      USE Jrf
      USE params
      USE PIetc
      USE PlPr
      USE power_mod
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE TSCgrap
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
      CHARACTER*46 TSCstring
      CHARACTER*47 Mystring
      INTEGER ips
      REAL*8    nRdotNorm(NPSIDIM)
      REAL*8    jRdotNorm(NPSIDIM)
!
!     Integrate Pray and Pql in LSC flux surfaces
      PrIntgrl(1) = PRayTot(1)
      PqIntgrl(1) = PqlTot (1)
      do ips = 2, Npsi
         PrIntgrl(ips) = PrIntgrl(ips-1) + PRayTot(ips)
         PqIntgrl(ips) = PqIntgrl(ips-1) + PqlTot(ips)
      enddo
!
!     Integrate Jrf and Jrf-plus in LSC flux surfaces
!     Divide by 2 pi R to make an area element !
      IrIntgrl(1) = js (1) * dVol(1) / (TWOPI * xmag )
      IpIntgrl(1) = jsp(1) * dVol(1) / (TWOPI * xmag )
      do ips = 2, npsi
         IrIntgrl(ips) = IrIntgrl(ips-1) + js(ips)*dVol(ips) /           &  
     &                                    (TWOPI * xmag )
         IpIntgrl(ips) = IpIntgrl(ips-1) +jsp(ips)*dVol(ips) /           &  
     &                                    (TWOPI * xmag )
      enddo
 
      do ips=1,Npsi
         if (PrIntgrl(ips) .le. 0.05_R8* PrIntgrl(Npsi) .or.             &  
     &       PrIntgrl(ips) .ge. 0.95_R8* PrIntgrl(Npsi)     ) then
             nRunDot(ips) = 0.00_R8
             jRunDot(ips) = 0.00_R8
         endif
         nRdotNorm(ips) = nRunDot(ips)/(NeAry(ips) + 1.E10_R8)
         jRdotNorm(ips) = jRunDot(ips)/(   js(ips) + 1.E01_R8)
      enddo
 
      if (PlFlg(RFDPSPL) .ge. TRUE) then
        call EZinit
        call GNinit(nTSCgraf)
!
        call Plt6Norm (PsiAry, PRayTot,  Npsi,   'Pray' , 1,' ')
!       call Plt6Norm (PsiAry, nRunDot,  Npsi ,'nRdot ' , 2,' ')
        call Plt6Norm (PsiAry, PqlTot,   Npsi,   'Pql ' , 2,' ')
! replc call Plt6Norm (PsiAry, PqlTot,   Npsi,   'Pql ' , 2,' ')
        call Plt6Norm (PsiAry, js,       Npsi ,  'Jql ' , 3,' ')
        call Plt6Norm (MidAry, PrIntgrl, Npsi, 'I Pray' , 4,' ')
! replc call Plt6Norm (PsiAry, jRunDot,  Npsi, 'jRdot ' , 5,' ')
! replc call Plt6Norm (MidAry, PqIntgrl, Npsi, 'I Pql ' , 5,' ')
        call Plt6Norm (MidAry, PqIntgrl, Npsi, 'I Pql ' , 5,' ')
!       call Plt6Norm (PsiAry, vRunIdx,  Npsi, 'ivRun ' , 5,' ')
        call Plt6Norm (MidAry, IrIntgrl, Npsi, 'I Jql ' , 6,' ')
!
!       TSCstring = ' Prf nRd Jrf w integral; LSC  '//EquTitle(81:96)
        TSCstring = ' Prf Pql Jrf w integral; LSC  '//EquTitle(81:96)
        Mystring = TSCstring // '$'
        call MkGrfLst(TSCstring)
        call EZwrit(50,700,Mystring,0,0)
        call EZfini(0,0)
        call GNfini
 
      endif
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
