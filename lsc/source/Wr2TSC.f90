!
      SUBROUTINE Wr2TSC(PelFnd, JrfFnd, PqlFnd)
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
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*46 TSCstring
      CHARACTER*47 Mystring
      CHARACTER*6  chrst6
      INTEGER  l, NumZeros, iion
      INTEGER iDBG, iamu, ichg
!     EQUIVALENCE (npsit, NpsiJ)
      REAL*8    PelFnd, JrfFnd, PqlFnd, PioFnd
      REAL*8    xlookup, yreturn, ydumy(NPSIDIM)
      REAL*8    powtsc (PPSI ), curtsc (PPSI), curtscp (PPSI),           &  
     &        powDql (PPSI ), dlnJdlnE(PPSI)
      REAL*8    powtscSm(PPSI)
      REAL*8                    dJdETRAN(PPSI)
      REAL*8    powtscI(PPSI ), curtscI(PPSI), curtscIp(PPSI),           &  
     &        powDqlI(PPSI )
!
      REAL*8    PowFac        , CurFac
      DATA    PowFac        , CurFac        /                            &  
     &            1.E-6_R8,     1.E-4_R8/
!             m^3 per cm^3  , m^2 per cm^2
!
      DATA    iDBG /0/
!
!     We write to TSC quantities of heating rate and current
!     driven based on the TSC \psi-grid.  These were accumulated on a
!     \psi-grid that is different.  It was made different so
!     that the heating rate and current driven stuff could have grids
!     concentrating on where the action was.
!     powtscSm  power to TSC Smoothed proportional to anecc*Jrf
!
!     We do the interpolation with linr1d (LINeaR 1D), which
!     does not introduce swings in the interpolant, such as where
!     the function becomes non zero after being flat at zero.
!
!cc   First integrate the power and current driven.  The points at which
!cc   integrals are known are between the points at which \psi is given.
!cc
!cc      int1      int2      int3      int4       int5  int6
!cc        |         |         |         |          |    |
!cc psi1      psi2      psi3      psi4       psi5      psi6
!cc   .         .         .         .          .         .
!
      call PJrfIgrl
 
      do 10 l = 1, NpsiJ - 1
         xlookup =  MidVec(l)
!
!     interpolate for power integral
         call linr1d(Npsi, MidAry, PrIntgrl,                             &  
     &               xlookup, yreturn)
         powtscI(l) =  abs( yreturn )
!        powtscI is the power in watts inside the l flux surface boundary
!
!     interpolate for current integral
         call linr1d(Npsi, MidAry, IrIntgrl,                             &  
     &               xlookup, yreturn)
         curtscI (l)=  yreturn
!        curtscI is the current in amps inside the l flux surface boundary
 
!     interpolate for current integral at enhanced E [ E plus dE]
         call linr1d(Npsi, MidAry, IpIntgrl,                             &  
     &               xlookup, yreturn)
         curtscIp(l)=  yreturn
!        curtscIp is the current in amps inside the l flux surface boundary
!
!     interpolate for the ql power deposited
         call linr1d(Npsi, MidAry, PqIntgrl,                             &  
     &               xlookup, yreturn)
         powDqlI(l) =  yreturn
 10   continue
 
         powtscI (NpsiJ) = PrIntgrl(Npsi)
         curtscI (NpsiJ) = IrIntgrl(Npsi)
         curtscIp(NpsiJ) = IpIntgrl(Npsi)
         powDqlI (NpsiJ) = PqIntgrl(Npsi)
 
         powtscI (NpsiJ-1) = PrIntgrl(Npsi)
         curtscI (NpsiJ-1) = IrIntgrl(Npsi)
         curtscIp(NpsiJ-1) = IpIntgrl(Npsi)
         powDqlI (NpsiJ-1) = PqIntgrl(Npsi)
 
!
!     divide by volume for power, area for current.
!     Note backward volume convention
           l = 1
           powtsc(l) = PowFac * (powtscI (l)             )/              &  
     &                           dVlVec(l)
           powDql(l) = PowFac * (powDqlI (l)             )/              &  
     &                           dVlVec(l)
           curtsc (l)= CurFac * (curtscI (l)             )/              &  
     &                         ( dVlVec(l) / (TWOPI * xmag) )
           curtscp(l)= CurFac * (curtscIp(l)             )/              &  
     &                         ( dVlVec(l) / (TWOPI * xmag) )
 
         do 20 l = 2, NpsiJ
 
           powtsc(l) = PowFac * (powtscI (l)-powtscI(l-1))/              &  
     &                           dVlVec(l)
           powDql(l) = PowFac * (powDqlI (l)-powDqlI(l-1))/              &  
     &                           dVlVec(l)
           curtsc (l)= CurFac * (curtscI (l)-curtscI (l-1))/             &  
     &                         ( dVlVec(l) / (TWOPI * xmag) )
           curtscp(l)= CurFac * (curtscIp(l)-curtscIp(l-1))/             &  
     &                         ( dVlVec(l) / (TWOPI * xmag) )
!          I have to mult by 2 pi R to make an area!
 
 20      continue
 
!
         do 30 l = 1, NpsiJ
           if(curtsc (l) .eq. 0.00_R8.or.                                &  
     &          dEdcAmnt .eq. 0.00_R8.or.                                &  
     &         EdcVec(l) .eq. 0.00_R8) then
             dlnJdlnE(l) = 0.00_R8
           else
             dlnJdlnE(l) = (curtscp(l) - curtsc (l))/curtsc (l)  *       &  
     &                                   EdcVec (l) / dEdcAmnt
           endif
!! TRANSP REQUEST:                                                           !
             dJdETRAN(l) = (curtscp(l) - curtsc (l))/ dEdcAmnt
 
 30      continue
!
!                                       Make sure we do not append information
      rewind (unit = nTSCread)
!
!                          powtsc is in watts/cc
!                          curtsc is in amps/cm**2
!                          powtscI is in watts
!                          curtscI is in amps
 
      PelFnd = powtscI(NpsiJ)
      JrfFnd = curtscI(NpsiJ)
      PqlFnd = powDqlI(NpsiJ)
 
      write(nTSCread,1011) (powtsc  (l),l=1,npsiJ-1)
      write(nTSCread,1012) (curtsc (l),l=1,npsiJ-1)
      write(nTSCread,1013) (dlnJdlnE(l),l=1,npsiJ-1)
      write(nTSCread,1014) (dJdETRAN(l),l=1,npsiJ-1)
!
!     compute the interpolated power for each ion, and then
!     write it out.  The outer loop on ions uses the array powtsc
!     over and over.
      PioFnd = 0._R8
      do iion = 1,nspc
         iamu = int (amass(iion)/1.6E-24_R8+ 0.5_R8)
         ichg = int (achrg(iion)         + 0.5_R8)
         write(chrst6,'(''Pamu'',i2)')iamu
         if      (iamu .eq. 1 .and. ichg .eq. 1) then
            chrst6='Pprot'
         else if (iamu .eq. 2 .and. ichg .eq. 1) then
            chrst6='Pdeut'
         else if (iamu .eq. 3 .and. ichg .eq. 1) then
            chrst6='Ptrit'
         else if (iamu .eq. 3 .and. ichg .eq. 2) then
            chrst6='Phe-3'
         else if (iamu .eq. 4 .and. ichg .eq. 2) then
            chrst6='Phe-4'
         else if (iamu .eq. 12 .or. ichg .eq. 6) then
            chrst6='Pcarb'
         else if (iamu .eq. 16 .or. ichg .eq. 8) then
            chrst6='Poxy '
         endif
         do l = 1, NpsiJ - 1
            powtscI(l) =  0.00_R8
         enddo
            powtscI (NpsiJ)   = 0.00_R8
            powtscI (NpsiJ-1) = 0.00_R8
               l = 1
               powtsc(l) = PowFac * (powtscI (l)             )/          &  
     &                              dVlVec(l)
            do l = 2, NpsiJ
               powtsc(l) = PowFac * (powtscI (l)-powtscI(l-1))/          &  
     &                               dVlVec(l)
             enddo
             if (iion .eq. 1) then
         write(nTSCread,1015) (powtsc  (l),l=1,npsiJ-1)
             else if (iion .eq. 2) then
         write(nTSCread,1016) (powtsc  (l),l=1,npsiJ-1)
             else if (iion .eq. nspc) then
         write(nTSCread,1017) (powtsc  (l),l=1,npsiJ-1)
             else
         write(nTSCread,1001) (powtsc  (l),l=1,npsiJ-1)
             endif
         PioFnd = PioFnd + powtscI(NpsiJ)
 
      enddo
 
      write(nTSCread,1003) npsiJ
      write(nTSCread,1001) (powDql (l),l=1,npsiJ-1)
      write(nTSCread,1004) PelFnd, PqlFnd, JrfFnd
      write(nTSCread,'(''EOF '')')
!
!                                       Make sure the reader does not hit EOF
      rewind (unit = nTSCread)
!
!                                       Try to trap an error where most of
!                                       the power and current values are zero.
!                                       This trap was installed because of
!                                       apparent intermittent returns of
!                                       zero current drive and heating.
!                                       This trap could be removed if a bug
!                                       is found.
      NumZeros = 0
      do 50 l=1,NpsiJ
        if (powtsc(l) .eq. 0.000_R8.and. curtsc(l) .eq. 0.000_R8)        &  
     &                   NumZeros = NumZeros + 1
 50   continue
      if (NumZeros .gt. ((npsiJ-1) - (npsiJ/6)) )                        &  
     &                   call LSCstop(' no rf current or heating')
!                                       End of zeros trap
!
!     .                                 Call FastFrac every time, but write it
!     .                                 out only if the proper flag is set.
!     .                                 New FastFrac logic at request of TRANSP
      call FastFrac(Npsi)
 
 1001 format(5e16.6)
 
 1011 format(5e16.6, t82,'Electron power')
 1012 format(5e16.6, t82,'Electron current')
 1013 format(5e16.6, t82,'d lnJ/d lnE')
 1014 format(5e16.6, t82,'dJ/dE')
 1015 format(5e16.6, t82,'1st  ion power')
 1016 format(5e16.6, t82,'2nd  ion power')
 1017 format(5e16.6, t82,'last ion power')
 
 1003 format('EOF',/,                                                    &  
     &       ' The above are powtsc in watts/cc ',/,                     &  
     &       ' and curtsc in amps/cm^2,',/,                              &  
     &       ' and d lnJ / d lnE which is dimensionless,',/,             &  
     &       ' and dJ / dE in amp/cm^2 / volt/m ,',/,                    &  
     &       ' with npsiJ = ',i5,/,                                      &  
     &       ' The following is PowDql, which is similar to powtsc.')
 1004 format('EOF',/,                                                    &  
     &' Total Power, Pel, Pql; and Total Current, Irf',/,                &  
     &       5e16.6 )
 
 
      do 60 l=1,Npsi,2
        ydumy(l) = vnormPos(l)
 60   continue
      do 61 l=2,Npsi,2
        ydumy(l) = vnormNeg(l)
 61   continue
      do 62 l=(Npsi*14)/15, Npsi
        ydumy(l)=0._R8
 62   continue
 63   continue
 
      if (iDBG .eq. TRUE) then
        call EZinit
        call Plt6Norm (xsv,    powtsc,   NpsiJ,  'Ptsc' , 1,' ')
        call Plt6Norm (MidVec, powtscI,  NpsiJ,'I Ptsc' , 2,' ')
        call Plt6Norm (MidAry, PrIntgrl, Npsi, 'I Plsc' , 3,' ')
        call Plt6Norm (xsv,    curtsc,   NpsiJ,  'Jtsc' , 4,' ')
        call Plt6Norm (MidVec, curtscI,  NpsiJ,'I Jtsc' , 5,' ')
        call Plt6Norm (MidAry, Irintgrl, Npsi, 'I Jlsc' , 6,' ')
        TSCstring = ' Prf Jrf TSC w integral; LSC '//EquTitle(81:96)
        Mystring = TSCstring // '$'
        call MkGrfLst(TSCstring)
        call EZwrit(50,700,Mystring,0,0)
        call EZfini(0,0)
      endif
 
      if(PlFlg(JRFPL).ge. TRUE      ) then
        call EZinit
        TSCstring = ' Ne Te Edc dJdE Itsc Run vs r'//EquTitle(81:96)
        Mystring = TSCstring // '$'
        call EZwrit(50,700,Mystring,0,0)
        call GNwrit(50,700,Mystring,0,0)
!
        call Plt6Norm (PsiAry, NeAry , Npsi,  'Ne '  ,  1,'.')
        call Plt6Norm (PsiAry, ydumy , Npsi,  'Run'  ,  2,'.')
        call Plt6Norm (MidVec, curtscI, NpsiJ,'Itsc' ,  3,' ')
        call Plt6Norm (PsiAry, TeAry , Npsi,  'Te '  ,  4,',')
        call Plt6Norm (PsiAry, EdcAry, Npsi,  'Edc ' ,  5,' ')
        call Plt6Norm (xsv,  dlnJdlnE, NpsiJ, 'LdJdE',  6,' ')
!
        call MkGrfLst(TSCstring)
        call EZfini(0,0)
      endif
 
      if(PlFlg(PITPRFPL) .ge. TRUE ) then
        call PitchPrf
 
      endif
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
