 
!     RdTSC.f(or)  begins                  ----------------------------|
!                                                                      |
!                                                                      |
 
      SUBROUTINE rdTSC
!     read equilibrium from TSC generated disk file
      USE Doflags
      USE Escan
      USE MKSetc
      USE params
      USE PIetc
      USE plcmx
      USE ProfBody
      USE RayWrk
      USE TSCgrap
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit REAL declarations:
!============
      INTEGER i, j, n, l, npsiPut
      REAL*8    voltPut(*)
      CHARACTER*8  name(10), timhd(2)
      REAL*8    psderiv(0:2, 0:2), psimat(0:3, 0:3)
      COMMON / plxxr / psderiv, psimat
      REAL*8   nePeak, TePeak, BoPeak, neEdge
      REAL*8   neRatio,TeRatio,BoRatio, IpRatio
      REAL*8   ZeroNe, ZeroTe
      REAL*8   PsiVal, ssqrt2, TBTemp
      DATA   ZeroNe, ZeroTe / 1.0E+09_R8, 1.0E-02_R8/
!     pe2Fac pi2Fac convert density in ^14 cm^-3 or ^20 m^-3
!     and frequency in GHz into
!     \omega_{pe}^2 / \omega^2 = pe2Fac * n_e[14]/f^2[GHz]
!     \omega_{pi}^2 / \omega^2 = pi2Fac * n_i[14]/f^2[GHz]
!     pe2Vec is a quantity that must be divided by fghz^2 to get:
!                                                     \omega_{pe}^2 / \omega^2
!
!     AelFac AioFac convert density,temperature,mass to thermal terms
!     OmcFac converts total magnetic field to GHz frequency
!     ceiFac convert B , ion density and charge to
!     \omega_{ce} \omega_{ci} /f^2[GHz]
!     TBTemp introduced by TB Terpstra to avoid a divide by zero with
!     the Alpha Open VMS systems in Feb 94
 
      pe2Fac    =( ECOULB * CLIGHT )**2/ ( PI * ELECMS ) * 1.E+4_R8
      pe2Fac14    =( ECOULB * CLIGHT )**2/ ( PI * ELECMS ) * 1.E-10_R8
      pi2Fac    = pe2Fac * (ELECMS/PROTMS)               * 1.E-4_R8
      AioFac    = 3._R8* pi2Fac * ECOULB/PROTMS/TWOPI**2 * 1.E-7_R8
      AelFac    = 0.75_R8* pe2Fac * ECOULB/ELECMS/TWOPI**2 * 1.E-3_R8
      OmcFac    =                 ECOULB/ELECMS/TWOPI    * 1.E+3_R8
      ceiFac    =( ECOULB/TWOPI )**2/ELECMS/PROTMS       * 1.E+2_R8
!
!
!                                       Make sure we do not hit EOF
      rewind (unit = nTSCwrit)
!
      read(nTSCwrit,2201,end=98,err=99) (equhd(i),i=1,10)
      read(nTSCwrit,2202,end=98,err=99) npsitm,nspc,kcycle,times
!
            write(timhd( 1),'('' time = '')')
            write(timhd( 2),'(f8.4)') times
            EquTitle = equhd(1)//equhd(2)//equhd(3)//                    &  
     &                 equhd(4)//equhd(5)//equhd(6)//                    &  
     &                 equhd(7)//equhd(8)//equhd(9)//equhd(10)//         &  
     &                 timhd(1)//timhd(2)
 
            NpsiJ    =   npsitm
            NpsiM1  =   npsitm-1
            if (NpsiJ .gt. PPSI)                                         &  
     &                    call LSCstop (' NpsiJ .gt. PPSI on read-in')
            if (nspc .gt. PIMP+1)                                        &  
     &                    call LSCstop (' nspc .gt. PIMP+1 on read-in')
            if (iError .gt. 0 ) return
!
      read(nTSCwrit,2203,end=98,err=99)( anecc(l),l=1,NpsiJ)
      read(nTSCwrit,2203,end=98,err=99)( tekev(l),l=1,NpsiJ)
      read(nTSCwrit,2203,end=98,err=99)((anicc(l,n),l=1,NpsiJ),n=1,nspc)
      read(nTSCwrit,2203,end=98,err=99)((tikev(l,n),l=1,NpsiJ),n=1,nspc)
            anecc(NpsiJ)  = min(anecc(NpsiJ),  ZeroNe)
            tekev(NpsiJ)  = min(tekev(NpsiJ),  ZeroTe)
            anicc(NpsiJ,1)= min(anicc(NpsiJ,1),ZeroNe)
            tikev(NpsiJ,1)= min(tikev(NpsiJ,1),ZeroTe)
      read(nTSCwrit,2203,end=98,err=99) (amass(n),n=1,nspc)
      read(nTSCwrit,2203,end=98,err=99) (achrg(n),n=1,nspc)
!
      read(nTSCwrit,2203,end=98,err=99) (voltlp(l),l=1,NpsiJ)
      read(nTSCwrit,2203,end=98,err=99) (rho(l),l=1,NpsiJ)
      read(nTSCwrit,2203,end=98,err=99) (vptemp(l),l=1,NpsiJ)
      read(nTSCwrit,2203,end=98,err=99) (xsv(l),l=1,NpsiJ)
            psiminx = xsv(1)
            psimaxx = xsv(NpsiJ)
      read(nTSCwrit,2203,end=98,err=99) (pary(l),l=1,NpsiJ)
      read(nTSCwrit,2203,end=98,err=99) (ppary(l),l=1,NpsiJ)
      read(nTSCwrit,2203,end=98,err=99) (gary(l),l=1,NpsiJ)
      read(nTSCwrit,2203,end=98,err=99) (gpary(l),l=1,NpsiJ)
      read(nTSCwrit,2204,end=98,err=99) nx,nz,isym,iplim
              if (nx  .gt. PNX)                                          &  
     &                call LSCstop (' nx .gt. PNX on read-in')
      read(nTSCwrit,2203,end=98,err=99)xary(1),xary(nx),zary(1),zary(nz)
!
            if (isym .eq. 0) then
!             call LSCwarn (' isym must be 1 for up/dn symmetry')
              if (nz  .gt. PNZ)                                          &  
     &                call LSCstop (' nz .gt. PNZ on read-in')
              NumR    =   nx
              NumZh   =  (nz+1)/2
              NumZ    =   nz
              Rmin    =   xary(1)
              Rmax    =   xary(NumR)
              Zmin    =   zary(1)
              Zmax    =   zary(NumZ)
            else if (isym .eq. 1) then
              NumR    =   nx
              NumZh   =   nz
              NumZ    = 2*nz - 1
              if (NumZ.gt. PNZ)                                          &  
     &                call LSCstop (' NumZ .gt. PNZ on read-in')
              zary(NumZ )=zary(nz)
              zary(NumZh)=zary(1)
              zary(1)    =-zary(NumZ)
              Rmin    =   xary(1)
              Rmax    =   xary(NumR)
              Zmin    =  -zary(NumZ)
              Zmax    =   zary(NumZ)
 
            else
              call LSCstop (' isym must be 1 for up/dn symmetry, or 0')
              return
            endif
!
      if (iError .gt. 0 ) return
!
 
!
!
      read(nTSCwrit,2203,end=98,err=99) psimin,psilim,xmag,zmag
      read(nTSCwrit,2203,end=98,err=99) rgzero,bgzero,apl
!
            Rmaj    =  rgzero
            Rmag    =  xmag
            Btesl   =  bgzero
            RBphi0  =  rgzero * bgzero
!
              l = 1
              dVlVec(l) = vptemp(l)* (+0.5_R8*(xsv(l+1) + xsv(l))        &  
     &        - psimin )
              iVlVec(l) = dVlVec(l)
              MidVec(l)= 0.5_R8*(xsv(l+1)+xsv(l))
            do 6 l=2,NpsiJ-1
              dVlVec(l) = vptemp(l)*( xsv(l+1) - xsv(l-1) )*0.5_R8
              iVlVec(l) = iVlVec(l-1) + dVlVec(l)
              MidVec(l)= 0.5_R8*(xsv(l+1)+xsv(l))
 6          continue
              l = NpsiJ
              dVlVec(l) = vptemp(l)* 0.5_R8*(xsv(l) - xsv(l-1))
              iVlVec(l) = iVlVec(l-1) + dVlVec(l)
              MidVec(l)= xsv(l)
            do 7 l=1,NpsiJ
              EdcVec(l) = voltlp(l)/(TWOPI * Rmaj)
 7          continue
            if(iEdc .eq. TRUE)then
               do 77 l = 1, NpsiJ
                  EdcVec(l) = Edcinp
 77            continue
            endif
 
      read(nTSCwrit,2203,end=98,err=99) psep(1),xsep(1),zsep(1)
      read(nTSCwrit,2203,end=98,err=99) psep(2),xsep(2),zsep(2)
      if( isym .eq. 1) then
      read(nTSCwrit,2203,end=98,err=99) ((PsiGrd(i, j) ,i=1,NumR),       &  
     &     j=NumZh,NumZ)
      else if (isym .eq. 0) then
      read(nTSCwrit,2203,end=98,err=99) ((PsiGrd(i, j) ,i=1,NumR),       &  
     &     j=    1,NumZ)
      else
        call LSCstop (' isym must be 1 for up/dn symmetry, or 0')
        return
      endif
      goto 100
!
!     Die gracefully on Read Error Here ............
!                                                  .
 98   call LSCstop(' EOF reached in RdTSC ')
      return
 99   call LSCstop(' read error  in RdTSC ')
      return
!                                                  .
!     Die gracefully on Read Error Here ............
!
 100  continue
!
!
 2201 format(10a8)
 2202 format(3i5,e16.6)
 2203 format(5e16.6)
 2204 format(4i5)
!
!
!           enmax = anecc(1)
!           Te0   = tekev(1)
!           Ti0   = tikev(1,1)
!
            do 1100 l=1,NpsiJ
            pi2Vec(l) = 0._R8
            AioVec(l) = 0._R8
            pe2Vec(l) = pe2Fac * (anecc(l)  /1.E14_R8)
            AelVec(l) = AelFac * (anecc(l)  /1.E14_R8) * tekev(l)
!
            do 1100 n=1,nspc
            TBTemp = ((amass(n)/PROTMS)*1.E24_R8)
            pi2Vec(l) = pi2Vec(l)                                        &  
     &                + pi2Fac * (anicc(l,n)/1.E14_R8) * achrg(n)**2     &  
     &                / TBTemp
!    ^                / (amass(n)/PROTMS*1.e24)  ! before Feb94
            AioVec(l) = AioVec(l)                                        &  
     &                + AioFac * (anicc(l,n)/1.E14_R8) * achrg(n)**2     &  
     &                / (TBTemp*TBTemp)             * tikev(l,n)
!    ^                / (amass(n)/PROTMS*1.e24)**2  * tikev(l,n) ! before Feb94
 1100       continue
!
            wcei2   =  0._R8
            do 1105 n=1,nspc
            TBTemp = ((amass(n)/PROTMS)*1.E24_R8)
            wcei2      = wcei2                                           &  
     &                 + ceiFac * anicc(1,n)/anecc(1) * achrg(n)**2      &  
     &                 / TBTemp
!    ^                 / (amass(n)/PROTMS*1.e24) ! before Feb94
 1105       continue
            wcei2   = wcei2 * Btesl**2
            capo2   = wcei2 / fghz**2
!
!
!
            pe2min = pe2Vec(NpsiJ - 1)
!
!                                       Make sure the update is at BOF
!
      rewind (unit = nTSCwrit)
!
!
      call plsgrd
!                                       plsgrd makes symmetry up-down
 
!     .                                 lcfs: last closed flux surface;
!     .                                 however, it is the 97% flux surface
      PsiVal = PsiLim                                                    &  
     &  -0.03_R8*(PsiLim-PsiMin)
!
!     .                                 Rmag not Rmaj
!     .                                 The 0.95 stays away avoids confusion
!     .                                 at the exact center where there can be
!     .                                 a tiny dip in psi (numerical noise).
      Rlcfs(1) = Rmag * 0.95_R8
      Zlcfs(1) = zmag
      call PsiSur ( PsiVal , Rlcfs , Zlcfs , ISIZElcfs , Nlcfs )
      Rlcfs(Nlcfs) = Rlcfs(1)
      Zlcfs(Nlcfs) = Zlcfs(1)
          RlcfsMax=Rmaj
          RlcfsMin=Rmaj
          ZlcfsMin=Zlcfs(1)
          ZlcfsMax=Zlcfs(1)
      do 1200 i=1,Nlcfs
          if(Rlcfs(i).gt.RlcfsMax)RlcfsMax=Rlcfs(i)
          if(Rlcfs(i).lt.RlcfsMin)RlcfsMin=Rlcfs(i)
          if(Zlcfs(i).gt.ZlcfsMax)ZlcfsMax=Zlcfs(i)
          if(Zlcfs(i).lt.ZlcfsMin)ZlcfsMin=Zlcfs(i)
 1200  continue
 
      RlcfsMax = RlcfsMax*1.01_R8
      RlcfsMin = RlcfsMin*0.99_R8
      ZlcfsMax = ZlcfsMax*1.01_R8
      ZlcfsMin = ZlcfsMin*1.01_R8
!     .                                 Keep the lcfs inside the grid!
      RlcfsMax = min(RlcfsMax, Rmax)
      RlcfsMin = max(RlcfsMin, Rmin)
      ZlcfsMax = min(ZlcfsMax, Zmax)
      ZlcfsMin = max(ZlcfsMin, Zmin)
 
      return
      ENTRY PutVinLSC(voltPut,npsiPut)
            do 2007 l=1,NpsiJ
              EdcVec(l) = voltPut(l+1)/(TWOPI * Rmaj)
 2007       continue
      return
 
      END
!
 
!                                                                      |
!                                                                      |
!     RdTSC.f(or)  ends                     ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
