!
!     ------------------------------------------------------------------
!
      SUBROUTINE JdepFromTSC(curtsc,djdets2,npsiTake)
      USE Jrf
      USE params
      USE PIetc
      USE ProfBody
      USE RayBins
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
      INTEGER l, ip, ips, iGotRuna, npsiTake, PassNumber, iFillJray
      REAL*8    curtsc(npsiTake), djdets2(npsiTake), curtscI(PPSI)
      REAL*8    xlookup, yreturn
      REAL*8    PowFac        , CurFac
      DATA    PowFac        , CurFac        /                            &  
     &            1.E-6_R8,     1.E-4_R8/
!             m^3 per cm^3  , m^2 per cm^2
      DATA iFillJray / 0 /
      REAL*8    ZERO
      DATA    ZERO/                                                      &  
     &         0.0_R8/
!
!
      Passnumber = 1
!
 1    continue
!
      if (PassNumber .eq. 1) then
        call GetEdc(ZERO)
!       call GetEdc(0.0)
      else
        call GetEdc(dEdcAmnt)
      endif
!                                       Fill array EdcAry
!     .                                 with 0.0 difference from the EdcVec
!     .                                 transferred from TSC; on second pass
!     .                                 fill with the dEdcAmnt
      do 10 ip = 1, npsi
        call jnorm(ip)
!     .                                 compute jrf
        call mkj(ip,js(ip), iGotRuna, iFillJray)
 
 10   continue
!
      IrIntgrl(1) = js (1) * dVol(1) / (TWOPI * xmag )
      do 17 ips = 2, npsi
         IrIntgrl(ips) = IrIntgrl(ips-1) + js(ips)*dVol(ips) /           &  
     &                                    (TWOPI * xmag )
!     Divide by 2 pi R to make an area element !
 17   continue
 
      do 20 l = 1, NpsiJ - 1
         xlookup =  MidVec(l)
!
!     interpolate for current integral
         call linr1d(Npsi, MidAry, IrIntgrl,                             &  
     &               xlookup, yreturn)
         curtscI (l)=  yreturn
 20   continue
 
      curtscI (NpsiJ) = IrIntgrl(Npsi)
      curtscI (NpsiJ-1) = IrIntgrl(Npsi)
 
      l = 1
      if(PassNumber .eq. 1)                                              &  
     &curtsc (l)= CurFac * (curtscI (l)             )/                   &  
     &                    ( dVlVec(l) / (TWOPI * xmag) )
      if(PassNumber .eq. 2)                                              &  
     &djdets2(l)= CurFac * (curtscI (l)             )/                   &  
     &                    ( dVlVec(l) / (TWOPI * xmag) )
      do 30 l = 2, NpsiJ - 1
        if(PassNumber .eq. 1)                                            &  
     &  curtsc (l)= CurFac * (curtscI (l)-curtscI (l-1))/                &  
     &                      ( dVlVec(l) / (TWOPI * xmag) )
!
!     .                                 djdets2 contains temporarily
!     .                                 the current at the incremented E
        if(PassNumber .eq. 2)                                            &  
     &  djdets2(l)= CurFac * (curtscI (l)-curtscI (l-1))/                &  
     &                      ( dVlVec(l) / (TWOPI * xmag) )
 
!       I have to mult by 2 pi R to make an area!
 
 30   continue
 
      if(PassNumber .eq. 1) then
        PassNumber = 2
        goto 1
      endif
 
!
!     !! Are NpsiJ and npsiTake right??
!
      do 40 l=npsiTake,2,-1
        curtsc (l)=curtsc (l-1)
        djdets2(l)=djdets2(l-1)
 40   continue
        l = 1
        curtsc (l)= 0.00_R8
        djdets2(l)=0.00_R8
 
      do 50 l=2,npsiTake
        djdets2(l) = (djdets2(l) - curtsc (l))/ dEdcAmnt
 50   continue
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
