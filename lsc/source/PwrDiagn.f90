!#include "f77_dcomplx.h"
!
!
      SUBROUTINE PwrDiagn
      USE Doflags
      USE FeBins
      USE params
      USE PlPr
      USE power_mod
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE tscunits
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nrLines
      INTEGER ivlocal, iplocal
      INTEGER ipstep, ThisPlot, NumPlts, iCross0
      CHARACTER*1 curve, dots, altdots, cros, bars
      CHARACTER*10 chip
      CHARACTER*24 MyString
      CHARACTER*30 TSCstring
      REAL*8    nparlocal
      REAL*8 AREAL
 
      INTEGER NumPofN, inlocal, iymjr, iymnr, nsmo
      REAL*8    NpLocMin, NpLocMax
      PARAMETER (NumPofN=101, NpLocMin=-10._R8, NpLocMax=10._R8)
      REAL*8    PofNpar(NumPofN), NparVec(NumPofN)
      REAL*8    PofNmin, PofNmax, NparLog(20), NparINV(NVELDIM)
      PARAMETER (nsmo=5)
      REAL*8    smofn(nsmo)
 
      DATA curve, dots,altdots, cros, bars,  NumPlts  /                  &  
     &       ' ',  '.',    ',',  'x',  '|',       12  /
      DATA smofn(1), smofn(2), smofn(3), smofn(4), smofn(5) /            &  
     &      1._R8,  2._R8,  3._R8,  2._R8,  1._R8/
      REAL*8    ZERO, ONE, ONEPT2
      DATA    ZERO, ONE, ONEPT2/                                         &  
     &         0.0_R8, 1.0_R8,    1.2_R8/
      ipstep = (npsi - 2)/NumPlts
      ivlocal = 1
 
 
      nrLines = 0
      if(PlFlg(RFDPL)  .ge. TRUE ) then
!                                       Compute power along rays
         call RfDamp
         if(PlFlg(RFDPL) .ge. TRUE) then
!           call mwpl('PRay$', PRay, DqlWind, GEOM)
            ivlocal = 1
            ThisPlot = 1
            do 20 iplocal = ipstep, ipstep*NumPlts, ipstep
!     .                                 if chip starts with a 'l' or 'L'
!     .                                 then the logarithm is plotted
               write(chip,'(''Pray('',i2,'')'')') iplocal
!
!     .                                 call the curve line first... this
!     .                                 sets the scale of the graph
               call Dql6Norm(Vpar(ivlocal),Pray(ivlocal,iplocal),nv,     &  
     &                           chip,bars, ThisPlot, NumPlts)
!
               ThisPlot = ThisPlot + 1
!
 20         continue
!
            call BrodCast(nv, wkv,    0._R8)
            do ivlocal = 1,nv
                  NparINV(ivZero)= 0.0_R8
               if (ivlocal .gt. ivZero )                                 &  
     &            NparINV(ivlocal) = 11._R8- 10._R8*vpar(ivlocal)
               if (ivlocal .lt. ivZero )                                 &  
     &            NparINV(ivlocal) =-11._R8- 10._R8*vpar(ivlocal)
               do iplocal = 1,npsi
                  wkv(ivlocal)=wkv(ivlocal) + Pray(ivlocal,iplocal)
               enddo
            enddo
!           EZnorm ( fin , fout , fmin , fmax , n , iCross0)
!            iCross0=1+1=2 --> crosses 0; iCross0=1+0=1 --> not cross 0
            call EZnorm ( wkv, wkv , wkv(nv+1),wkv(nv+2), nv, iCross0 )
! bar graph of Pray(iv)
            call EZsets (    600,    1000, 125, 375,                     &  
     &                       -ONE,    ONE,ZERO, ONEPT2, 1)
!    ^                       -1.,      1., 0.0, 1.2, 1)
            call GNsets (    500,    1000,   0, 380,                     &  
     &                       -ONE,    ONE,ZERO, ONEPT2, 1)
            write (MyString,'( ''Pray(v/c)/'', 1pe8.1,''$'')') wkv(nv+2)
            call EZwrit(600, 400, MyString , 0, 0 )
            call GNtitl(          MyString        )
            call EZaxes (2,5,1,6)
            call EZbars(vpar,wkv,nv,'y')
            call GNbars(vpar,wkv,nv,'y')
!
! continuous graph of Pray(iv)
            call EZsets (    600,    1000, 475, 725,                     &  
     &                       -ONE,    ONE,ZERO, ONEPT2, 1)
!    ^                       -1.,      1., 0.0, 1.2, 1)
            call GNsets (    500,    1000, 380, 760,                     &  
     &                       -ONE,    ONE,ZERO, ONEPT2, 1)
            write (MyString,'( ''Pray(v/c) curve'',''$'')')
            call EZwrit(600, 750, MyString , 0, 0 )
            call GNtitl(          MyString        )
            call EZaxes (2,5,1,6)
            call EZcurv(vpar,wkv,nv)
            call GNcurv(vpar,wkv,nv)
! bar graph of Pray(npar)
            call BrodCast(NumPofN, PofNpar, 0.0_R8)
            call ugrid(NparVec,NumPofN, NpLocMin, NpLocMax)
            do ivlocal = 1,ivZero-1
               nparlocal= max(1._R8/vpar(ivlocal), NpLocMin)
               inlocal = int((nparlocal-NpLocMin)                        &  
     &                         / 20._R8*AREAL(NumPofN-1)) + 1
               PofNpar(inlocal) = PofNpar(inlocal) + wkv(ivlocal)
            enddo
            do ivlocal = ivZero+1,nv
               nparlocal= min(1._R8/vpar(ivlocal), NpLocMax)
               inlocal = int((nparlocal-NpLocMin)                        &  
     &                         / 20._R8*AREAL(NumPofN-1)) + 1
               PofNpar(inlocal) = PofNpar(inlocal) + wkv(ivlocal)
!               call EZbars(nparlocal,wkv(ivlocal),1,'y')
            enddo
            call vecnorm(NumPofN, PofNpar)
            call vecMnMx(PofNpar,NumPofN, PofNmin, PofNmax)
            call EZrnd2  (PofNmax, PofNmax, iymjr, iymnr)
            call EZsets (  100,     500, 125, 375,                       &  
     &                NpLocMin,NpLocMax,ZERO, PofNmax, 1)
!    ^                NpLocMin,NpLocMax, 0.0, PofNmax, 1)
            call GNsets (    0,     500,   0, 380,                       &  
     &                NpLocMin,NpLocMax,ZERO, PofNmax, 1)
            write (MyString,'(''Pray(npar)'',''$'')')
            call EZwrit(100, 400, MyString , 0, 0 )
            call GNtitl(          MyString        )
            call EZaxes (2,5,iymjr,iymnr)
            call EZbars(NparVec,PofNpar,NumPofN,'y')
            call GNbars(NparVec,PofNpar,NumPofN,'y')
! continuous graph of Pray(npar)
            call vecnorm(nsmo,smofn)
            call smooth(PofNpar,NumPofN,NumPofN,1,nsmo,smofn)
            call EZsets (  100,     500, 475, 725,                       &  
     &                NpLocMin,NpLocMax,ZERO, PofNmax, 1)
!    ^                NpLocMin,NpLocMax, 0.0, PofNmax, 1)
            call GNsets (    0,     500, 380, 760,                       &  
     &                NpLocMin,NpLocMax,ZERO, PofNmax, 1)
            write (MyString,'(''Pray(npar) curve'',''$'')')
            call EZwrit(100, 750, MyString , 0, 0 )
            call GNtitl(          MyString        )
            call EZaxes (2,5,iymjr,iymnr)
            call EZcurv(NparVec,PofNpar,NumPofN)
            call GNcurv(NparVec,PofNpar,NumPofN)
            write(TSCstring,'( '' Pray(iv) vs v/c '')')
! close frame
            call MkGrfLst( TSCstring )
            call EZfini(0,0)
            call GNfini
         endif
 
      endif
!
      if (PlFlg(RFDPL)   .ge. TRUE    ) then
!                                       Compute power from QL heating
         call RfHeat
         call JdepCalc
         call PJrfIgrl
 
      endif
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
