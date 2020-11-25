!#include "f77_dcomplx.h"
!
!     ------------------------------------------------------------------
!
      SUBROUTINE JrfDiagn
      USE Doflags
      USE FeBins
      USE Jrf
      USE params
      USE PlPr
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE tscunits
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ivlocal, iplocal
      INTEGER ipstep, ThisPlot, NumPlts, iCross0
      CHARACTER*1 curve, dots, altdots, cros, bars
      CHARACTER*10 chip
      CHARACTER*24 MyString
      CHARACTER*30 TSCstring
      REAL*8    nparlocal
      REAL*8 AREAL
 
      INTEGER NumJofN, inlocal, iymjr, iymnr, nsmo
      REAL*8    NpLocMin, NpLocMax
      PARAMETER (NumJofN=101, NpLocMin=-10._R8, NpLocMax=10._R8)
      REAL*8    JofNpar(NumJofN), NparVec(NumJofN)
      REAL*8    JofNmin, JofNmax, NparLog(20)
      PARAMETER (nsmo=5)
      REAL*8    smofn(nsmo)
      REAL*8                                                             &  
     &        ZERO, ONE, ONEPT2
      DATA    ZERO, ONE, ONEPT2/                                         &  
     &         0.0_R8, 1.0_R8,    1.2_R8/
      DATA curve, dots,altdots, cros, bars,  NumPlts  /                  &  
     &       ' ',  '.',    ',',  'x',  '|',       12  /
      DATA smofn(1), smofn(2), smofn(3), smofn(4), smofn(5) /            &  
     &      1._R8,  2._R8,  3._R8,  2._R8,  1._R8/
      ipstep = (npsi - 2)/NumPlts
      ivlocal = 1
 
      if(PlFlg(RFDPL) .ge. TRUE)then
            ivlocal = 1
            ThisPlot = 1
            do 20 iplocal = ipstep, ipstep*NumPlts, ipstep
!     .                                 if chip starts with a 'l' or 'L'
!     .                                 then the logarithm is plotted
               write(chip,'(''Jray('',i2,'')'')') iplocal
!
!     .                                 call the curve line first... this
!     .                                 sets the scale of the graph
               call Dql6Norm(Vpar(ivlocal),Jray(ivlocal,iplocal),nv,     &  
     &                           chip,bars, ThisPlot, NumPlts)
!
               ThisPlot = ThisPlot + 1
!
 20         continue
!
            call BrodCast(nv, wkv,    0._R8)
            do ivlocal = 1,nv
               do iplocal = 1,npsi
                  wkv(ivlocal)=wkv(ivlocal) + Jray(ivlocal,iplocal)
               enddo
            enddo
!           EZnorm ( fin , fout , fmin , fmax , n , iCross0)
!            iCross0=1+1=2 --> crosses 0; iCross0=1+0=1 --> not cross 0
            call EZnorm ( wkv, wkv , wkv(nv+1),wkv(nv+2), nv, iCross0 )
! bar graph of Pray(iv)
            call EZsets (    600,    1000, 125, 375,                     &  
     &                       -ONE,    ONE,-ONEPT2, ONEPT2, 1)
            call GNsets (    500,    1000,   0, 380,                     &  
     &                       -ONE,    ONE,-ONEPT2, ONEPT2, 1)
            write (MyString,'( ''Jray(v/c)/'', 1pe8.1,''$'')') wkv(nv+2)
            call EZwrit(600, 400, MyString , 0, 0 )
            call GNtitl(          MyString        )
            call EZaxes (2,5,1,6)
            call EZbars(vpar,wkv,nv,'y')
            call GNbars(vpar,wkv,nv,'y')
!
! continuous graph of Jray(iv)
            call EZsets (    600,    1000, 475, 725,                     &  
     &                       -ONE,    ONE,-ONEPT2, ONEPT2, 1)
            call GNsets (    500,    1000, 380, 760,                     &  
     &                       -ONE,    ONE,-ONEPT2, ONEPT2, 1)
            write (MyString,'( ''Jray(v/c) curve'',''$'')')
            call EZwrit(600, 750, MyString , 0, 0 )
            call GNtitl(          MyString        )
            call EZaxes (2,5,1,6)
            call EZcurv(vpar,wkv,nv)
            call GNcurv(vpar,wkv,nv)
! bar graph of Jray(npar)
            call BrodCast(NumJofN, JofNpar, 0.0_R8)
            call ugrid(NparVec,NumJofN, NpLocMin, NpLocMax)
            do ivlocal = 1,ivZero-1
               nparlocal= max(1._R8/vpar(ivlocal), NpLocMin)
               inlocal = int((nparlocal-NpLocMin)                        &  
     &                         / 20._R8*AREAL(NumJofN-1)) + 1
               JofNpar(inlocal) = JofNpar(inlocal) + wkv(ivlocal)
            enddo
            do ivlocal = ivZero+1,nv
               nparlocal= min(1._R8/vpar(ivlocal), NpLocMax)
               inlocal = int((nparlocal-NpLocMin)                        &  
     &                         / 20._R8*AREAL(NumJofN-1)) + 1
               JofNpar(inlocal) = JofNpar(inlocal) + wkv(ivlocal)
!               call EZbars(nparlocal,wkv(ivlocal),1,'y')
            enddo
            call vecnorm(NumJofN, JofNpar)
            call vecMnMx(JofNpar,NumJofN, JofNmin, JofNmax)
            call EZrnd2  (JofNmax, JofNmax, iymjr, iymnr)
            call EZsets (  100,     500, 125, 375,                       &  
     &                NpLocMin,NpLocMax,-JofNmax, JofNmax, 1)
            call GNsets (    0,     500,   0, 380,                       &  
     &                NpLocMin,NpLocMax,-JofNmax, JofNmax, 1)
            write (MyString,'(''Jray(npar)'',''$'')')
            call EZwrit(100, 400, MyString , 0, 0 )
            call GNtitl(          MyString        )
            call EZaxes (2,5,iymjr,iymnr)
            call EZbars(NparVec,JofNpar,NumJofN,'y')
            call GNbars(NparVec,JofNpar,NumJofN,'y')
! continuous graph of Pray(npar)
            call vecnorm(nsmo,smofn)
            call smooth(JofNpar,NumJofN,NumJofN,1,nsmo,smofn)
            call EZsets (  100,     500, 475, 725,                       &  
     &                NpLocMin,NpLocMax,-JofNmax, JofNmax, 1)
            call GNsets (    0,     500, 380, 760,                       &  
     &                NpLocMin,NpLocMax,-JofNmax, JofNmax, 1)
            write (MyString,'(''Jray(npar) curve'',''$'')')
            call EZwrit(100, 750, MyString , 0, 0 )
            call GNtitl(          MyString        )
            call EZaxes (2,5,iymjr,iymnr)
            call EZcurv(NparVec,JofNpar,NumJofN)
            call GNcurv(NparVec,JofNpar,NumJofN)
            write(TSCstring,'( '' Jray(iv) vs v/c '')')
! close frame
            call MkGrfLst( TSCstring )
            call EZfini(0,0)
            call GNfini
         endif
      return
      END
!     .                                                                |
!     .                                                                |
!     jrf.F ends   -----------------------------------------------------
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
