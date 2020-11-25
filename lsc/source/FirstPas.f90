!#include "f77_dcomplx.h"
!     -----------------------------------------------------------------
      SUBROUTINE FirstPas
      USE params
      USE RayBins
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER NmPts, ThisPlot, TotlPlts, PltsOnPg, PltIndx
      INTEGER ir, iz, ip, iv, i
      PARAMETER (PltsOnPg =6)
      INTEGER xlow(PltsOnPg) , xhigh(PltsOnPg)
      INTEGER ylow(PltsOnPg) , yhigh(PltsOnPg)
      INTEGER ydrop(PltsOnPg), ywrit(PltsOnPg), ytitl
 
      INTEGER ixlo(PltsOnPg) , ixhig(PltsOnPg)
      INTEGER iylo(PltsOnPg) , iyhig(PltsOnPg)
      INTEGER iydro(PltsOnPg), iywri(PltsOnPg), iytit
 
      CHARACTER*20 MyString
      DATA    (xlow(i) , i=1,PltsOnPg)  /  50,400, 750, 50,400, 750 /
      DATA    (xhigh(i), i=1,PltsOnPg)  / 300,650,1000,300,650,1000 /
      DATA    (ylow(i) , i=1,PltsOnPg) /  400,400, 400, 50, 50,  50 /
      DATA    (yhigh(i), i=1,PltsOnPg) /  650,650, 650,300,300, 300 /
      DATA    (ywrit(i), i=1,PltsOnPg) /  625,625, 625,275,275, 275 /
      DATA    (ydrop(i), i=1,PltsOnPg) /  600,600, 600,250,250, 250 /
      DATA     ytitl                   /  700 /
 
      DATA    (ixlo(i) , i=1,PltsOnPg)  /   0,340, 680,  0,340, 680 /
      DATA    (ixhig(i), i=1,PltsOnPg)  / 340,680,1020,340,680,1020 /
      DATA    (iylo(i) , i=1,PltsOnPg) /  380,380, 380,  0,  0,   0 /
      DATA    (iyhig(i), i=1,PltsOnPg) /  680,680, 680,300,300, 300 /
      DATA    (iywri(i), i=1,PltsOnPg) /  700,700, 700,320,320, 320 /
      DATA    (iydro(i), i=1,PltsOnPg) /  675,675, 675,295,295, 295 /
      DATA     iytit                   /  750 /
 
 
      REAL*8     RtPsi(NZONDIM), RelPower(NZONDIM)
      INTEGER  iPin1, iPout1, iPin2, iPout2, iPin3, iPout3, iPctr
      INTEGER  iZin1, iZout1, iZin2, iZout2, iZin3, iZout3
      REAL*8      rpsi
      REAL*8      Pin1,  Pout1,  Pin2,  Pout2, Pin3, Pout3, Pmaxi
      DATA      Pin1,  Pout1,  Pin2,  Pout2, Pin3, Pout3, Pmaxi /        &  
     &          1.15_R8,   1.12_R8,  1.09_R8,   1.06_R8, 1.03_R8,        &  
     &          1.00_R8,  1.5_R8/
!     iPctr     integer Psi counter  1: in 1st pass; 2: out 1st pass
!     iPctr     integer Psi counter  3: in 2nd pass; 4: out 2nd pass
!     iPctr     integer Psi counter  5: in 3rd pass; 6: out 3rd pass or later
 
      REAL*8                                                             &  
     &        ZERO, ONE, RE81, RE82
      DATA ZERO, ONE/                                                    &  
     &      0.0_R8, 1.0_R8/
      REAL*8 AREAL
 
      rpsi = AREAL ( npsi )
      TotlPlts = nrays
 
      do 20 ir=1,nrays
      ThisPlot = ir
 
      PltIndx = mod(ThisPlot-1, PltsOnPg) + 1
      if (ThisPlot .eq. 1 ) then
        call EZinit
        call GNinit(nTSCgraf)
      endif
!
 
      iPctr = 1
      iPin1 = npsi
      do 10 iz = 1, nzones
      NmPts = iz
      ip = izind(iz,ir)
      iv = ivind(iz,ir)
 
      RtPsi(iz) = sqrt ( AREAL(ip) / rpsi )
      RelPower(iz) = power(iz,ir)/power(1,ir)
      if ( RelPower(iz) .le. 0.01_R8.or. iv .eq. 0 ) go to 11
!
      if (iPctr .eq. 1 ) then
      if ( ip .le. iPin1 ) then
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pin1
        iPin1 = ip
        iZin1 = iz
        iPctr = 1
      else
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pout1
        iPout1= ip
        iZout1= iz
        iPctr = 2
      endif
 
      else if (iPctr .eq. 2 ) then
      if ( ip .ge. iPout1 ) then
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pout1
        iPout1 = ip
        iZout1 = iz
        iPctr  = 2
      else
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pin2
        iPin2  = ip
        iZin2  = iz
        iPctr = 3
      endif
 
      else if (iPctr .eq. 3 ) then
      if ( ip .le. iPin2 ) then
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pin2
        iPin2 = ip
        iZin2 = iz
        iPctr = 3
      else
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pout2
        iPout2= ip
        iZout2= iz
        iPctr = 4
      endif
 
      else if (iPctr .eq. 4 ) then
      if ( ip .ge. iPout2 ) then
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pout2
        iPout2 = ip
        iZout2 = iz
        iPctr  = 4
      else
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pin3
        iPin3  = ip
        iZin3  = iz
        iPctr = 5
      endif
 
      else if (iPctr .eq. 5 ) then
      if ( ip .le. iPin3 ) then
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pin3
        iPin3 = ip
        iZin3 = iz
        iPctr = 5
      else
        if(RelPower(iz) .ge. 0.99_R8) RelPower(iz) = Pout3
        iPout3= ip
        iZout3= iz
        iPctr = 6
      endif
 
      endif
 
 10   continue
 11   continue
 
      call EZsets ( xlow(PltIndx), xhigh(PltIndx),                       &  
     &              ylow(PltIndx), yhigh(PltIndx),                       &  
     &             ZERO, ONE , ZERO,  Pmaxi, 1)
!    ^              0.0, 1.0 ,  0.0,  Pmaxi, 1)
 
      call GNsets ( ixlo(PltIndx), ixhig(PltIndx),                       &  
     &              iylo(PltIndx), iyhig(PltIndx),                       &  
     &             ZERO, ONE , ZERO,  Pmaxi, 1)
 
      if ( PltIndx .eq. PltsOnPg .or. ThisPlot .eq. TotlPlts ) then
        call EZwrit(xlow(1),ytitl,                                       &  
     &    ' Q-L absorption vs root psi for 6 rays$',0,0)
        call GNwrit(ixlo(1),iytit,                                       &  
     &    ' Q-L absorption vs root psi for 6 rays$',0,0)
 
        call MkGrfLst(' Q-L absrpn vs root psi 6 rays')
      endif
 
!     call EZaxes ( 1,10,1,15 )
      call EZaxes ( 1, 5,3,1)
 
      write(MyString,'('' Ray '',i2,''$'')') ThisPlot
      call EZwrit(xlow(PltIndx),ywrit(PltIndx),  MyString , 0,0 )
      call GNwrit(ixlo(PltIndx),iywri(PltIndx),  MyString , 0,0 )
      write(MyString,'('' Nll '',f5.2,''$'')') npar(1,ThisPlot)
      call EZwrit(xlow(PltIndx),ydrop(PltIndx),  MyString , 0,0 )
      call GNwrit(ixlo(PltIndx),iydro(PltIndx),  MyString , 0,0 )
 
      call EZpnts (RtPsi,RelPower,NmPts)
      call GNpnts (RtPsi,RelPower,NmPts)
!
!
      if ( PltIndx .eq. PltsOnPg .or. ThisPlot .eq. TotlPlts ) then
        call EZwrit(xlow(1),ytitl,                                       &  
     &    ' Q-L absorption vs root psi for 6 rays$',0,0)
        call GNwrit(ixlo(1),iytit,                                       &  
     &    ' Q-L absorption vs root psi for 6 rays$',0,0)
 
        call MkGrfLst(' Q-L absrpn vs root psi 6 rays')
 
        call EZfini(0,0)
        call GNfini
      endif
 
 20   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
