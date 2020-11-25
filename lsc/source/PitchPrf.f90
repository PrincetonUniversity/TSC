!#include "f77_dcomplx.h"
!
!     -----------------------------------------------------------------
!
      SUBROUTINE PitchPrf
      USE dielec
      USE params
      USE power_mod
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,ii,ix,ixi,iy,iyi,j, NUMPTS, NUMSMO
      PARAMETER (NUMPTS=200, NUMSMO=8)
      REAL*8                                                             &  
     &        r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      REAL*8                                                             &  
     &        x(NUMPTS),Pitch(NUMPTS),PitchDeg(NUMPTS),PowerDep(NUMPTS)
      REAL*8                                                             &  
     &        Bpol, Btot, Enh(NUMPTS)
      REAL*8                                                             &  
     &        Dens(NUMPTS), TempAry(NUMPTS)
      REAL*8                                                             &  
     &        PowerMax, PitchMin, PitchMax, DensMax,                     &  
     &        temp1,    temp2,    temp3
 
      REAL*8                                                             &  
     &        ZERO, ONE, ONEPT5, TEN
      DATA    ZERO, ONE, ONEPT5, TEN/                                    &  
     &         0.0_R8, 1.0_R8,    1.5_R8,10.0_R8/
      REAL*8 AREAL
 
!
      do i=1,NUMPTS
        x(i) = RlcfsMin +                                                &  
     &        (RlcfsMax-RlcfsMin)/AREAL(NUMPTS-1)*AREAL(i-1)
        r = x(i)
        z = 0.0_R8
        call plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
!
        if (pe2 .le. 0.0_R8) then
           Enh(i) = 1.0_R8
           Bpol   = 0.0_R8
        else
            call eps (r , z , 0._R8, 0._R8)
            Bpol = sqrt ( Br*Br + Bz*Bz )
            Btot = sqrt ( Bpol*Bpol + RBphi*RBphi/r /r  )
            Enh(i) = ( r /RlcfsMax ) * ( 1._R8                           &  
     &                       -sqrt(abs(1._R8-Epar))*                     &  
     &                       Bpol/Btot/sqrt(abs(Eper))  )                   
            if(Enh(i) .le. 0.1_R8) Enh(i) = 0.1_R8
            Enh(i) = 1.00_R8/ Enh(i)
        endif
!
        Pitch(i)    = Bz/RBphi*r
        PitchDeg(i) = atan(Pitch(i))* (180._R8/3.1415926_R8)
        ii = 0
        PowerDep(i) = 0.0_R8
        Dens(i)     = pe2
        if (pe2 .gt. 0.00_R8) then
          do j=2,Npsi
            if(PsiAry(j) .ge. psi .and. PsiAry(j-1) .lt. psi ) then
              PowerDep(i) = PowerDep(i) + PrayTot(j)
              ii = ii+1
            endif
          enddo
              PowerDep(i) = PowerDep(i)/(AREAL(ii)+.001_R8)
          continue
        endif
 
      enddo
!     .
      do 21 i=1,NUMPTS
        TempAry(i)=PowerDep(i)
 21   continue
!
      do 30 i=1,NUMPTS
        temp1 = 0.00_R8
        do 25 j=-NUMSMO,NUMSMO,1
          ii = i+j
          if (ii .lt. 1     ) ii = 1
          if (ii .gt. NUMPTS) ii = NUMPTS
          temp1 = temp1+TempAry(ii)
 25     continue
        PowerDep(i) = temp1 / AREAL(2*NUMSMO + 1)
 30   continue
 
      PowerMax = 0.00_R8
      PitchMin = 0.00_R8
      PitchMax = 0.00_R8
      do 40 i=1,NUMPTS
        if (PitchMin .gt.    Pitch(i)) PitchMin = Pitch(i)
        if (PitchMax .lt.    Pitch(i)) PitchMax = Pitch(i)
        if (DensMax  .lt.    Dens(i) ) DensMax  = Dens(i)
        if (PowerMax .lt. PowerDep(i)) PowerMax = PowerDep(i)
 40   continue
 
      do 50 i=1,NUMPTS
        PowerDep(i)=PowerDep(i)/(PowerMax + 1.0E-20_R8)*0.98_R8
        Dens(i)    =Dens(i)    /(DensMax  + 1.0E-20_R8)*0.98_R8
 50   continue
 
      PitchMax = max(PitchMax, abs(PitchMin))
 
      call EZrnd2 (RlcfsMax,temp2,i,ii)
      temp1 = temp2-RlcfsMin
      call EZrnd2 (temp1,temp1,ix,ixi)
      temp1=temp2-temp1
!      temp2 = 2.0
!      temp1 = 1.2
      ix =4
      ixi=1
      call EZinit
      call EZrndu(PitchMax,PitchMax, iy)
      iyi = 1
      PitchMin = - PitchMax
      call EZsets (100,450,475,700,temp1,temp2,PitchMin,PitchMax,1)
      call EZwrit (125, 725,'Bz/Bphi vs Rmaj$',0,0)
 
      call GNsets (  0,500,380,760,temp1,temp2,PitchMin,PitchMax,1)
      call GNtitl (         'Bz/Bphi vs Rmaj$'    )
 
      call EZcurv ( x  , Pitch, NUMPTS)
      call EZaxes (ix,ixi,iy,iyi)
      call GNcurv ( x  , Pitch, NUMPTS)
!
      PitchMax = 10._R8
      PitchMin =-10._R8
      iy=2
      iyi=2
      call EZsets (100,450,100,325,temp1,temp2,PitchMin,PitchMax,1)
      call EZwrit (125, 350,'Pitch in deg vs Rmaj$',0,0)
      call EZcurv ( x  , PitchDeg, NUMPTS)
      call EZaxes (ix,ixi,iy,iyi)
 
      call GNsets (  0,500,  0,380,temp1,temp2,PitchMin,PitchMax,1)
      call GNtitl (         'Pitch in deg vs Rmaj$'    )
      call GNcurv ( x  , PitchDeg, NUMPTS)
!
      call EZsets (600,950,100,325,temp1,temp2, ZERO   ,ONEPT5,1)
      call EZwrit (625,        350,'Ne_ Prf... vs Rmaj$',0,0)
      call EZaxes (ix,ixi,   3,1)
      call EZcurv ( x  , Dens    , NUMPTS )
      call EZpnts ( x  , PowerDep, NUMPTS )
 
      call GNsets (500,1000, 0,380,temp1,temp2, ZERO   ,ONEPT5,1)
      call GNtitl (                'Ne_ Prf... vs Rmaj$'    )
      call GNcurv ( x  , Dens    , NUMPTS )
      call GNpnts ( x  , PowerDep, NUMPTS )
 
      call EZsets (600,950,475,700,temp1,temp2, ZERO   ,TEN,1)
      call EZwrit (625,        725,'Npar enhance vs Rmaj$',0,0)
      call EZaxes (ix,ixi,   2,5)
      call EZpnts ( x  , Enh     , NUMPTS )
 
      call GNsets (500,1000,380,760,temp1,temp2, ZERO   ,TEN,1)
      call GNtitl (                'Npar enhance vs Rmaj$'    )
      call GNpnts ( x  , Enh     , NUMPTS )
!
!
      call MkGrfLst(' Pitch Pray Ne Npar gain vs R ')
      call EZfini(0,0)
      call GNfini
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
