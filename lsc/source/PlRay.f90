!#include "f77_dcomplx.h"
!     PlRay     begins                      ---------------------------|
!     -                                                                |
!     -                                                                |
 
      SUBROUTINE PlRay(iFillPsi)
      USE dielec
      USE Doflags
      USE FeBins
      USE params
      USE PIetc
      USE RayBins
      USE RayWrk
      USE TSCgrap
      USE tscunits
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL psisur, plasma2d
!
!     Screen size is 1024 points wide by 780 points high
!     Nominal location of graphs in the x or horizontal direction
!     -       is   50-300    400-650    750-1000
!     Nominal location of graphs in the y or vertical   direction
!     -       is   50-300    400-650
!
      CHARACTER*10 FileNa
      CHARACTER*30 MyString
      INTEGER i, igraph, inc, incin, incx, incxi,                        &  
     &        ipol, ipoli, j, NPSISUR, Npts, iFillPsi
      INTEGER iDoAcces
      REAL*8                                                             &  
     &        ALPHA, Bpol, cosa, csx, csy, Cval, dR, dS, dZ,             &  
     &        niceMx,                                                    &  
     &        parnx, pernx, PsiVal, Qval,                                &  
     &        rmaxx, rminn, sina,  woc, Xpol, Xtor, zmaxx
      REAL*8                                                             &  
     &                   r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      PARAMETER (NPSISUR = 3)
!
      DATA  FileNa      /                                                &  
     &     'surf.dat'   /
      DATA  iDoAcces /1/
      INTEGER ipsisur, npsipts(NPSISUR), icross(4), nfound
      REAL*8     rpsi(NPSISUR,NPLTDIM), zpsi(NPSISUR,NPLTDIM)
      REAL*8     RRary(NPLTDIM), ZZary(NPLTDIM), ep
 
      DATA   ALPHA /                                                     &  
     &    0.523599_R8/
!     ALPHA = PI/6 = 0.52 for 30 degrees
!     parnx = 10 --> very heavy electron damping
!     pernx = 150--> synchronism with 40 kev D ions
      DATA parnx,pernx,igraph,inc,incin,incx,incxi /                     &  
     &       10._R8, 150._R8,     0,  2,    2,   2,    5 /
      REAL*8                                                             &  
     &        ZERO, ONE, RE81, RE82, ONE80
      DATA ZERO, ONE, ONE80/                                             &  
     &      0.0_R8, 1.0_R8, 180._R8/
      REAL*8 AREAL

      if ( iplt .eq. NPLTDIM )                                           &  
     & write(nLSCcom2,                                                   &  
     &         '(''  , PLOT ARRAY IS FILLED '' )'                        &  
     &          )
 
!     -                                 Fill Arrays of The Psi Surfaces  IF
      if (iFillPsi .eq. 1) then
      PsiVal = PsiLim                                                    &  
     &  -0.05_R8*(PsiLim-PsiMin)
!    ^  -0.03*(PsiLim-PsiMin) ! old reliable
!     .                                 Rmag not Rmaj
!     .                                 The 0.95 stays away avoids confusion
!     .                                 at the exact center where there can be
!     .                                 a tiny dip in psi (numerical noise).
      RRary(1) = Rmag * 0.95_R8
      ZZary(1) = 0._R8
      call PsiSur ( PsiVal , RRary , ZZary , NPLTDIM , Npts )
      RRary(Npts) = RRary(1)
      ZZary(Npts) = ZZary(1)
          rmaxx=Rmaj
          rminn=Rmaj
          zmaxx=0._R8
          ipsisur=1
          npsipts(ipsisur)=Npts
      do 5 i=1,Npts
          if(rrary(i).gt.rmaxx)rmaxx=rrary(i)
          if(rrary(i).lt.rminn)rminn=rrary(i)
          if(abs(zzary(i)).gt.zmaxx)zmaxx=abs(zzary(i))
          rpsi(ipsisur,i) = RRary(I) - Rmaj
          zpsi(ipsisur,i) = ZZary(I)
5     continue
 
      PsiVal = PsiLim                                                    &  
     &  -0.30_R8*(PsiLim-PsiMin)
!                                       Rmag not Rmaj
      RRary(1) = Rmag * 0.95_R8
      ZZary(1) = 0._R8
      call PsiSur ( PsiVal , RRary , ZZary , NPLTDIM , Npts )
      RRary(Npts) = RRary(1)
      ZZary(Npts) = ZZary(1)
          ipsisur=2
          npsipts(ipsisur)=Npts
      DO 6 I=1,Npts
          rpsi(ipsisur,i) = RRary(I) - Rmaj
          zpsi(ipsisur,i) = ZZary(I)
 6    continue
 
      PsiVal = PsiLim                                                    &  
     &  -0.70_R8*(PsiLim-PsiMin)
!                                       Rmag not Rmaj
!     RRary(1) = Rmag * 0.90            0.90 might be needed for GA D-III-D
      RRary(1) = Rmag * 0.90_R8
      ZZary(1) = 0._R8
 13   call     plasma2d (RRary(1),ZZary(1), psi,                         &  
     &                   Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      if (psi .gt. PsiVal .and.                                          &  
     &    RRary(1) .eq. 0.90_R8*Rmag ) then
        RRary(1) =  0.5_R8*(RRary(1) + Rmag)
        goto 13
      endif
 
      call PsiSur ( PsiVal , RRary , ZZary , NPLTDIM , Npts )
           if (iError .gt. 0)  return
 
      RRary(Npts) = RRary(1)
      ZZary(Npts) = ZZary(1)
          ipsisur=3
          npsipts(ipsisur)=Npts
      DO 7 I=1,Npts
          rpsi(ipsisur,i) = RRary(I) - Rmaj
          zpsi(ipsisur,i) = ZZary(I)
 7    continue
 
!     q(Psi)    = RBphi(Psi)/TWOPI * Integral dSpol/(R**2 Bpol)
!     Curnt(MA) =   10/FourPi * Integral dSpol * Bpol
      Qval      = 0._R8
      Cval      = 0._R8
      ipsisur   = 1
      Npts      = npsipts(ipsisur)
      Do 8    I = 2, Npts
      dR  = rpsi(ipsisur,i) - rpsi(ipsisur,(i-1))
      dZ  = zpsi(ipsisur,i) - zpsi(ipsisur,(i-1))
      dS  = sqrt ( dR**2 + dZ**2 )
       R  = Rmaj + ( rpsi(ipsisur,i) + rpsi(ipsisur,i-1) )/2._R8
       Z  =        ( zpsi(ipsisur,i) + zpsi(ipsisur,i-1) )/2._R8
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      Bpol = SQRT ( Br**2 + Bz**2 )
      Qval = Qval + dS/(Bpol*R**2 )
      Cval = Cval + dS* Bpol
8     Continue
      Qval = Qval * RBphi  /    TWOPI
      Cval = Cval * 10/(2._R8*TWOPI)
      qlim = Qval
      q0   = 1._R8
!     -                                 finished filling surface arrays
 
 
!     -                                 write the file now
!            open (nTSCunus,status='unknown',file = FileNa)
!            do 6001 i=1,NPSISUR
!               write(nTSCunus,5600) npsipts(i)
!               do 6000 j=1,npsipts(i)
!                  write(nTSCunus,5601) rpsi(i,j), zpsi(i,j)
! 6000          continue
! 6001       continue
! 5600   format(i3)
! 5601   format(2e15.8)
! 5602   format(3e15.8)
!
!            write(nTSCunus,5602) rminn, rmaxx, zmaxx
!            write(nTSCunus,5602) qlim,  q0   , Cval
!            close(nTSCunus)
      endif
!     .                                 ENDIF on finding psi surfaces ENDIF
!
!
!    -                                  The psi surfaces used in plotting
!    -                                  are now defined.
!
!    -                                  First page of ray plots:
!    -                                  Isometric of ray in the torus:
!    -
      call EZinit
      xtor = zmaxx + rmaxx * R1O2
      if ( rmaxx*R3O2 .gt. xtor ) xtor = rmaxx*R3O2
      call EZsets ( 620, 970,   35,385 ,                                 &  
     &            -xtor,xtor,-xtor,xtor, 1 )
!     call GNsets ( 620, 970,   35,385 ,
!    ^            -xtor,xtor,-xtor,xtor, 1 )
      do 20 j=1,12
      csx = cos (TWOPI/12._R8*AREAL(j-1) + PIO4 )
      csy = cos (TWOPI/12._R8*AREAL(j-1) - PIO4 )
      ipsisur=1
      Npts   =npsipts(ipsisur)
      do 10  i=1,Npts
      xp(i) = r3o2*(rmaj+rpsi(ipsisur,i))*csx                              
10    yp(i) =-r1o2*(rmaj+rpsi(ipsisur,i))*csy                            &
     &                  +zpsi(ipsisur,i)
      call ezcurv (xp,yp,npts)
!     call GNcurv (xp,yp,npts)
20    continue
      cosa = cos(alpha)
      sina = sin(alpha)
      do 50   i=1,iplt
      xp(i) =          cosa*yrr(i)*(cos(ypp(i)) - sin(ypp(i)))
50    yp(i) = yzz(i) - sina*yrr(i)*(cos(ypp(i)) + sin(ypp(i)))
      call EZcurv (xp,yp,iplt)
      write(MyString, '('' Ray '',i2,''$'')') iray
      call EZwrit(850,45, MyString, 0,0)
      write(MyString, '(''MA:'',f6.2,''$'')') Cval
      call EZwrit(650,335, MyString, 0,0)
      write(MyString, '(''Bt:'',f6.2,''$'')') Btesl
      call EZwrit(820,335, MyString, 0,0)
100   continue
!     -                                 Ray path in poloidal cross section:
      xpol = max( zmaxx , rmaj-rminn , rmaxx-rmaj )
      call EZrnd2 ( xpol, xpol, ipol, ipoli )
      call EZsets ( 700,  950,  425,  675,                               &  
     &            -xpol, xpol,-xpol, xpol, 1  )
      call GNsets ( 680, 1020,  380,  720,                               &  
     &            -xpol, xpol,-xpol, xpol, 1  )
      ipsisur=1
      Npts=npsipts(ipsisur)
      do 105 i=1,Npts
         xxp(i) = rpsi(ipsisur,i)
         yyp(i) = zpsi(ipsisur,i)
105   continue
      call EZpnts( xxp,yyp,Npts)
      call GNpnts( xxp,yyp,Npts)
      ipsisur=2
      Npts=npsipts(ipsisur)
      do 106 i=1,Npts
         xxp(i) = rpsi(ipsisur,i)
         yyp(i) = zpsi(ipsisur,i)
106   continue
      call EZpnts (xxp,yyp,Npts)
      call GNpnts (xxp,yyp,Npts)
      ipsisur=3
      Npts=npsipts(ipsisur)
      do 107 i=1,Npts
         xxp(i) = rpsi(ipsisur,i)
         yyp(i) = zpsi(ipsisur,i)
107   continue
      call EZpnts( xxp,yyp,Npts )
      call GNpnts( xxp,yyp,Npts )
 
      do 110 i=1,iplt
      xp(i) = yrr(i) - rmaj
110   yp(i) = yzz(i)
      call EZcurv ( xp , yp , iplt )
      call GNcurv ( xp , yp , iplt )
      call EZcros( (Rmag-Rmaj), zmag, 1 )
      call EZaxes (ipol, ipoli, ipol, ipoli)
      call Find4Dmp(dlnPdsX(1,iray), izone, icross ,nfound)
!
      do 112 i=1,nfound
      call EZcros(  (RofRay(icross(i))-Rmaj),                            &  
     &               ZofRay(icross(i)),        1       )
 112  continue
!
!     -                                 Nper vs sqrt of psi:
      call EZsets ( 100, 350,  75, 325,                                  &  
     &             ZERO, ONE,ZERO, pernx, 1 )
!    ^              0.0, 1.0, 0.0, pernx, 1 )
      call GNsets (   0, 340,   0, 380,                                  &  
     &             ZERO, ONE,ZERO, pernx, 1 )
      call EZwrit ( 100, 335,' Nper vs rt-psi $',0,0   )
      call GNtitl (          ' Nper vs rt-psi $'       )
      call EZwrit ( 100, 300,' ~150 -->NB damp$',0,0   )
!     call GNwrit ( 100, 300,' ~150 -->NB damp$',0,0   )
      woc = sqrt(woc2)
      do 120 i=1,iplt
        xxp(i) = yps(i)
         xp(i) = ypr(i)
         yp(i) = ypl(i)
 120  continue
!
!     -         Temp fix for counting problems
                xxp(1)=xxp(2)
                 xp(1)= xp(2)
                 yp(1)= yp(2)
!     -         Temp fix for counting problems
!
      call EZcurv (xxp, xp, iplt )
      call GNcurv (xxp, xp, iplt )
      call EZaxes ( incx,incxi,inc, 10/inc )
!
!     -                                 Npar vs sqrt of psi:
      call EZsets ( 100, 350, 425, 675,                                  &  
     &             ZERO, ONE, -parnx, parnx, 1 )
!    ^              0.0, 1.0, -parnx, parnx, 1 )
      call GNsets (   0, 340, 380, 760,                                  &  
     &             ZERO, ONE, -parnx, parnx, 1 )
      call EZcurv (xxp, yp, iplt )
      call GNcurv (xxp, yp, iplt )
      write(MyString, '(''N_ll:'',f5.2,''   N_po:$'')')yp(1)
      call EZwrit ( 100, 710,MyString,0,0  )
      call GNtitl (          ' Npar vs rt-psi $'       )
      write(MyString, '(''N_to:'',f5.2,                                  &  
     &                  ''  '',f5.1,''$'')')enpar,enpol
      call EZwrit ( 100, 685,MyString,0,0  )
 
      call EZaxes ( incx,incxi,inc, 10/inc )
 
!     -                                 Npar of good damping vs sqrt of psi:
!
!     .                                 Indicates where damping is expected
!     .                                 Nll = sqrt(511kev/Tekev) / 3.5
!     .                                 w/kv less than 4 ... some damping
!     .                                 w/kv less than 3 ... pretty strong
!     w/kv      zet^3 exp(-zet^2)
!      1.             .21000
!      1.7            .41000   (max!)
!      2.             .37000
!      3.             .11000
!      3.5            .03300   (use this value for Landau damping indication)
!      4.             .00700
!      5.             .00010
      do 130 i=1,NpsiJ
         xxp(i) = ( PsiVec(i)-PsiVec(1) )/( PsiVec(NpsiJ)-PsiVec(1) )
         if (xxp(i) .lt. 0.00_R8) then
           write(MyString,'(''PsiNorm <0 in PlRay, i='',i3)')i
           call LSCwarn(MyString)
         endif
         xxp(i) = sqrt( abs(xxp(i)) )
          yp(i) = 6.5_R8/sqrt(TeeVec(i)+0.01_R8)
         yyp(i) =-6.5_R8/sqrt(TeeVec(i)+0.01_R8)
         if (yp(i) .gt. parnx) go to 131
 130  continue
 
 131  call EZpnts (xxp, yp, i )
      call EZpnts (xxp,yyp, i )
      call GNpnts (xxp, yp, i )
      call GNpnts (xxp,yyp, i )
 
!
!     -                                 Quantities having to do with damping
!     -                                 and polarization
!     -
      if (scatKdeg .le. 0.001_R8) then
!
         continue
!        May 2000, code taken out, placed in Nth-Nrho-snip.F
!        That code was for NJ Fisch for alpha channeling
!
      else if (scatKdeg .gt. 0.001_R8) then
 
!     -                                 Incident and scattered theta plots
!     -                                 if we have scattering.
      do 300 i=1,iscatplt
         xp(i) = AREAL(i)
         yp(i) = inciThet(i)
        yyp(i) = scatThet(i)
 300  continue
      nicemx = xp(iscatplt)
      if(nicemx .gt. AREAL(NPLTDIM)/4.0_R8) then
        nicemx = AREAL(NPLTDIM)
      else if (nicemx .lt. 10._R8) then
        nicemx = 10._R8
        i      = 1
        j      = 10
      else
        call EZrnd2(nicemx, nicemx, i,j)
      endif
      call EZsets ( 425, 600,  75, 325,                                  &  
     &           ZERO, nicemx, ZERO, ONE80, 1 )
!    ^            0.0, nicemx, 0.0 , 180., 1 )
!     call EZaxes (i,j,2,3)
      call EZaxes (i,j,2,1)
      call EZpnts (xp, yp,iscatplt)
      call EZcros (xp,yyp,iscatplt)
!     call EZwrit ( 425, 335, 'inThet. scatThet crs$',0,0)
      call EZwrit ( 425, 335, 'Thet i/o ./x$',0,0)
 
      call EZsets ( 425, 600, 425, 675,                                  &  
     &             ZERO, ONE80,ZERO,ONE80, 1)
!    ^              0.0, 180., 0.0 , 180. ,1)
!     call EZaxes (2,3,2,3)
      call EZaxes (2,1,2,1)
      call Ezcros (yyp, yp,iscatplt)
!     call EZwrit ( 425, 685,'Scatt(x) for Inci(y) $',0,0  )
      call EZwrit ( 425, 710,'Scat(x) for $',0,0  )
      call EZwrit ( 425, 685,'Inci(y)$',0,0  )
 
      endif
 
      igraph = igraph + 1
      call MkGrfLst(' Ray plots in r and k spaces  ')
      call EZfini(0,0)
!     -    EZfini(0,0)  first 0 means more graphs coming; second 0 ignored
!
           call GNfini
      if (iDoAcces .eq. 1) then
        iDoAcces = 0
        igraph = igraph + 1
        ep = .005_R8
        call AccesiSB(RlcfsMin-ep,  RlcfsMax+ep,                         &  
     &              -(ZlcfsMax+ep),+ZlcfsMax+ep,  Rmag)
        call MkGrfLst   (' Cntrs  psi, access, enhance')
      endif
      return
      END
!     .                                                                |
!     .                                                                |
!     PlRay     ends                        ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
