      subroutine kcoppi( chiemks, chiimks, chiaux, chiohm,               &  
     &  ane0,anel,xmag,tfluxb,tflux,ptot,forme,formi,vprime,gps,         &  
     &  q95,arad,gzero,zeff,alphar,a121,a122,a123,                       &  
     &  a124,a126,tfluxs,lsaw,dpsi,a3003,a3011,a3004,a3005,a3006,lhmode,  &  
     &  glob132)
!.......................................................................
!.....coppi/tang profile consistency model
!         REF: S.C.Jardin, M.G.Bell,N.Pomphrey, Nucl Fusion 33 p371(1993)
!
!.......................................................................
!
!      DESCRIPTION OF RETURNED QUANTITIES:
!
!
!     chiemks   - electron thermal conductivity in m**2/sec
!     chiimks   - ion      thermal conductivity in m**2/sec
!     chiaux    - relative size of auxialliary heated part (diagnostic)
!     chiohm    - relative size of ohmic heated part (diagnostic)
!
!     DESCRIPTION OF INPUT QUANTITIES:
!
!     ane0      - central electron density (1/m**3)
!     anel       - electron density at this surface (1/m**3)
!     xmag      - magnetic axis location (m)
!     tflux     - toroidal flux (W)
!     tfluxb    - toroidal flux at p/v boundary (W)
!     tfluxs    - toroidal flux at sawtooth surface [eg. q=1] (W)
!     ptot      - total input power to plasma, ohmic+auxialliary (Watts)
!     forme     - integrated power to electrons inside surface tflux
!     formi     - integrated power to ions inside surface tflux
!     forme     - normalized heating power (including Ohmic heating) to
!     formi       electrons (forme) and ions(formi) inside surface tflux
!                 [see P(tflux)/P(tfluxb) in eq 26 of Ref]
!     vprime    - differential volume wrt toroidal flux (m**3/W)
!     gps       - [grad(tflux)]**2  [W/m]**2
!     q95       - safety factor at the 95% flux surface
!     arad      - approximate minor radius
!     gzero     - R times the vacuum toroidal field (T-m)
!     zeff      - effective charge state at this surface
!     alphar    - density exponent
!     a121      - transport multiplier for auxialliary heated part
!     a122      - transport multiplier for ohmic part
!     a123      - constant added to q95 (0.5)
!     a124      - enhancement factor inside sawtooth surface (2.0)
!     a126      - ratio of chi-i to chi-e (2.0)
!                - if negative, then abs[a126] is chi-i in m**2/sec
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER lsaw,lhmode,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 chiimks,chiaux,chiohm,ane0,anel,xmag,tfluxb,tflux,ptot
      REAL*8 forme,formi,vprime,gps,q95,arad,gzero,zeff,alphar,a121
      REAL*8 a122,a123,a124,a126,tfluxs,dpsi,a3003,a3011,a3004
      REAL*8 a3005,a3006,glob132,chiemks,rho,chieref,chiiref,pi
      REAL*8 q95pct,alphaq,rps,formfs,formis,befoe,befoi,chi0
      REAL*8 chiauxs,chiohms,fsaw,fhmode,fitb,coef1,coef2,coef3
      REAL*8 xphi,xrho,xrho1,ratio
!============
      dimension rho(20),chieref(20),chiiref(20)
      data rho/                                                          &  
     &3.556259_R8,7.099242_R8,1.061974E+01_R8,1.411509E+01_R8,           &
     &1.754770E+01_R8,                                                            &  
     &2.091405E+01_R8,                                                   &  
     &2.421049E+01_R8,2.743326E+01_R8,3.059613E+01_R8,3.372893E+01_R8,   &  
     & 3.681868E+01_R8,                                                  &  
     &3.983904E+01_R8,                                                   &  
     &4.276131E+01_R8,4.555545E+01_R8,4.819492E+01_R8,5.065723E+01_R8,   &  
     & 5.292394E+01_R8,                                                  &  
     &5.498439E+01_R8,                                                   &  
     &5.683379E+01_R8,5.846668E+01_R8/
      data chiiref/                                                      &  
!    +9.501198E+05,3.196480E+05,7.977266E+04,4.581731E+04,4.430783E+04,
     &1.000000E+05_R8,1.000000E+05_R8,7.977266E+04_R8,4.581731E+04_R8,   &  
     & 4.430783E+04_R8,                                                  &  
     &5.804945E+04_R8,                                                   &  
     &6.606377E+04_R8,5.000729E+04_R8,3.289179E+04_R8,2.068061E+04_R8,   &  
     & 1.396197E+04_R8,                                                  &  
     &1.082066E+04_R8,                                                   &  
     &9.403654E+03_R8,8.890910E+03_R8,9.418563E+03_R8,1.148820E+04_R8,   &  
     & 1.537367E+04_R8,                                                  &  
     &1.994862E+04_R8,                                                   &  
!    +2.300151E+04,1.650565E+03/
     &2.300151E+04_R8,2.000000E+04_R8/
      data chieref/                                                      &  
!    +8.416930E+04,1.585817E+05,2.369681E+05,4.197183E+05,8.024711E+05,
!    +6.043049E+05,
     &8.416930E+04_R8,1.585817E+05_R8,2.369681E+05_R8,4.197183E+05_R8,   &  
     & 4.000000E+05_R8,                                                  &  
     &4.000000E+05_R8,                                                   &  
     &2.677928E+05_R8,2.589918E+05_R8,1.564898E+05_R8,1.001236E+05_R8,   &  
     & 9.462502E+04_R8,                                                  &  
     &9.019772E+04_R8,                                                   &  
     &8.586827E+04_R8,7.713179E+04_R8,6.122165E+04_R8,5.160074E+04_R8,   &  
     & 3.784617E+04_R8,                                                  &  
     &4.034664E+04_R8,                                                   &  
     &5.417601E+04_R8,9.484745E+04_R8/
      data pi/3.1415926535_R8/
!============      
      q95pct = q95
      if(q95pct.gt.6.0_R8) q95pct=6.0_R8
      if(q95pct.lt.2.0_R8) q95pct=2.0_R8
      alphaq = (q95pct + a123)
      if(alphaq .lt. 2.0_R8) alphaq = 2.0_R8
      rps = 1.E-4_R8
      formfs = forme
      formis = formi
      if(formfs.gt.1.0_R8) formfs = 1.0_R8
      if(formfs.lt.rps) formfs = rps
      if(formis.gt.1.0_R8) formis = 1.0_R8
      if(formis.lt.rps) formis = rps
!
      befoe = 8*pi**2*formfs*(ane0/(anel*vprime))                        &  
     & *xmag*tfluxb*exp(.667_R8*alphaq*tflux/tfluxb)
!
      befoi = 8*pi**2*formis*(ane0/(anel*vprime))                        &  
     & *xmag*tfluxb*exp(.667_R8*alphaq*tflux/tfluxb)
!
      chiaux = a121*(7.5E8_R8)*(ptot/ane0)**0.6_R8                       &  
     &       /((gzero*q95pct)**(0.8_R8)*arad**0.2_R8)
!     chiaux = 0.5*27.6/glob132
      chiohm = a122*1.25E20_R8/ane0*arad*gzero**0.3_R8*zeff**0.2_R8      &  
     &    *(1.0_R8+0.25_R8*alphar)/(xmag**2.2_R8*q95pct**1.6_R8)
!
      chi0 = sqrt(chiohm**2 + chiaux**2)
      chiauxs = chiaux
      chiohms = chiohm
!
!
!.....enhance electron transport for sawtooth model
      fsaw = 1._R8
!     if(tflux .le. tfluxs) fsaw = a124
      if(tfluxs .gt. 0.0_R8) then
      if(tflux .eq. (lsaw+2)*dpsi) fsaw = 0.05_R8*a124 + 0.95_R8
      if(tflux .eq. (lsaw+3)*dpsi) fsaw = 0.025_R8*a124 + 0.975_R8
      if(tflux .eq. (lsaw+1)*dpsi) fsaw = 0.1_R8*a124 + 0.9_R8
      if(tflux .eq. (lsaw)*dpsi) fsaw = 0.3_R8*a124 + 0.7_R8
      if(tflux .eq. (lsaw-1)*dpsi) fsaw = 0.5_R8*a124 + 0.5_R8
      if(tflux .eq. (lsaw-2)*dpsi) fsaw = 0.7_R8*a124 + 0.3_R8
      if(tflux .eq. (lsaw-3)*dpsi) fsaw = 0.9_R8*a124 + 0.1_R8
      if(tflux .lt. (lsaw-3)*dpsi) fsaw = a124
      endif
!
!.....reduce electron and ion transport near plasma edge for H-mode
      fhmode = 1._R8
!     if(tflux .ge. tfluxb*a3011 .and. a3003 .ne. 0.) fhmode = a3003
      if(a3003 .ne. 0) then
      if(tflux .eq. (lhmode-5)*dpsi) fhmode=0.025_R8*a3003+0.975_R8
      if(tflux .eq. (lhmode-4)*dpsi) fhmode=0.050_R8*a3003+0.950_R8
      if(tflux .eq. (lhmode-3)*dpsi) fhmode=0.075_R8*a3003+0.925_R8
      if(tflux .eq. (lhmode-2)*dpsi) fhmode=0.100_R8*a3003+0.900_R8
      if(tflux .eq. (lhmode-1)*dpsi) fhmode=0.300_R8*a3003+0.700_R8
      if(tflux .eq. (lhmode)*dpsi) fhmode=0.500_R8*a3003+0.500_R8
      if(tflux .eq. (lhmode+1)*dpsi) fhmode=0.700_R8*a3003+0.300_R8
      if(tflux .eq. (lhmode+2)*dpsi) fhmode=0.900_R8*a3003+0.100_R8
      if(tflux .eq. (lhmode+3)*dpsi) fhmode=0.925_R8*a3003+0.075_R8
      if(tflux .eq. (lhmode+4)*dpsi) fhmode=0.950_R8*a3003+0.050_R8
      if(tflux .eq. (lhmode+5)*dpsi) fhmode=0.975_R8*a3003+0.025_R8
      if(tflux .gt. (lhmode+5)*dpsi) fhmode=a3003
!     if(tflux .eq. (lhmode)*dpsi) fhmode=0.500*a3003+0.500
!     if(tflux .eq. (lhmode+1)*dpsi) fhmode=0.600*a3003+0.400
!     if(tflux .eq. (lhmode+2)*dpsi) fhmode=0.700*a3003+0.300
!     if(tflux .eq. (lhmode+3)*dpsi) fhmode=0.800*a3003+0.200
!     if(tflux .eq. (lhmode+4)*dpsi) fhmode=0.900*a3003+0.100
!     if(tflux .eq. (lhmode+5)*dpsi) fhmode=0.937*a3003+0.063
!     if(tflux .eq. (lhmode+6)*dpsi) fhmode=0.975*a3003+0.025
!     if(tflux .gt. (lhmode+6)*dpsi) fhmode=a3003
!     if(tflux .eq. (lhmode-1)*dpsi) fhmode=0.650*a3003+0.350
!     if(tflux .eq. (lhmode)*dpsi) fhmode=0.750*a3003+0.250
!     if(tflux .eq. (lhmode+1)*dpsi) fhmode=0.850*a3003+0.150
!     if(tflux .eq. (lhmode+2)*dpsi) fhmode=0.950*a3003+0.050
!     if(tflux .gt. (lhmode+2)*dpsi) fhmode=a3003
      endif
!
!...reduce electron and ion confinement over a region for ITBs
      fitb=1._R8
      if(tflux .ge. a3004*tfluxb .and. tflux .le. a3005*tfluxb           &  
     & .and. a3006 .ne. 0._R8) then
      fitb=1._R8+ a3006*(1._R8-(tflux/(a3005*tfluxb)))**3.0_R8
      endif
!     if(tflux .ge. a3011*tfluxb .and. tflux .le. a3005*tfluxb
!    + .and. a3006 .ne. 0.) fitb=a3006
!     if(tflux .ge. a3011*tfluxb .and. tflux .le. a3005*tfluxb
!     if(a3006 .ne. 0.) then
!     tfluxave=(a3011+a3005)*tfluxb/2.
!     deltf=(a3005-a3011)*tfluxb/2.
!     fitb=1.-(1.-a3006)*exp(-(tflux-tfluxave)**2/(deltf**2))
!     endif
!
!
!......no convection,  factor of abs(acoef(126)) less ion conduction
      coef1 = 0._R8
!     coef2 = -befoe*chi0*fsaw
      coef2 = -befoe*chi0*fsaw*fhmode*fitb
      if(a126.gt.0) then
      coef3 = -abs(a126)*befoi*chi0*fsaw*fhmode*fitb
      else
!
!....change made 11/4/97....scj
      coef3 = min(a126*gps , coef2)
      endif
!
!.....electron thermal conductivity in m**2/sec
      if(gps.ne.0) then
      chiemks = -coef2/gps
      chiimks = -coef3/gps
      endif
!
      do 9090 i=1,19
      xphi = (tflux/tfluxb)**a3003
      xrho = rho(i)/rho(20)
      xrho1 = rho(i+1)/rho(20)
      if(xphi .gt. xrho .and. xphi .le. xrho1) then
      chiemks = chieref(i) + ((chieref(i+1)-chieref(i))/(xrho1-xrho))    &  
     & *(xphi-xrho)
      chiimks = chiiref(i) + ((chiiref(i+1)-chiiref(i))/(xrho1-xrho))    &  
     & *(xphi-xrho)
      go to 9091
      endif
 9090 continue
 9091 continue
      if((tflux/tfluxb) .le. 0.0_R8) chiemks = chieref(1)
      if((tflux/tfluxb) .ge. 1.0_R8) chiemks = chieref(20)
      if((tflux/tfluxb) .le. 0.0_R8) chiimks = chiiref(1)
      if((tflux/tfluxb) .ge. 1.0_R8) chiimks = chiiref(20)
      ratio = 28.6_R8/glob132
      chiemks = ratio*a3005*chiemks/1.E4_R8
      chiimks = ratio*a3006*chiimks/1.E4_R8
!
      if(gps .ne. 0._R8) then
      coef2 = -gps*chiemks
      coef3 = -gps*chiimks
      endif
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
