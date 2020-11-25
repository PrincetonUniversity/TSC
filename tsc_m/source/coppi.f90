      subroutine coppi( chiemks, chiimks, chiaux, chiohm,                &  
     &  ane0,anel,xmag,tfluxb,tflux,ptot,forme,formi,vprime,gps,         &  
     &  q95,arad,gzero,zeff,alphar,a121,a122,a123,                       &  
     &  a124,a126,a3015,tfluxs,lsaw,dpsi,a3003,a3011,a3004,a3005,a3006,  &
     &  lhmode,glob132)
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
!     a126      - ratio of chi-i to chi-e (2.0)
!                - if negative, then abs[a126] is chi-i in m**2/sec
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER lsaw,lhmode
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 chiimks,chiaux,chiohm,ane0,anel,xmag,tfluxb,tflux,ptot
      REAL*8 forme,formi,vprime,gps,q95,arad,gzero,zeff,alphar,a121
      REAL*8 a122,a123,a124,a126,tfluxs,dpsi,a3003,a3011,a3004,a3015
      REAL*8 a3005,a3006,glob132,chiemks,pi,q95pct,alphaq,rps
      REAL*8 formfs,formis,befoe,befoi,chi0,chiauxs,chiohms,fhmode
      REAL*8 fitb,coef1,coef2,coef3,xchi
!============
      data pi/3.1415926535_R8/
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
!....this was changed back on 5/26/2010 to enable thermal quench simulations
!     befoi = 8*pi**2*formis*(ane0/(anel*vprime))                        &
      befoi = 8*pi**2*formfs*(ane0/(anel*vprime))                        &
     & *xmag*tfluxb*exp(.667_R8*alphaq*tflux/tfluxb)
!
      chiaux = a121*(7.5E8_R8)*(ptot/ane0)**0.6_R8                       &  
     &       /((gzero*q95pct)**(0.8_R8)*arad**0.2_R8)
      chiohm = a122*1.25E20_R8/ane0*arad*gzero**0.3_R8*zeff**0.2_R8      &  
     &    *(1.0_R8+0.25_R8*alphar)/(xmag**2.2_R8*q95pct**1.6_R8)
!
      chi0 = sqrt(chiohm**2 + chiaux**2)
      chiauxs = chiaux
      chiohms = chiohm
!
!
!.....reduce electron and ion transport near plasma edge for H-mode
      fhmode = 1._R8
!	Changed on 01/19/2011 (fmp) to allow steeper raise of thermal conductivity.
!	This is to avoid too high boostrap current near the edge.
      if(tflux .ge. tfluxb*a3011 .and. a3003 .ne. 0._R8) fhmode =       &
     &        a3003*exp(a3015*((tflux-(real(lhmode-1._R8)*dpsi))        &
     &         /(tfluxb-(real(lhmode-1._R8)*dpsi)))**2)
!

!...reduce electron and ion confinement over a region for ITBs
      fitb=1._R8
      if(tflux .ge. a3004*tfluxb .and. tflux .le. a3005*tfluxb           &  
     & .and. a3006 .ne. 0._R8) then
      fitb=1._R8+ a3006*(1._R8-(tflux/(a3005*tfluxb)))**2
      endif
!
!......no convection,  factor of abs(acoef(126)) less ion conduction
      coef1 = 0._R8
      coef2 = -befoe*chi0 *fhmode*fitb
      if(a126.gt.0) then
      coef3 = -abs(a126)*befoi*chi0 *fhmode*fitb
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
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
