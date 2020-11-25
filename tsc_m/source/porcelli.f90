      subroutine porcelli(r1,ak1,ad1,a,s1, deltaq,                       &  
     &    betap1, cp, betai0, betaalpha, etpara,                         &  
     &    va, r, rhoi, betai1, rp, rn, tau, omegastari,                  &  
     &    valphen, omegastara, gammap, dwcore,dwcorec,                   &  
     &    dw,dwc,dwc2, ratio, scrit,dwbussac,dwelong,dwko,dwfast,        &  
     &    dwideal)
!
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 ak1,ad1,a,s1,deltaq,betap1,cp,betai0,betaalpha,etpara
      REAL*8 va,r,rhoi,betai1,rp,rn,tau,omegastari,valphen
      REAL*8 omegastara,gammap,dwcore,dwcorec,dw,dwc,dwc2,ratio
      REAL*8 scrit,dwbussac,dwelong,dwko,dwfast,dwideal,r1,ch,crho
      REAL*8 cstar,pi,amu0,taueta,taua,s,eps,rhohat,alpha,cmhd
      REAL*8 betac,cel
!============
      data ch,crho,cstar / 0.4_R8, 1.0_R8, 3.0_R8/
      data pi / 3.1415926535_R8/
      amu0 = 4*pi*1.E-7_R8
!
!.....input variables:
!
!     r1   average radius of q=1 surface
!     ak1  average ellipticity of q=1 surface
!     ad1  average triangularity of q=1 surface
!     a    average minor radius of plasma boundary
!     s1   magnetic shear at q=1 surface
!     deltaq = 1 - q_0
!     betap1 = poloidal beta inside q=1 surface (see eq 12)
!     cp   special integral of pressure inside q=1 surface
!     betai0 = peak ion toroidal beta
!     betaalpha  special integral of pressure    (see eq 27)
!     etpara  eta parallel on axis
!     va   Alfven speed on axis
!     r    major radius
!     rhoi = ion Larmor radius
!     betai1 = ion toroidal beta at the q=1 surface
!     rp = pressure scale length
!     rn = density scale length
!     tau = te/ti
!     omegastari ... drift parameter for the ions
!     valphen    ... Alfven velocity
!     omegastara ... alpha-particle prec frequency
 
!.....output variables
!
!     scrit  shear at the q=1 surface
!     dwcore
!     dwcorec
!     dw
!     dwc
!     dwc2
!     ratio           Equation 15:   ratio needs to be .le. 1 to trigger a sawtooth
!     dwbussac
!     dwelong
!     dwko
!     dwfast
!
!
      taueta = r1**2*amu0/etpara
      taua = sqrt(3._R8*r)/va
      s = taueta/taua
      eps = r1 / r
      rhohat = rhoi / r1
!
!
!.....Critical shear at rational surface
      alpha = 1.5_R8*cstar**(-7._R8/6._R8)*(tau/(1._R8+tau))**(7._R8/    &  
     & 12._R8)
      scrit = alpha * sqrt(s**(1._R8/3._R8)*rhohat)*                     &  
     &        (betai1*r**2/r1**2)**(7._R8/12._R8)*(r1/rn)*(r1/rp)**      &  
     & (1._R8/6._R8)
!
!.....Ideal MHD term
      cmhd = 3._R8*pi/2
      betac = 0.3_R8*(1._R8- (5*r1/(3._R8*a)))
      dwbussac = -cmhd*eps**2*(betap1**2 - betac**2)
!
!.....Elongation term
      cel = (pi/3._R8)*deltaq**2
!
!.....Ideal MHD term: [Eriksson and Wahlberg, Phys Plasmas, 9 (2002) 1606]
      dwideal = .5_R8*pi*eps**2/deltaq *                                 &  
     &   (- .75_R8*(ak1-1._R8)*betap1*(1._R8-2._R8*ad1/eps)              &  
     &   + deltaq*(13._R8/48._R8- 3._R8*betap1**2 + .5_R8*(ak1-1)*       &  
     &  (13._R8*betap1**2-0.25_R8*betap1 - 1._R8+ ad1/(6._R8*eps)*       &  
     & (4._R8*betap1-7._R8))))
!
      dwelong = -cel*((ak1-1._R8)/2._R8)**2
!
!.....Kruskal-Oberman term
      dwko = 0.6_R8*cp*sqrt(eps)*betai0 / s1
!
!.....core delta-W
      dwcore = dwbussac + dwelong + dwko
!.....critical value
      dwcorec = ch*omegastara*taua
!
!.....Fast particle term
      dwfast = eps**1.5_R8*betaalpha/s1
!
!.....Total delta-W
      dw = dwcore + dwfast
!
!.....critical value
      dwc = 0.5_R8*omegastari*taua
      dwc2=  crho*rhohat
!
!.....final criteria ...eq (15)
      ratio = omegastari/(cstar*gammap)
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
