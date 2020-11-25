      MODULE RayWrk
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
      INTEGER, PRIVATE, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     RayWrk --- COMMON block for Working data.--was ComWrk-- ---------|
!                                                                      |
!                                                                      |
      REAL*8                                                             &  
     &        enpar,  enpol,  fghz,                                      &  
     &        HstpLH, nparmax,nparmin, npolmin, npolmax, omega,          &  
     &        ScatKdeg, thet0,                                           &  
     &        TotPwr
      CHARACTER*8 couplers(NGRPDIM)
      CHARACTER*96 EquTitle
      REAL*8    centers(NGRPDIM), widths(NGRPDIM)
      REAL*8    powers(NGRPDIM),  phaseDeg(NGRPDIM)
      REAL*8    ntor(NTORDIM), Spec(NTORDIM), npol   (NPOLDIM)
      REAL*8    scatThet(NPLTDIM), inciThet(NPLTDIM)
      REAL*8                                                             &  
     &        begin, btesl, capo2,                                       &  
     &        ecyc , ecyc2, epsq , icyc , ipsq ,                         &  
     &        kappa, qlim , q0   ,                                       &  
     &        Rmag, rmaj , rmax , rmin ,                                 &  
     &        wcei2, woc2 , woc4 ,                                       &  
     &        zeff , zmin , zmax ,                                       &  
     &        psilim, psimin
      INTEGER                                                            &  
     &        TurnNegs, lfast  , lstop  ,                                &  
     &        mstpLH  , nfreq  , nstep,                                  &  
     &        ntors, npols, iscatplt
      INTEGER nGrps

      DATA    lfast, nfreq / 0,   100/
      DATA    nstep /20000/
      DATA    ntors, npols /NTORDIM,       1/
      DATA    PhaseDeg(1), PhaseDeg(2), PhaseDeg(3) /                    &  
     &           90._R8,   180._R8,  135._R8/
      DATA    couplers(1), couplers(2), couplers(3)          /           &  
     &        'PBXMFAST' , 'PBXMSLOW' , 'TOKDEVAR'           /
      DATA    enpol,  fghz, HstpLH /                                     &  
     &           0._R8,   4.6_R8, .005_R8/
      DATA    nparmax, nparmin, thet0, ScatKdeg, TurnNegs /              &  
     &           5.5_R8,     2.5_R8,   0.0_R8,     0.00_R8,       0  /
      DATA    npolmax, npolmin /                                         &  
     &           +1.0_R8,    -1.0_R8/
      DATA    nGrps  , centers(1), centers(2), centers(3) /              &  
     &           3   ,      4.0_R8,      4.0_R8,      4.0_R8/
      DATA             widths(1) , widths(2) , widths(3)  /              &  
     &                      1.0_R8,      1.0_R8,      1.0_R8/
      DATA             powers(1) , powers(2) , powers(3)  /              &  
     &                      1.0_R8,      0.1_R8,      0.1_R8/





!     COMMON /BkCray/ couplers, EquTitle
!     COMMON /BkRray/
!    ^        enpar,  enpol,  fghz,
!    ^        HstpLH, nparmax,nparmin, npolmin, npolmax, omega,
!    ^        ScatKdeg, thet0,
!    ^        scatThet, inciThet,
!    ^        centers, widths, powers, phaseDeg,
!    ^        ntor, Spec, npol,
!    ^        TotPwr
 
!     COMMON /BkRwrk/
!    ^        begin, btesl, capo2,
!    ^        ecyc , ecyc2, epsq , icyc , ipsq ,
!    ^        kappa, qlim , q0   ,
!    ^        Rmag, rmaj , rmax , rmin ,
!    ^        wcei2, woc2 , woc4 ,
!    ^        zeff , zmin , zmax ,
!    ^        psilim, psimin
 
!     COMMON /BkIray/
!    ^        TurnNegs, lfast  , lstop  ,
!    ^        mstpLH  , nfreq  , nstep  ,
!    ^        ntors, npols, iscatplt,
!    ^        nGrps
 
 
!      enpar    n_{\parallel} launched for ray being worked at the moment
!      enpol    n_{poloidal}  launched for all rays(zero is good enough)
!      fghz     frequency in GHz
!      omega    RF frequency (radians/sec)
!      phaseDeg phasing in degrees between waveguides; only for nGrps < 1
!      h        step size
!      thet0    angle of launch 0=> outside midplane, .25=> top
!      lfast    1 if fast wave is launched; 0 for slow wave
!      nfreq    stock is taken every nfreq steps in ray
!      nstep    number of ray steps allowed
!
!      begin    value of path length to begin ray (0 at start)
!      btesl    B_{T}  field at nominal major radius
!      capo2    \omega_{ce}\omega_{ci} / \omega^2
!               no LH resonance at any density if .lt. 1
!      curnt    toroidal current I_p in MA
!      ecyc     \omega_{ce} / \omega
!      ecyc2    ecyc^2
!      epsq     \omega_{pe}^2 / \omega^2
!      icyc     \omega_{pi} / \omega
!      iscatplt index of the scatter event; used for filling inciThet, scatThet
!      ipsq     \sum_i \omega_{pi}^2 / \omega^2
!      kappa    elongation: \kappa
!      lstop
!      m        step counter in integration
!      psilim   \psi_{lim}  Flux in webers per radian
!      psimin   \psi_{min}
!      qlim     q_{lim}
!      q0       q_{0}
!      rmaj
!      Rmag     magnetic axis....same as xmag in other commons
!      rmax     outer radius of flux grid
!      rmin     inner radius of flux grid
!      scatKdeg degrees by which nperp may be rotated on bounce according to
!               exp{- (dTheta_deg/scatKdeg)^2 }
!               most of the rotations are less than the scatKdeg
!      scatThet scattered theta in degrees (1_180); nrays of these
!      inciThet incident theta in degrees; no scattering means same answer
!      TotPwr   Total power lauched in the whole spectrum, watts
!                = \sum_{iray=1}^{nrays} power(1,iray)
!      wcei2    \omega_{ce} \omega_{cH} \sum_i ( n_i Z_i^2 / n_e) m_H/m_i
!      woc2     \omega^2/c^2
!      woc4     woc2^2
!      zeff     Z_{eff}
!      zmax     upper extent of flux grid
!      zmin     lower extent of flux grid
!
!                                                                      |
!                                                                      |
!     RayWrk --- COMMON block for Working data.------------------------|
 
 
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE RayWrk
