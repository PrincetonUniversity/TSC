C     RayWrk --- COMMON block for Working data.--was ComWrk-- ---------|
C                                                                      |
C                                                                      |
      REAL*8
     ^        enpar,  enpol,  fghz,
     ^        HstpLH, nparmax,nparmin, npolmin, npolmax, omega,
     ^        ScatKdeg, thet0,
     ^        TotPwr
      CHARACTER*8 couplers(NGRPDIM)
      CHARACTER*96 EquTitle
      REAL*8    centers(NGRPDIM), widths(NGRPDIM)
      REAL*8    powers(NGRPDIM),  phaseDeg(NGRPDIM)
      REAL*8    ntor(NTORDIM), Spec(NTORDIM), npol   (NPOLDIM)
      REAL*8    scatThet(NPLTDIM), inciThet(NPLTDIM)
      REAL*8
     ^        begin, btesl, capo2,
     ^        ecyc , ecyc2, epsq , icyc , ipsq ,
     ^        kappa, qlim , q0   ,
     ^        Rmag, rmaj , rmax , rmin ,
     ^        wcei2, woc2 , woc4 ,
     ^        zeff , zmin , zmax ,
     ^        psilim, psimin
      INTEGER
     ^        TurnNegs, lfast  , lstop  ,
     ^        mstpLH  , nfreq  , nstep,
     ^        ntors, npols, iscatplt
      INTEGER nGrps
      COMMON /BkCray/ couplers, EquTitle
      COMMON /BkRray/
     ^        enpar,  enpol,  fghz,
     ^        HstpLH, nparmax,nparmin, npolmin, npolmax, omega,
     ^        ScatKdeg, thet0,
     ^        scatThet, inciThet,
     ^        centers, widths, powers, phaseDeg,
     ^        ntor, Spec, npol,
     ^        TotPwr
 
      COMMON /BkRwrk/
     ^        begin, btesl, capo2,
     ^        ecyc , ecyc2, epsq , icyc , ipsq ,
     ^        kappa, qlim , q0   ,
     ^        Rmag, rmaj , rmax , rmin ,
     ^        wcei2, woc2 , woc4 ,
     ^        zeff , zmin , zmax ,
     ^        psilim, psimin
 
      COMMON /BkIray/
     ^        TurnNegs, lfast  , lstop  ,
     ^        mstpLH  , nfreq  , nstep  ,
     ^        ntors, npols, iscatplt,
     ^        nGrps
 
 
C      enpar    n_{\parallel} launched for ray being worked at the moment
C      enpol    n_{poloidal}  launched for all rays(zero is good enough)
c      fghz     frequency in GHz
c      omega    RF frequency (radians/sec)
c      phaseDeg phasing in degrees between waveguides; only for nGrps < 1
c      h        step size
c      thet0    angle of launch 0=> outside midplane, .25=> top
c      lfast    1 if fast wave is launched; 0 for slow wave
c      nfreq    stock is taken every nfreq steps in ray
c      nstep    number of ray steps allowed
c
c      begin    value of path length to begin ray (0 at start)
C      btesl    B_{T}  field at nominal major radius
C      capo2    \omega_{ce}\omega_{ci} / \omega^2
C               no LH resonance at any density if .lt. 1
C      curnt    toroidal current I_p in MA
C      ecyc     \omega_{ce} / \omega
C      ecyc2    ecyc^2
C      epsq     \omega_{pe}^2 / \omega^2
C      icyc     \omega_{pi} / \omega
c      iscatplt index of the scatter event; used for filling inciThet, scatThet
C      ipsq     \sum_i \omega_{pi}^2 / \omega^2
C      kappa    elongation: \kappa
c      lstop
C      m        step counter in integration
C      psilim   \psi_{lim}  Flux in webers per radian
C      psimin   \psi_{min}
C      qlim     q_{lim}
C      q0       q_{0}
C      rmaj
c      Rmag     magnetic axis....same as xmag in other commons
C      rmax     outer radius of flux grid
C      rmin     inner radius of flux grid
c      scatKdeg degrees by which nperp may be rotated on bounce according to
c               exp{- (dTheta_deg/scatKdeg)^2 }
c               most of the rotations are less than the scatKdeg
c      scatThet scattered theta in degrees (1_180); nrays of these
c      inciThet incident theta in degrees; no scattering means same answer
c      TotPwr   Total power lauched in the whole spectrum, watts
c                = \sum_{iray=1}^{nrays} power(1,iray)
C      wcei2    \omega_{ce} \omega_{cH} \sum_i ( n_i Z_i^2 / n_e) m_H/m_i
C      woc2     \omega^2/c^2
C      woc4     woc2^2
C      zeff     Z_{eff}
C      zmax     upper extent of flux grid
C      zmin     lower extent of flux grid
c
C                                                                      |
C                                                                      |
C     RayWrk --- COMMON block for Working data.------------------------|
 
 
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
