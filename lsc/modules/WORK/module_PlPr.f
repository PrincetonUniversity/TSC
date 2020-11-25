c     PlPr.inc ---------------------------------------------------------
c     Plot flags; Print flags; window information
      INTEGER nPlFlg, nPrFlg, NDIAG
      PARAMETER (nPlFlg =  8, nPrFlg =  3)
      PARAMETER (NDIAG = 10)
      INTEGER
     ^     RAYPL,    SPECPL,   RFDPL,    RFDPSPL,  DAMPL,    JRFPL,
     ^     PITPRFPL, DQLFEPL
      COMMON /PlotICom/
     ^     RAYPL,    SPECPL,   RFDPL,    RFDPSPL,  DAMPL,    JRFPL,
     ^     PITPRFPL, DQLFEPL
c
      INTEGER
     ^     RAYWR,    FSTFRCWR, NPAPWRWR
      COMMON /PrinICom/
     ^     RAYWR,    FSTFRCWR, NPAPWRWR
 
 
 
c     set PlFlg (----PL  ) = TRUE         to : Plot the  quantity
c     set PrFlg (----WR  ) = TRUE         to : Write the quantity
 
c PLOT FLAGS
c  1             RAYPL     Ray position in r,z,phi and in kperp,kparallel; enhancement
c  2             SPECPL    Launched spectrum vs npar and v
c  3             RFDPL     Pray,Jray at 12 psi's, 6 per page; total Pray Jray vs npar & vpar
c  4             RFDPSPL   Pray Pql Jrf and integrals vs rtPsi
c  5             DAMPL     ql absn vs rt psi, 6 rays per page
c  6             JRFPL     Ne Te Run Edc Itsc dJ/dE E/J
c  7             PITPRFPL  Bz/Bphi, pitch, n// enhance, Ne&Prf; vs Rmajor
c  8             DQLFEPL
 
c WRITE FLAGS
c  1             RAYWR     ray stuff
c  2             FSTFRCWR  fractions of fast particles
c  3             NPAPWRWR  n-parallels and powers in rays vs index
 
      INTEGER
     ^        FeWind(NWINDIM), DqlWind(NWINDIM),
     ^        PlFlg(nPlFlg), PrFlg(nPrFlg), idiag(NDIAG)
      COMMON /PlPrICom/
     ^        FeWind         , DqlWind,
     ^        PlFlg, PrFlg, idiag
c
c     .          idiag     array of numbers stating at which ramp up
c     .                    iteration requested diagnostics will be
c     .                    performed
c     PlPr.inc ends ----------------------------------------------------
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
