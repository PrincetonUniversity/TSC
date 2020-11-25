      MODULE PlPr
      USE PARAMS
      USE EMPARAMS
      USE Ramppwr, ONLY : NRAMPDIM
      IMPLICIT NONE
!     PlPr.inc ---------------------------------------------------------
!     Plot flags; Print flags; window information
      INTEGER nPlFlg, nPrFlg, NDIAG
      PARAMETER (nPlFlg =  8, nPrFlg =  3)
      PARAMETER (NDIAG = 10)
      INTEGER                                                            &  
     &     RAYPL,    SPECPL,   RFDPL,    RFDPSPL,  DAMPL,    JRFPL,      &  
     &     PITPRFPL, DQLFEPL

      DATA                                                               &  
     &     RAYPL,    SPECPL,   RFDPL,    RFDPSPL,  DAMPL,    JRFPL,      &  
     &     PITPRFPL, DQLFEPL                                             &  
     &     /                                                             &  
     &     1,        2,        3,        4,        5,        6,          &  
     &     7,        8                                                   &  
     &     /

!     COMMON /PlotICom/
!    ^     RAYPL,    SPECPL,   RFDPL,    RFDPSPL,  DAMPL,    JRFPL,
!    ^     PITPRFPL, DQLFEPL
!
      INTEGER                                                            &  
     &     RAYWR,    FSTFRCWR, NPAPWRWR
!     COMMON /PrinICom/
!    ^     RAYWR,    FSTFRCWR, NPAPWRWR
 
 
 
!     set PlFlg (----PL  ) = TRUE         to : Plot the  quantity
!     set PrFlg (----WR  ) = TRUE         to : Write the quantity
 
! PLOT FLAGS
!  1             RAYPL     Ray position in r,z,phi and in kperp,kparallel; enhancement
!  2             SPECPL    Launched spectrum vs npar and v
!  3             RFDPL     Pray,Jray at 12 psi's, 6 per page; total Pray Jray vs npar & vpar
!  4             RFDPSPL   Pray Pql Jrf and integrals vs rtPsi
!  5             DAMPL     ql absn vs rt psi, 6 rays per page
!  6             JRFPL     Ne Te Run Edc Itsc dJ/dE E/J
!  7             PITPRFPL  Bz/Bphi, pitch, n// enhance, Ne&Prf; vs Rmajor
!  8             DQLFEPL
 
! WRITE FLAGS
!  1             RAYWR     ray stuff
!  2             FSTFRCWR  fractions of fast particles
!  3             NPAPWRWR  n-parallels and powers in rays vs index
 
      INTEGER                                                            &  
     &        FeWind(NWINDIM), DqlWind(NWINDIM),                         &  
     &        PlFlg(nPlFlg), PrFlg(nPrFlg), idiag(NDIAG)
      DATA PlFlg      /                                                  &  
     &     nPlFlg * 0 /
      DATA                                                               &  
     &     RAYWR,    FSTFRCWR, NPAPWRWR                                  &  
     &     /                                                             &  
     &     1,        2,        3                                         &  
     &     /
      DATA PrFlg      /                                                  &  
     &     nPrFlg * 0 /

      DATA idiag(1),idiag(2),idiag(3),idiag(4),idiag(5) /                &  
     &     NRAMPDIM,      0 ,      0 ,      0 ,      0  /

      DATA idiag(6),idiag(7),idiag(8),idiag(9),idiag(10)/                &  
     &           0 ,      0 ,      0 ,      0 ,      0  /
      DATA FeWind               /                                        &  
     &     1, NVELDIM, 1, NPSIDIM /
      DATA DqlWind              /                                        &  
     &     1, NVELDIM, 1, NPSIDIM /


!     COMMON /PlPrICom/
!    ^        FeWind         , DqlWind,
!    ^        PlFlg, PrFlg, idiag
!
!     .          idiag     array of numbers stating at which ramp up
!     .                    iteration requested diagnostics will be
!     .                    performed
!     PlPr.inc ends ----------------------------------------------------
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE PlPr
