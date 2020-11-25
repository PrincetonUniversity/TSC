      MODULE Xray
      USE PARAMS
      USE EMPARAMS
      USE XPARAMS
      IMPLICIT NONE
c
c     .         --------------------------------------------------------
c     Xray.inc
c     E_photon        in MeV
c     E_incvec        in MeV
c     E_min E_max
      CHARACTER*2 FoilCode
      INTEGER nMUbins, nEbins, iAbsXray
      INTEGER Eint(NVELDIM)
      REAL*8
     ^        mu_min, mu_max, E_min, E_max, E_ph_min, dE_ph,
     ^        dFoilTCM
      REAL*8    dmu_inv
      REAL*8
     ^        sigtot(NMUDIM, NENDIM), E_incvec(NENDIM), XmnFac(NENDIM),
     ^        muvec(NMUDIM), Efrac(NVELDIM), Eofv(NVELDIM),
     ^        dEofv(NVELDIM)
 
      REAL*8
     ^        inten(NMUDIM, NPSIDIM)
 
!     COMMON / Xray1a / nMUbins, nEbins, iAbsXray
!     COMMON / Xray2b / Eint
!     COMMON / Xray3c / sigtot, E_incvec, XmnFac, muvec, inten, Efrac
!     COMMON / Xray4d / E_ph_min, dE_ph
!     COMMON / Xray5e / dFoilTCM
!     COMMON / Xray6f / Eofv, dEofv
!     COMMON / Xray7g / mu_min, mu_max, dmu_inv, E_min, E_max
!     COMMON / Xray8h / FoilCode
c                       ! dFoilTCM; or; dcu: thickness of foil (cm)
c                       ! FoilCode; or; FC: CU,AG,TA,MO, or 00 for no foil
c     Xray.inc
c     .         --------------------------------------------------------
c
 
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE Xray
