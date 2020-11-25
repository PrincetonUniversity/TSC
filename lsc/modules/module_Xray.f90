      MODULE Xray
      USE PARAMS
      USE EMPARAMS
      USE XPARAMS
      IMPLICIT NONE
!
!     .         --------------------------------------------------------
!     Xray.inc
!     E_photon        in MeV
!     E_incvec        in MeV
!     E_min E_max
      CHARACTER*2 FoilCode
      INTEGER nMUbins, nEbins, iAbsXray
      INTEGER Eint(NVELDIM)
      REAL*8                                                             &  
     &        mu_min, mu_max, E_min, E_max, E_ph_min, dE_ph,             &  
     &        dFoilTCM
      REAL*8    dmu_inv
      REAL*8                                                             &  
     &        sigtot(NMUDIM, NENDIM), E_incvec(NENDIM), XmnFac(NENDIM),  &  
     &        muvec(NMUDIM), Efrac(NVELDIM), Eofv(NVELDIM),              &  
     &        dEofv(NVELDIM)
 
      REAL*8                                                             &  
     &        inten(NMUDIM, NPSIDIM)
 
!     COMMON / Xray1a / nMUbins, nEbins, iAbsXray
!     COMMON / Xray2b / Eint
!     COMMON / Xray3c / sigtot, E_incvec, XmnFac, muvec, inten, Efrac
!     COMMON / Xray4d / E_ph_min, dE_ph
!     COMMON / Xray5e / dFoilTCM
!     COMMON / Xray6f / Eofv, dEofv
!     COMMON / Xray7g / mu_min, mu_max, dmu_inv, E_min, E_max
!     COMMON / Xray8h / FoilCode
!                       ! dFoilTCM; or; dcu: thickness of foil (cm)
!                       ! FoilCode; or; FC: CU,AG,TA,MO, or 00 for no foil
!     Xray.inc
!     .         --------------------------------------------------------
!
!cj   DATA mu_min, mu_max / -1._R8, 1._R8/
!cj   DATA nMUbins / NMUDIM /
!cj   DATA nEbins  / NENDIM /
!cj   DATA E_ph_min / 1.E-6_R8/
!cj   DATA dE_ph / .01_R8/
!cj   DATA E_min, E_max / 0.01_R8, 0.2_R8/
!cj   DATA FoilCode / 'AG' /
!cj   DATA dFoilTCM / 0.000_R8/
!cj   DATA iAbsXray / 1 /



 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE Xray
