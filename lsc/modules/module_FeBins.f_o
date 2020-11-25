      MODULE FeBins
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
      INTEGER, PRIVATE, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

c     FeBins.inc ------------------------------------------------------|
c                                                                      |
c                                                                      |
c     dvsym   REAL    VARRAY
c             delta v's figured symmetrically
c     dvisym  REAL    VARRAY
c             inverse of the delta v's figured symmetrically
c     nv      INTEGER
c             actual number of parallel velocity zones (input)
c     fe      REAL    VQLARRAY
c             (iv, ipsi, iitr)
c             electron distribution function indexed
c             on parallel velocity and space.  Velocity is normalized to
c             c, speed of light.  Integral of fe over v gives density in
c             cm^{-3}.  Jardin data file passes density in this form.
c     FstFracN(NPSIDIM)
c             The fraction of electrons that can be called fast particles
c             for their location because of interaction with the wave.
c             \integral_-1^1 ( fe(v/c,\psi) - FeMaxw(v/c,\psi) ) dv/c
c             divided by
c             \integral_1^1 FeMaxw(v/c,\psi) dv/c
c     FstFracE(NPSIDIM)
c             The fraction of electron energy due to fast particles
c             defined analogously to FstFracN but with v^2 weight.
c     dfdv    (iv, ipsi, 1) un smoothed df/dv
c     dfdv    (iv, ipsi, 2)    smoothed df/dv
c     .
c     TailPeps  Regarding a ficticious fast electron tail in both
c     TailNeps  directions, the fraction eps(ilon) for Pressure, Ne and Te
c     TailTeps  such that
c     .                  TailNeps = TailPeps * TailTeps
c     .                  TailTeps = T_thermal/T_fast
c     .                  TailPeps = (n_fast T_fast)/(n_thermal T_thermal)
c     .         and
c     .              f(v) = (2 pi v_t^2)^-.5 1. n_e exp[- v^2/(2 v_t^2)]; all v
c     .         f_fast(v) = (2 pi v_f^2)^-.5 1. n_f exp[- v^2/(2 v_f^2)]; all v
c     .
c     .         where n_f/n_e = TailNeps
c     .               v_t/v_f = TailTeps^.5
c     .               n_f v_f^2 / [ n_t v_t^2 ] = TailPeps
c     .
c     TailVtrn  is the transition velocity relative to v_t at which the fast
c               electron tail becomes more important
c                        TailVtrn = v/v_t | transition
c     .                           = sqrt[2 ln(1/Neps 1/Teps^.50)/(1 - Teps^2)]
      INTEGER ivZero, nv, iITR, iSMO
      DATA    iSMO,      nv /
     ^           1, NVELDIM /

      REAL*8    WeghtItr,
     ^        Vmin, Vmax, VthNorm,
     ^        TeUnits, dv, fe0, nu0
      DATA TeUnits /1000._R8/
      DATA Vmin, Vmax, WeghtItr /
     ^      -1._R8,  1._R8,     0.20_R8/

c    ^        frv1minus, frv2minus, v1minus, v2minus,
c    ^        frv1plus,   frv2plus, v1plus,  v2plus, epsvgr,
c    ^                             , betaZ, lnLambda
      REAL*8    TailPeps, TailNeps, TailTeps, TailVtrn
      REAL*8
     ^        Vpar(NVELDIM), fe(NVELDIM, NPSIDIM, NITRDIM),
     ^        Vtherm(NPSIDIM), FeNorm(NPSIDIM),
     ^        VperpSq(NPSIDIM), dvsym(NVELDIM), dvisym(NVELDIM),
     ^     dvip(NVELDIM), dvplus(NVELDIM),
     ^     dfdv(NVELDIM, NPSIDIM,2), nu0psi(NPSIDIM),
     ^     FstFracN(NPSIDIM), FstFracE(NPSIDIM)
!     COMMON /feIbin/
!    ^        ivZero, nv, iITR, iSMO
!     COMMON /feRbin/ WeghtItr,
!    ^        Vmin, Vmax, VthNorm,
!    ^        TeUnits, dv, fe0, nu0,
!    ^        Vpar       , fe                , Vtherm         ,
!    ^         FeNorm,     VperpSq,
!    ^        dvsym, dvisym, dvplus, dvip,
!    ^     dfdv, nu0psi,
!    ^     FstFracN, FstFracE,
!    ^ TailPeps, TailNeps, TailTeps, TailVtrn
 
c     ivZero: index for which Vpar(ivZero) = 0.
c     iITR:   index for iteration on fe, either 1 or 2...NITRDIM
c     iSMO:   index for smoothing of Dql, either 1 (un) or 2 (smoothed)
c             note that Dql(v,p,1) goes with df/dv(v,p,2) and that
c                       Dql(v,p,2) goes with df/dv(v,p,1)
c     .                                 iSMO = 2 --> Do smoothing as
c     .                                 EValeo did in original code
c     .                                 using smoothed Dql=Dql(iv,ip,2)
c     .                                 iSMO = 1 -->
c     .                                        Use unsmoothed Dql(.,.,1)
c     .                                 and smoothed dfdv(iv,ip) for
c                                              Pql & Jrf
c
c     nv: total number of velocity grid points
c     Vpar:  parallel velocity / speed of light
c     vthermal:  thermal velocity / speed of light
c     fe:  electron velocity distribution function normalized so that
c     integral fe dVpar = electron number density cm-3
c     TeUnits:  normalization for electron temperature in keV
c                                                                      |
c                                                                      |
c     FeBins.inc ------------------------------------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE FeBins
