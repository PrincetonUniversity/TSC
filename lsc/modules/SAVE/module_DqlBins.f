      MODULE DqlBins
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
c     DqlBins ---------------------------------------------------------|
c                                                                      |
c                                                                      |
c     File: DqlBins,inc starts
c     Dql     REAL    VQLARRAY
c             (iv, ipsi)
c             electron quasilinear diffusion coefficient indexed
c             on parallel velocity and space
c     Smoothing provided by nsmoo and nwmw accorting to
c     qlsm = (normalization) exp [ -(iv - imax)^2/(2 nsmw^2) ]
c     Values of qlsm at edge of smoothing vector are:
c     nsmw=1; nsmoo=3, 1/root(e); nsmoo=5, 1/e^2; nsmoo=7, 1/e^4.5 = 1/100
c
c     How much smoothing? Suppose NV=1999, 1000 bins from zero to light speed.
c     Suppose T= 1--3 kV; n_par = 2--5
c      v_T/c   40 --  80 bin
c     3 v_T/c 120 -- 240 bin
c     1/n_par 200 -- 500 bin
c     nsmw     ?? 30   ?? 10 min
c     nsmoo    ?? 75   ?? 25
      INTEGER DqlBox(4), nsmoo, nsmw, nsmsym
      REAL*8
     ^        DcollNorm, nuNorm, DqlHite, Pwrnorm, DqlNorm
      REAL*8
     ^        Dql(NVELDIM, NPSIDIM,2),
     ^        Dcoll(NVELDIM, NPSIDIM), nuColl(NVELDIM, NPSIDIM),
     ^        qlsm(NVELDIM)

      DATA   nsmoo,    nsmw /
     ^     NSMODEF, NSMWDEF /

!     COMMON /DqlCom0/
!    ^        DcollNorm, nuNorm, DqlHite, Pwrnorm, DqlNorm,
!    ^        Dql                  ,
!    ^        Dcoll                  , nuColl ,
!    ^        qlsm
 
!     COMMON /DqlCom1/
!    ^        DqlBox, nsmoo, nsmw, nsmsym
c     Dql is the quasilinear diffusion coefficient
c     Dql(v,psi,1) is the UNsmoothed version
c     Dql(v,psi,2) is the   SMOOTHED version
c     Dcoll is the collisional diffusion coefficient
c     nuColl is the collisional drag coefficient
c                                                                      |
c                                                                      |
c     DqlBins ---------------------------------------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE DqlBins
