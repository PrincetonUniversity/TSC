c     ---------------------------------- TSC.inc ---------------------|
c                                                                     |
c                                                                     |
      CHARACTER*8  equhd(PWORDS)
      INTEGER
     ^     iplim, isym, kcycle, npsitm, nspc,
     ^     nx, nz, COEFDIM
      INTEGER
     ^     NpsiJ, NpsiM1, NumR, NumZ, NumZh, ISIZElcfs, Nlcfs
      PARAMETER(COEFDIM = 4, ISIZElcfs=800)
      COMMON / TSCi0 /
     ^     iplim, isym, kcycle, npsitm, nspc,
     ^     nx, nz,
     ^     NpsiJ, NpsiM1, NumR, NumZ, NumZh, Nlcfs, equhd
      REAL*8
     ^     PsiGrd(PNX, (2 * PNZ - 1))
      REAL*8
     ^     PsiInc, RBphi0,
     ^     pe2Fac, pi2Fac, pe2Fac14, AioFac, AelFac, ceiFac, OmcFac
      REAL*8
     ^     achrg(PIMP+1), amass(PIMP+1),
     ^     anecc(PPSI), anicc(PPSI,PIMP+1),
     ^     apl, bgzero, gary(PPSI), gpary(PPSI),
     ^     pary(PPSI), ppary(PPSI), psep(2), xsep(2), zsep(2),
     ^     rgzero, rho(PPSI), tekev(PPSI),
     ^     tikev(PPSI,PIMP+1), times, voltlp(PPSI),
     ^     vptemp(PPSI), xary(PNX), xmag, zmag, xsv(PPSI),
     ^     zary(2*PNZ-1)
      REAL*8
     ^     dVlVec(PPSI ), iVlVec(PPSI ),
     ^     pi2Vec(PPSI ), AioVec(PPSI ), AelVec(PPSI ),
     ^     pe2Vec(PPSI ), RBpVec(PPSI ), TeeVec(PPSI ),
     ^     VprVec(PPSI ), PsiVec(PPSI ), EdcVec(PPSI ), MidVec(PPSI)
      REAL*8
     ^     pi2Coefs(COEFDIM), AioCoefs(COEFDIM ),   AelCoefs(COEFDIM),
     ^     pe2Coefs(COEFDIM), RBphiCoefs(COEFDIM), TeeCoefs(COEFDIM)
 
      REAL*8 Rlcfs(ISIZElcfs), Zlcfs(ISIZElcfs),
     ^     RlcfsMax, RlcfsMin, ZlcfsMin,     ZlcfsMax
 
      EQUIVALENCE (RBpVec(1), gary(1)),   (TeeVec(1), tekev(1))
      EQUIVALENCE (VprVec(1), vptemp(1)), (PsiVec, xsv)
 
      COMMON / TSCr0 /
     ^     PsiInc, RBphi0,
     ^     pe2Fac, pi2Fac, pe2Fac14, AioFac, AelFac, ceiFac, OmcFac
      COMMON / TSCr1 /
     ^     PsiGrd,
     ^     achrg, amass,
     ^     anecc, anicc,
     ^     apl, bgzero, gpary,
     ^     pary, ppary, psep, xsep, zsep,
     ^     rgzero, rho,
     ^     tikev, times, voltlp,
     ^     vptemp, xary, xmag, zmag,
     ^     zary,   dVlVec, iVlVec,
     ^     pi2Vec, AioVec, AelVec,
     ^     pe2Vec, RBpVec, TeeVec, PsiVec, EdcVec, MidVec,
     ^     pi2Coefs, AioCoefs, AelCoefs,
     ^     pe2Coefs, RBphiCoefs, TeeCoefs
 
      COMMON / TSClcfs /
     ^     Rlcfs, Zlcfs,
     ^     RlcfsMax, RlcfsMin, ZlcfsMin, ZlcfsMax
c     lcfs: last closed flux surface
c     This needs to be established in the equilibrium read because
c     it can happen that the rmax rmin passed are much different than
c     the min/max of the plasma.
c
c                                                                     |
c                                                                     |
c     ---------------------------------- TSC.inc ---------------------|
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
