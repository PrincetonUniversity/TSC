      SUBROUTINE buf_ini( lu, idrss, idrse )
      USE CLINAM

      IMPLICIT NONE
      INTEGER :: lu
      CHARACTER*(*) :: idrss, idrse

      SELECT CASE (idrss) 

      CASE ("ibgcm1")

      read(lu)   ibgcm1      ,                                           &  
     &  nx      , nz     , nxm    , nzm    , nlhcd   , nlhcd2  ,iplot  ,  &  
     &  ncycle , kcycle , ires   , nimax  , nzs     ,  irayt , iglobea ,  &  
     &  lujayp   , luwire , luwirer, nxp    , nzp    , ncoil  ,          &  
     &  iprnt  , iplt   , iplt2  , nskip2 , nskipl , nskipr , nis     ,  &  
     &  numsaw , npert  , neqmax , nebeam , ngrmax , itrmod , idata   ,  &  
     &  iimp   , ilte   , nmult  , ilhcd   ,                             &  
     &  ineg   , irst1  , irst2  , intact , icol   , imag   , jmag    ,  &  
     &  idtmin , jdtmin , ipest  , itpass , isym   , ifunc   , idiv   ,  &  
     &  isaw   , nwire  , nzero  , ninc   , igone  , nc0    , isurf   ,  &  
     &  ifunct , ndiv   , ipres  , idens  , nh     , nsepmax, ifk     ,  &  
     &  itot, icnt(:20),  nnqueue, nnsize , nobs   , igroupp, nloop   ,  &  
     &  nreccdf, imovie , npts   , iplim  , nlim   , numfb  , nrecord ,  &  
     &  lrswtch,irunaway, nplrec , lenscr , icube  , iflux  , ibootst ,  &  
     &  nqflag , npsi   , npsit  , npsim  , ivplbp , ivplvi , ivplfr  ,  &  
     &  nopest , itpest , ialpha , iplate , nplate , ibalsw ,nskipsfi,   &  
     &  npsitmx,  jsym  , nslhrt , nslhpc , ifstrt , iicrh  , ifwcd  ,   &  
     &  ifncm1aa, iecrh, ieccd, idefnpe


      CASE ("ibgcm1aa")

      read(lu) ibgcm1aa,                                                 &  
     &  itrmode, ibaloon(1:ppsi), icall, numpf,                          &  
     &  nreg, iseries(:pngroup), indxd1(:pnfeed), idelay(:pnfeed),       &  
     &  lgroup(:pngroup),  jnreg(:pngroup),                              &  
     &  ncnt, itcv, nwprnt, nwmax, kwprnt,                               &  
     &  iexvc(:penx,:penz), nfeedv(:ptpts,:pnfeed),                      &  
     &  icoil(:penx,:penz), iexv(:penx,:penz),iexs(:penx,:penz),         &  
     &  igroupc(:pncoil), ilcoil(:pncoil), ilwire(:pnwire),              &  
     &  iresgs , ngroup ,ngroups ,nogroup(:pngroup),nogroups(:pngroup),  &  
     &  ngroupt,nogroupt(:pngroup),icmaxs(:pngroup),                     &  
     &  icoils(:pnwire,:pngrps),                                         &  
     &  izsep(:pnsep), jzsep(:pnsep), izsepo(:pnsep), jzsepo(:pnsep),    &  
     &  ngrvw1(:pnwire), iwire(:pnwire), jwire(:pnwire), nseg(:pnplat) ,  &  
     &  ngrvw2(:pnwire),ifeed  (:pnwire), igroupw(:pnwire),              &  
     &  ngrvw3(:pnwire),ilima  (:pnlim) ,   jlima(:pnlim) ,              &  
     &  ngrvw4(:pnwire), nfeedo(:pnfeed),           nfeed2 (:pnfeed),    &  
     &  ngrvc1(:pncoil),iobs  (:pobs),           npltobs(:pobs),         &  
     &  ngrvc2(:pncoil),jobs  (:pobs),           iflagfb(:pnfeed),       &  
     &  ngrvc3(:pncoil), nrfb(:pnfeed),           ipext(:pnfeed),        &  
     &  ngrvc4(:pncoil),                                                 &  
     &  nsturb, iptch(:pkw), mloopx, nn, iwtype, nsegs(:pkw), iltch,     &  
     &  numicsuf, impnum(:18),ichargst(:18),                             &  
     &  ifncm1

      CASE ("ibgc1a")

      read(lu) ibgc1a      ,                                             &  
     &  neqdskg,           iminsep, imaxsep, jminsep, jmaxsep, istart  ,  &  
     &  nlhcdin,  n13    , nskip  , nspc   , nthep2 , nthep3 , nskipsf ,  &  
     &  ivpljp , ivplvc , ivplvt , icplet , icplgf , icplwf , icplpr  ,  &  
     &  icplbv , icpluv , iminn  , imaxx  , jminn  , jmaxx  , nthe    ,  &  
     &  ifrst(:50), nterm , iwall  , nlscgnu, ilim   , jlim , numsep  ,  &  
     &  iskipsf, itimets, ntpts  ,  nin    , nout   , no167a , nmiss  ,  &  
     &  neqdska  , neqdsk , nsprsi , nsprso , nsc2   , nsc3   , nenin ,  &  
     &  neqin  , neqou  , nsc1   , i95    , iwayne , kkm95  , kkm95m ,   &  
     &  kmax   , kmaxo  , icirc  , i90    , kkm90  , kkm90m , isvd    ,  &  
     &  icplxp , ndivlu , mframe , iffac , itevv  , noplot(:100),nlhfk,  &  
     &  iminy(:250)      , imaxy(:250)      , jminy(:250)      ,         &  
     &  jmaxy(:250)    , itemp  ,  irfp  , irippl , ntfcoil, iqtrubmax,  &  
     &  npitch , iripmod, neqplt , jminny(:pnx) , jmaxxy(:pnx),          &  
     &        ifnc1a

      CASE DEFAULT

      print *, "Error in buf_ini: match not found", idrss

      END SELECT

      return
      END SUBROUTINE buf_ini
