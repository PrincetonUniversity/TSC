      MODULE CLINAM
      USE PARAM

      IMPLICIT NONE
      INTEGER, PRIVATE, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
!============
! idecl:  explicitize implicit REAL declarations:
!============
!
      character*8 ifileoy,ifiloyi,ifiloyo,                               &  
     &          ifiltcv,                                                 &  
     &          ijwname,                                                 &  
     &          ikeep,ieqdska
      character*86 s100(30)
      character*100 suffix
      character*120 :: iouttsc, ilhcdou, ograph, ieqflou, ieqflin,       &  
     &                 isprsi, isprso, isprsp, ieqplt, ieqdsk, ifiledi
      INTEGER :: numargs
      CHARACTER*120 :: filename
!     -----------------------
!
!l    common / comp / ....................................... 86-04-09 .
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  nnp
!
!l    common / com1 / ....................................... 86-05-19 .
      INTEGER ::  ibgcm1      ,                                          &  
     & nx      , nz     , nxm    , nzm    , nlhcd   , nlhcd2  , iplot  ,  &  
     &  ncycle , kcycle , ires   , nimax  , nzs     ,  irayt , iglobea ,  &  
     & lujayp   , luwire , luwirer, nxp    , nzp    , ncoil  ,           &  
     &   iprnt  , iplt   , iplt2  , nskip2 , nskipl , nskipr , nis     ,  &  
     &   numsaw , npert  , neqmax , nebeam , ngrmax , itrmod , idata   ,  &  
     &   iimp   , ilte   , nmult  , ilhcd   ,                            &  
     &   ineg   , irst1  , irst2  , intact , icol   , imag   , jmag    ,  &  
     &   idtmin , jdtmin , ipest  , itpass , isym   , ifunc   , idiv   ,  &  
     &   isaw   , nwire  , nzero  , ninc   , igone  , nc0    , isurf   ,  &  
     &   ifunct , ndiv   , ipres  , idens  , nh     , nsepmax, ifk     ,  &  
     &   itot,    nnqueue, nnsize , nobs   , igroupp, nloop   ,          &  
     &   nreccdf, imovie , npts   , iplim  , nlim   , numfb  , nrecord ,  &  
     &   lrswtch,irunaway, nplrec , lenscr , icube  , iflux  , ibootst ,  &  
     &   nqflag , npsi   , npsit  , npsim  , ivplbp , ivplvi , ivplfr  ,  &  
     &   nopest , itpest , ialpha , iplate , nplate , ibalsw ,nskipsfi,  &  
     &   npsitmx,  jsym  , nslhrt , nslhpc , ifstrt , iicrh  , ifwcd  ,  &  
     &   ifncm1aa, iecrh, ieccd, idefnpe
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  icnt

      INTEGER ::  ibgcm1c
      CHARACTER*120 :: ifileou, ifilein, ifiles1, ifileen, ifile7a,      &  
     &               ifiles2,ifiles3
      CHARACTER*10 :: idate, isuffix,idumm1,idumm2,itime,imach
      CHARACTER*10, ALLOCATABLE, DIMENSION(:) ::  name
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  lmix 
      INTEGER ::  ifncm1c

      INTEGER ::  ibgcm1aa,                                              &  
     &            itrmode, icall, numpf, nreg, ncnt, itcv, nwprnt,       &  
     &            nwmax, kwprnt, iresgs , ngroup , ngroups, ngroupt,     &  
     &            nsturb, mloopx, nn, iwtype,  iltch,  numicsuf 
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  ibaloon
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  iseries, lgroup, jnreg,     &  
     &                                 nogroupt,                         &  
     &                               icmaxs, nogroup, nogroups 
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  indxd1, idelay, nfeedo,     &  
     &                                 nfeed2,                           &  
     &                              nrfb, ipext, iflagfb
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  igroupc, ilcoil
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  izsep, jzsep, izsepo,       &  
     &                                 jzsepo
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  ngrvw1, iwire, jwire,       &  
     &                                 ngrvw2, ifeed,                    &  
     &                              igroupw, ngrvw3, ngrvw4, ngrvc1,     &  
     &                              ngrvc2, ngrvc3, ngrvc4, ilwire
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  iobs, jobs, npltobs
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  ilima, jlima
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  nseg
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  iptch, nsegs
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  impnum, ichargst
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  iexvc, icoil, iexv, iexs
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  nfeedv
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  icoils
      INTEGER :: ifncm1

      INTEGER :: ibgc1a      ,                                           &  
     &   neqdskg,          iminsep, imaxsep, jminsep, jmaxsep, istart  ,  &  
     &   nlhcdin,  n13   , nskip  , nspc   , nthep2 , nthep3 , nskipsf ,  &  
     &   ivpljp , ivplvc , ivplvt , icplet , icplgf , icplwf , icplpr  ,  &  
     &   icplbv , icpluv , iminn  , imaxx  , jminn  , jmaxx  , nthe    ,  &  
     &   nterm , iwall  , nlscgnu, ilim   , jlim  , numsep  ,            &  
     &   iskipsf, itimets, ntpts  ,  nin    , nout   , no167a , nmiss  ,  &  
     &   neqdska, neqdsk , nsprsi , nsprso , nsc2   , nsc3   , nenin   ,  &  
     &   neqin  , neqou  , nsc1   , i95    , iwayne , kkm95  , kkm95m ,  &  
     &   kmax   , kmaxo  , icirc  , i90    , kkm90  , kkm90m , isvd    ,  &  
     &   icplxp , ndivlu , mframe , iffac  , itevv  , nlhfk,             &  
     &   itemp  ,  irfp  , irippl , ntfcoil, iqtrubmax,                  &  
     &   npitch , iripmod, neqplt 
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  ifrst
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  noplot
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  iminy, imaxy, jminy, jmaxy
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  jminny, jmaxxy
      INTEGER :: ifnc1a
!
!l    common / com2 / ....................................... 86-05-19 .
      REAL*8 ::  begcm2     ,                                            &  
     &  bsitot , cditot , alx    , alz    , pi      , pmmin  , pmmax  ,  &  
     &   vol   , uint  , sqrt2  , facimp , sumdb  , sumb   , sf      ,   &  
     &   tpi,    tlastg, aminor,  dtfac,                                 &  
     &   dtmax  , amu    , eta    , dtsf   , time   , dt     , reboun  ,  &  
     &   amux   , dtmin  , dtold  , ekin   , ekino   ,                   &  
     &   renum  , snum   , amach  , amuz   , erate  , qadia  , apldia  ,  &  
     &   ener3  , ener4  , ener5  , ener6  , ener7  , rlinv  , qprof2p ,  &  
     &   eintp  , eintt  , ekinp  , ekint  , enmjs  , usdd   , usdt    ,  &  
     &   usdv   , usde   , usdi   , usdp   , usdr   , zmag   , p0s     ,  &  
     &   tzerfbs, tzerfb , udsd   , udst   , udsv   , udse   , udsi    ,  &  
     &   udsp   , udsr   , facinc , gmin   , pplash , pplaso , udsh    ,  &  
     &   usdh   , dtmins , tmaxs  , psi2b  , psiob  , psib   , dpsi    ,  &  
     &  bpowers , erates , psidif , etah   , ai2    , ai3    , helic   ,  &  
     & crashtime, hypermult, transmult, hyperfrac, teflat_time
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  acoef , redef  
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  global
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  sawtime
      REAL*8 :: fincm2aa

      REAL*8 ::  begc2aa,                                                &  
     &   ai4    , apl    , xmag   , psimin , ts     , alp    , enermjo ,  &  
     &   w1lim  ,alhitot , aqp    , aip    , heact  , beta   , eqrate  ,  &
     &   tped   ,qadd    , fhmodei, pwidthc, chiped , fraci(8),           &
     &   nflag  ,expn1   , expn2  , firitb , secitb , fracn0  , newden ,  &
     &   psirat , anoflux, fracox , fracca , enermj , toldfv , wdotmw  ,  &  
     &   ain99  , ain95  , ain90  , vtran  , freqlh , fraclst ,          &  
     &  psilimo ,          zmag2  , cwire  , tdel   , tdel2  , tdelw   ,  &  
     &            plimiter,eqrat2 , vsec   , curmin , curmax , ener9   ,  &  
     &   ffac   , sum1e  , sum2e  , sum3e  , sum4e  , psaboun, resis   ,  &  
     &   tauekg , torl   , r1mix  , smallr , tfmult , betapol ,          &  
     &   ali2   , aliga  , apn    , apc    , atc    , smallt , qboun   ,  &  
     &   ene2min, ene2max, th, thm, etamin , volts  , curs   , powers  ,  &  
     &   fwitot, ecitot, fluxlink,fincm2
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  cisum

      REAL*8 ::  begc2a      ,                                           &  
     &   tcuodtp, deex   , deez   , dxisq  ,                             &  
     &   xcurf  , zcurf  , dxdz , tcurdtp, xplas  , rebounn ,            &  
     &   rebouno, tfmax  , dzisq  , phirpg , rdpsi  , rdpsis , rdthe   ,  &  
     &   rdthes , dthe   , hfluxp , diatot , teform , tiform , chiauxs ,  &  
     &   chiohms, arad   , tevmax , rfbval , pfbval , apld0  , phalo   ,  &  
     &   psilim , alphap , p0     , tcuro  , alphag , gzero  , delpsi  ,  &  
     &   gp1    , gp2    , qzero  , q0     , ai3d0  , ai3o   , fbchi   ,  &  
     &   tcurdtpcl,cplasfs, p0fs  , aplo   , ai6    , tcuros , dshell  ,  &  
     &   alfeed , blfeed , x1sep  , x2sep  , z1sep  , z2sep  , tzero   ,  &  
     &   xdist  , xlim   , zlim   , r0     , alphar , e0     , alphae  ,  &  
     &   tpwr   , cpwr   , psi95  , tjphi  , dtjphi , hfluxpo, hfluxav ,  &  
     &  ebsav(9), amgas  , zgas   , ajpress, bpowsmw, thalos , whalos ,  &  
     &    pfwmc ,  plhmc ,  zeffb , uintold,  fsaw,   tfluxs , fbchii
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  cmomo, cmom
      REAL*8 ::   finc2a

      REAL*8 ::  begc2b      ,                                           &  
     &                     pcurmin, pcurmax, tcurmin, tcurmax, prad    ,  &  
     &                     arps   , acps   , avps   , ait    , q95pct  ,  &  
     &   small  , dpar   , alts   , arts   , acts   , avts   , ucor    ,  &  
     &   etav   , tevv   , psisep , times  , dts    , xlim2  , abeam   ,  &  
     &   dbeam  , zplas  , dtmaxs , enerst2, betar  , qsaw   , resgap  ,  &  
     &   xzeric , axic   , zzeric , bzic   , tauems , zimp   , phi     ,  &  
     &   phi2   , eshelw , eratew , zeff   , ccon   , ebounx , fasym   ,  &  
     &   vzero  , pzero  ,  tolds , tolds2 , vsec2  , pold   , pold2   ,  &  
     &   vloopp , delg   , gprfp  , betaj  , paux   , pohmic , palpha  ,  &  
     &   xmagw  , zmagw  , psis1  , psis2  , psis3  , ellmom , palphap ,  &  
     &  xsepcal , zsepcal, psepcal, gvmin  , gvmax  , templmn , templmx,  &  
     &   gcmin  , gcmax  , poh0   , dt1    , dt2    , dt3    , zzero   ,  &  
     &   betat0 , betat  , rvolav , el90   , del90  , rtot   , psi90  ,  &  
     &   tav1   , tav2   , tav3   , tav4   , tav5   , sum7e  , enerpb ,  &  
     &   timedos, tmovie , dtmovie, gvmin0 , gvmax0 , powerr , poweri ,  &  
     &   energyi, energyr, resid  , rzerw  , azerw  , ezerw  , dzerw  ,  &  
     &   tpergl , tave1  , tave2  , tave3  , tave4  , tave5  , deltmks,  &  
     &   sum1o  , sum2o  , sum3o  , sum4o  , sum7o  , sumdo  ,           &  
     &                finc2b
!
      REAL*8 ::  begc2c,                                                 &  
     &   ebeamkev , ambeam , fracpar ,                                   &  
     &   rvolavi, rvolavh,alphades,alphafrs,alphabes, toldaux,           &  
     &   gpmin  , gpmax  , gemin  , gemax  , gpmint , gpmaxt , gemint  ,  &  
     &   gemaxt , ripmax , rtfcoil,                                      &  
     &                rmajor,rminor,shape3,                              &  
     &                shape4,shape5,shape6,shape7
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xsep, zsep, psep
      REAL*8, ALLOCATABLE, DIMENSION(:) ::   engrupi, gvsuma, delv,      &  
     &                                 engroup
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  globmin, globmax
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  xplot, zplot
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  pltsav
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  apld0sv
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  fluxu, fluxl
      REAL*8 ::                                                          &  
     &                hfluxtot, el95   , del95  , pscrape ,              &  
     &                eps3co , eps4co ,                                  &  
     &                fraclos,  fraclosr, ripmult ,                      &  
     &                alh,dlh,a1lh,a2lh,                                 &  
     &   zion   , aion   , cprof  ,                                      &  
     &                finc2c
!
!l    common / com3 / ....................................... 86-05-17 .
      REAL*8 ::  begcm3
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xary, xsqoj, ajey, xarh
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  zary, vzy
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  w11, w12
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  psi, ajphi, vd, gs, g, q,   &  
     &                                 psio,                             &  
     &                                psizer, rjcd, rjcd2
      REAL*8, ALLOCATABLE, DIMENSION(:) ::   cgroup,sgroup
      REAL*8 ::  fincm3

      REAL*8 ::  begc3a
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::                             &  
     &                              wo, go, uo, abig, rsurf, b, bo,      &  
     &                              uary, omeg, etay, bsqsv, rdens,      &  
     &                              rjdb
      REAL*8 ::  finc3a

      REAL*8 ::  begc3b
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  w, u, qo, r, ro, qsurf,    &  
     &                                 pr, pro,                          &  
     &                              roj
      REAL*8 ::  finc3b

      REAL*8 ::  begc3c
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  go2, pso2, uo2, bmagy,     &  
     &                                 ptsave,                           &  
     &                              gtsave
      REAL*8 ::  finc3c

      REAL*8 ::  begc3d
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  psibl, psibr, psitl, psitr,   &  
     &                                 psisl,                            &  
     &                           psisr, psiblo, psibro, psiobl,          &  
     &                           psiobr, psidbl, psidbr, omderr,         &  
     &                           omderl, avbndl, avbndr
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  psibb, psibt, psitb, psitt,   &  
     &                                 psisb,                            &  
     &                           psist, psibbo, psibto, psiobb, psiobt,  &  
     &                           psidbb, psidbt, omderb, omdert,         &  
     &                           avbndb, avbndt
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  vn, psiold
      REAL*8 ::  finc3d
!
!l    common / com4 / ....................................... 86-05-18 .
      REAL*8 ::  begcm4
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  gcurfb , gcurfbo, gvgpmax,   &  
     &                                 gvgpmin,                          &  
     &                              gpgpmax, gpgpmin, gigpmax, gigpmin,  &  
     &                              gvolt0
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  fsave1
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  sumfeed, savfeed
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  xcon0, zcon0
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  fadj
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  amats
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  bvecs
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  shpsum, shpdif, shpdifo
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xveh
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  wmega, phin
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  sigmasvd
      REAL*8 :: g1old, g2old, g3old, g4old, g5old, g6old, ckstrt,        &  
     &          fincm4

      REAL*8 ::  begc4a 
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  fluxu0, fluxl0, xobs, zobs
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  tfbons, tfbon, tfbofs,       &  
     &                                 tfbof,                            &  
     &                             fbfac , fbcon, fbfac1, tzerfba,       &  
     &                             fbfacd, tperfba, fbfaci 
      REAL*8 ::  finc4a

      REAL*8 ::  begc4b
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  gcur
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  pcur, ppres, tpro, tpros,    &  
     &                                 gzerov,                           &  
     &                            beamp,rnorm, fact, fbchia,facd  ,      &  
     &                            vloopv, xmagz, zmagz, tevv0, ffac0,    &  
     &                            rzerv, azerv, zeffv, fbchiia
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  resgs     
      REAL*8, ALLOCATABLE, DIMENSION(:) ::   tpest 
      REAL*8 ::  finc4ba

      REAL*8 :: begc4ba
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  ezerv, dzerv, plhamp,        &  
     &                                 betarv,                           &  
     &                            alpharv, frcparv, alhd,  dlhd,         &  
     &                            a1lhd, a2lhd, aclhd,  dclhd,           &  
     &                            a1clhd, a2clhd, afwd, dfwd,            &  
     &                            a1fwd, a2fwd, acfwd, dcfwd,            &  
     &                            a1cfwd, a2cfwd, thalov, whalov,        &  
     &                            picrh, fwcd, heactv, tpedv, qaddv,    &
     &                            fhmodeiv, pwidthcv, chipedv,           &
     &                            nflagv,expn1v,expn2v,firitbv,          &
     &                            secitbv, fracn0v, newdenv,             &
     &                            pecrh, eccd, aecd, decd, a1ecd, a2ecd
      REAL*8 :: finc4b

      REAL*8 :: begc4c
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  atnvw1, atnvw2, atnvw3,    &  
     &                                 atnvw4,                           &  
     &                      atnvc1, atnvc2, atnvc3, atnvc4,fraciv 
      REAL*8 :: finc4c
!
!l    common / com5 / ....................................... 86-05-16 .
      REAL*8 :: begcm5
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xcoil, zcoil, ccoil, dcoil,   &  
     &                                 ccoils,                           &  
     &                             aturnsc, rscoil, rscoils, aindc,      &  
     &                             ccoil0 , worksp 
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xwire, zwire, cwics, rswire,  &  
     &                                  rswires,                         &  
     &                             cwires, cwire0, cwiradi, ajwire,      &  
     &                             aturnsw, aindw, resave
      REAL*8, ALLOCATABLE, DIMENSION(:) ::   xlima, zlima
      REAL*8 :: fincm5
!
      REAL*8 ::  begc5a
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  as0 ,as1 ,as2 ,as3 ,as4 ,    &  
     &                           cs0 ,cs1 ,cs2 ,cs3 ,cs4 ,               &  
     &                           dsi0,dsi1,dsi2,dsi3,dsi4,               &  
     &                           dse0,dse1,dse2,dse3,dse4,               &  
     &                           gam ,as  ,cs  ,dsi ,dse ,               &  
     &                           vary,anue,vp  ,vp2 ,vpg ,               &  
     &                           xmja,adpo,adeo,adno,adio,               &  
     &                           adi ,adp ,ade ,adn ,adpi,               &  
     &                           ado ,adoo,d2pi,d2pe,d2s ,               &
     &                           qamin,bsqar,bsqiar,bpolar,              &  
     &                           flux0,flux10,flux20,flux30,             &  
     &                           chiitima,chietema,diffnema,             &  
     &                           chiitimao,chietemao,diffnemao,          &  
     &                           dchii_dtip,dchii_dtep,                  &
     &                           dchie_dtip,dchie_dtep,                  &
     &                           dti_dpsi,dte_dpsi,                      &
     &                           velthi11, velthi21, velthi31,           &  
     &                           ratioi,ratioe 
      REAL*8 :: finc5a

      REAL*8 :: begc5b
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  equila, etpara, deltaa,      &  
     &                                 gxmja, gxmja2,                    &  
     &                           xmja2,  savei , savee,  srave, polflx,  &  
     &                           torflx, qprof,  temflx, volflx,tprof,   &  
     &                           gprof , gemflx, xsv2 ,  xsv   ,ftrap,   &  
     &                           savia , avhyp,  sraveb, vlooph,alphanu,  &  
     &                           sraveedg, dperpa, aplosts, aplostsr,    &  
     &                           savea,  savebre, savecyc, saveimp,      &  
     &                           savefw, anhe,   savebm, savibm, savelh,  &  
     &                           savilh, rimp1,  chiesec, alphade,       &  
     &                           alphapr, sighot, chiisec, anheb,ajaveq,  &  
     &                           densum, gja2,   rminora, rmajora,       &  
     &                           elonga, ajave, ajbtsq, ajbsq, savifw
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  qprof2
      REAL*8 :: finc5b

      REAL*8 :: begc5c
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  asv, bsv, csv, dsv
      REAL*8 :: finc5c
!
!
      REAL*8 ::  begcap
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  vgain  , voltmax , tramp   ,  &  
     &                             voltn  , diold   , vgain2  ,          &  
     &                             rewire , diol2   , vgain3  ,          &  
     &                             thrumx , thrtp   , thrtm   ,          &  
     &                             thrur  , thrps   , vfeed   ,          &  
     &                             save1  , save2   , save3   ,          &  
     &                             save4  , save5   , save6   ,          &  
     &                             thruc   
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  x2ave
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  cwirin 
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  cwold, gcursv, cdot
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  gmary  
      REAL*8 :: dtcorr ,  rctim  ,thrp   
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xrplate, xlplate, zrplate,   &  
     &                                 zlplate
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  strike, fplate, dsep
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  hplate, dplate
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  xsega, zsega
      REAL*8 :: fincap
!
!
      REAL*8 :: begc6a
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  rsceg  , ceg  , cegp1  ,     &  
     &                                 cegm1  ,                          &  
     &                              ceg0 ,  ceg0m1 ,  aintceg,           &  
     &                              flxeg,  flxegm1,  flxegm2,           &  
     &                              veg  ,  vegm1  ,  aturneg,           &  
     &                              voheg,  vohegm1,  ceg0m2,            &  
     &                              teineg ,tauseg ,                     &  
     &                              gainpeg,gaindeg,gainieg,             &  
     &                              vmineg,vmaxeg,                       &  
     &                              dcur   ,dcurm1 , dcurm2 
      REAL*8, ALLOCATABLE, DIMENSION(:) ::   vegplot
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  gvolt
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  ameg, ameginv
      REAL*8 :: dtsave,                                                  &  
     &          finc6a
!
      REAL*8 :: begcot,                                                  &  
     &          tlast 
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  tempc, dxcoil, dzcoil, fcu,   &  
     &                                 fss,                              &  
     &                             ccics
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  templ
      REAL*8 :: fincot
!
      REAL*8 :: begufc
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  ufdata
      REAL*8 :: finufc
  
      REAL*8 ::  delpmx,delgmx
!
!     added for quasineutrality
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  nbeami
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  q_snbi
      INTEGER :: nspec_beam
!
!     ------------------------------------------------------------------
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE CLINAM
