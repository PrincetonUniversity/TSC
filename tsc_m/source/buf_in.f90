      SUBROUTINE buf_in( lu, idrss, idrse )
      USE CLINAM
      USE SAPROP
      USE NONCOR
      USE SPECIE
      USE RADTAB
      USE FVV1
      USE TCVCOM
      USE RUNAWAY
      USE WALLCL

      IMPLICIT NONE
      INTEGER :: lu, i, j, k
      CHARACTER*(*) :: idrss, idrse
!     CHARACTER*8 :: idrss, idrse
      CHARACTER*8 :: recnam

!     write(6,*) "reading: ", trim(idrss)

      recnam=trim(idrss)
!     SELECT CASE (trim(idrss)) 
      SELECT CASE (trim(recnam)) 

      CASE ("begcm2")

      read(lu)  begcm2     ,                                              &  
     &  bsitot , cditot , alx  , alz  , pi   , pmmin  , pmmax  ,          &  
     &   vol , uint , sqrt2 , facimp , sumdb , sumb  , sf      ,          &  
     &   tpi,tlastg,aminor,                                               &  
     &   acoef(:4999+pncoil),redef(:4999+pncoil),dtfac,                   &  
     &   dtmax , amu  , eta  , dtsf , time , dt   , reboun  ,             &  
     &   amux , global(:pglob) , dtmin , dtold , ekin , ekino  ,          &  
     &   renum , snum , amach , amuz , erate , qadia , apldia  ,          &  
     &   ener3, ener4 , ener5  , ener6  , ener7  , rlinv  , qprof2p ,     &  
     &   eintp  , eintt  , ekinp  , ekint  , enmjs  , usdd   , usdt    ,  &  
     &   usdv   , usde   , usdi   , usdp   , usdr   , zmag   , p0s     ,  &  
     &   tzerfbs, tzerfb , udsd   , udst   , udsv   , udse   , udsi    ,  &  
     &   udsp   , udsr   , facinc , gmin   , pplash , pplaso , udsh    ,  &  
     &   usdh   , dtmins , tmaxs  , psi2b  , psiob  , psib   , dpsi    ,  &  
     &  bpowers , erates , psidif , etah   , ai2    , ai3    , helic   ,  &  
     & crashtime, sawtime(:psaw), hypermult, transmult, hyperfrac,        &  
     &  fincm2aa
 

      CASE ("begc2aa")

      read(lu) begc2aa,                                                   &  
     &   ai4    , apl    , xmag   , psimin , ts     , alp    , enermjo ,  &  
     &   w1lim  ,alhitot , aqp    , aip    , heact  , beta   , eqrate  ,  &  
     &   psirat , anoflux, fracox , fracca , enermj , toldfv , wdotmw  ,  &  
     &   ain99 , ain95  , ain90  , vtran  ,cisum(:5), freqlh , fraclst ,  &  
     &  psilimo ,          zmag2  , cwire  , tdel   , tdel2  , tdelw   ,  &  
     &            plimiter,eqrat2 , vsec   , curmin , curmax , ener9   ,  &  
     &   ffac   , sum1e  , sum2e  , sum3e  , sum4e  , psaboun, resis   ,  &  
     &   tauekg , torl   , r1mix  , smallr , tfmult , betapol ,           &  
     &   ali2   , aliga  , apn    , apc    , atc    , smallt , qboun   ,  &  
     &   ene2min, ene2max, th, thm, etamin , volts  , curs   , powers  ,  &  
     &   fwitot , ecitot, fluxlink,tped,  qadd,  fhmodei, pwidthc, chiped, fraci(:8), &
     &   nflag, expn1, expn2, firitb, secitb, fracn0, newden,      fincm2

      CASE ("begc2a")

      read(lu)  begc2a      ,                                             &  
     &   tcuodtp, cmomo(:2,:2)      , deex   , deez   , dxisq  ,          &  
     &   xcurf , zcurf , cmom(:2,:2) ,dxdz , tcurdtp, xplas  , rebounn ,  &  
     &   rebouno, tfmax  , dzisq  , phirpg , rdpsi  , rdpsis , rdthe   ,  &  
     &   rdthes , dthe   , hfluxp , diatot , teform , tiform , chiauxs ,  &  
     &   chiohms, arad   , tevmax , rfbval , pfbval , apld0  , phalo   ,  &  
     &   psilim , alphap , p0     , tcuro  , alphag , gzero  , delpsi  ,  &  
     &   gp1    , gp2    , qzero  , q0     , ai3d0  , ai3o   , fbchi   ,  &  
     &   tcurdtpcl, cplasfs, p0fs   , aplo , ai6   , tcuros , dshell  ,   &  
     &   alfeed , blfeed , x1sep  , x2sep  , z1sep  , z2sep  , tzero   ,  &  
     &   xdist  , xlim   , zlim   , r0     , alphar , e0     , alphae  ,  &  
     &   tpwr   , cpwr   , psi95  , tjphi  , dtjphi , hfluxpo, hfluxav ,  &  
     &  ebsav(:9), amgas  , zgas   , ajpress, bpowsmw, thalos , whalos ,  &  
     &    pfwmc ,  plhmc ,  zeffb , uintold,  fsaw,   tfluxs, fbchii,     &  
     &                finc2a

      CASE ("begc2b")

      read(lu) begc2b      ,                                              &  
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
     &   betat0 , betat  , rvolav , el90   , del90  , rtot   , psi90  ,   &  
     &   tav1   , tav2   , tav3   , tav4   , tav5   , sum7e  , enerpb ,   &  
     &   timedos, tmovie , dtmovie, gvmin0 , gvmax0 , powerr , poweri ,   &  
     &   energyi, energyr, resid  , rzerw  , azerw  , ezerw  , dzerw  ,   &  
     &   tpergl , tave1  , tave2  , tave3  , tave4  , tave5  , deltmks,   &  
     &   sum1o  , sum2o  , sum3o  , sum4o  , sum7o  , sumdo  ,            &  
     &                finc2b

      CASE ("begc2c")
      
      read(lu)  begc2c,                                                   &  
     &   ebeamkev , ambeam , fracpar ,                                    &  
     &   rvolavi, rvolavh,alphades,alphafrs,alphabes, toldaux,            &  
     &   gpmin  , gpmax  , gemin  , gemax  , gpmint , gpmaxt , gemint  ,  &  
     &   gemaxt , ripmax , rtfcoil,                                       &  
     &                rmajor,rminor,shape3,                               &  
     &                shape4,shape5,shape6,shape7,                        &  
     & xsep(:pnsep),zsep(:pnsep),psep(:pnsep),                            &  
     &  engrupi(:pngroup),                                                &  
     &   xplot(:6,:2*pnparts),zplot(:6,:2*pnparts), engroup(:pngroup) ,   &  
     &                gvsuma(:pngroup), delv(:pngroup) ,                  &  
     &                globmin(:pgls),globmax(:pgls),                      &  
     &                pltsav(:2*pnsave*pglobp),apld0sv(:2*pnsave),            &
     &                fluxu (:pobs,:2*pnsave),   fluxl (:pobs,:2*pnsave),     &
     &                hfluxtot, el95   , del95  , pscrape ,               &  
     &                eps3co , eps4co ,                                   &  
     &                fraclos,  fraclosr, ripmult ,                       &  
     &                alh,dlh,a1lh,a2lh,                                  &  
     &   zion   , aion   , cprof  ,                                       &  
     &                finc2c

      CASE ("begcm3") 

      read(lu) begcm3,                                                   &  
     &         xary (:penx),     zary      (:penz), vzy     (:penz),     &  
     &         xsqoj(:penx),     ajey (:penx)     , xarh(:penx)    ,     &  
     &         w11(:penx+2*penz),     w12(:penx+2*penz)    ,             &  
     &         psi (:penx,:penz),ajphi(:penx,:penz),                     &  
     &         vd   (:penx,:penz),                                       &  
     &         gs   (:penx,:penz),g    (:penx,:penz), q  (:penx,:penz),  &  
     &         psio(:penx,:penz),psizer(:penx,:penz),rjcd(:penx,:penz),  &  
     &         rjcd2(:penx,:penz),                                       &  
     &   cgroup(:pngroup),sgroup(:pngroup),                              &  
     &                fincm3

      CASE ("begc3a")

      read(lu) begc3a,                                                   &  
     &         wo    (:pnx,:pnz), go    (:pnx,:pnz),                     &  
     &         uo    (:pnx,:pnz), abig  (:pnx,:pnz), rsurf (:pnx,:pnz),  &  
     &         b     (:pnx,:pnz), bo   (:pnx,:pnz), uary(:pnx,:pnz),     &  
     &         omeg  (:pnx,:pnz), etay  (:pnx,:pnz),                     &  
     &         bsqsv (:pnx,:pnz), rdens (:pnx,:pnz), rjdb(:pnx,:pnz),    &  
     &                finc3a

      CASE ("begc3b")

      read(lu) begc3b,                                                   &  
     &         w     (:pnx,:pnz), u     (:pnx,:pnz), qo    (:pnx,:pnz),  &  
     &         r     (:pnx,:pnz), ro    (:pnx,:pnz), qsurf (:pnx,:pnz),  &  
     &         pr    (:pnx,:pnz), pro   (:pnx,:pnz), roj   (:pnx,:pnz),  &  
     &                finc3b

      CASE ("begc3c")

      read(lu)  begc3c,                                                   &  
     &          go2   (:pnx,:pnz), pso2  (:pnx,:pnz), uo2   (:pnx,:pnz),  &  
     &          bmagy (:pnx,:pnz), ptsave(:pnx,:pnz), gtsave(:pnx,:pnz),  &  
     &                finc3c

      CASE ("begc3d")

      read(lu) begc3d,                                                    &  
     &      psibl(:penz)  ,psibr(:penz)  ,psibb(:penx)  ,psibt(:penx)  ,  &  
     &      psitl(:penz)  ,psitr(:penz)  ,psitb(:penx)  ,psitt(:penx)  ,  &  
     &      psisl(:penz)  ,psisr(:penz)  ,psisb(:penx)  ,psist(:penx)  ,  &  
     &      psiblo(:penz) ,psibro(:penz) ,psibbo(:penx) ,psibto(:penx) ,  &  
     &      psiobl(:penz) ,psiobr(:penz) ,psiobb(:penx) ,psiobt(:penx) ,  &  
     &      psidbl(:penz) ,psidbr(:penz) ,psidbb(:penx) ,psidbt(:penx) ,  &  
     &      omderr(:penz) ,omderl(:penz) ,omderb(:penx) ,omdert(:penx) ,  &  
     &      avbndl(:penz) ,avbndr(:penz) ,avbndb(:penx) ,avbndt(:penx) ,  &  
     &      vn(:pboun)    ,psiold(:pboun),                                &  
     &                finc3d

      CASE ("begcm4")

      read(lu) begcm4,                                                    &  
     &                gcurfb (:pngroup),       gcurfbo(:pngroup),         &  
     &                gvgpmax(:pngroup),       gvgpmin(:pngroup),         &  
     &                gpgpmax(:pngroup),       gpgpmin(:pngroup),         &  
     &                gigpmax(:pngroup),       gigpmin(:pngroup),         &  
     &                     gvolt0(:pngroup),                              &  
     &                fsave1(:pnfeed,:pdelay),                            &  
     &                sumfeed(:pnfeed),savfeed(:pnfeed),                  &  
     &                 xcon0(:ptpts,:10), zcon0(:ptpts,:10),              &  
     &   fadj(:pngroup,:pncoil),amats(:20,:20,:pdelay),                   &  
     &   bvecs(:20,:pdelay),shpsum(:20),shpdif(:20),shpdifo(:20),         &  
     &   xveh(:201),g1old,g2old,g3old,g4old,g5old,g6old,                  &  
     &    wmega(:200), phin(:200), ckstrt,sigmasvd(:pngroup),             &  
     &                fincm4

      CASE ("begc4a")
      read(lu) begc4a,                                                    &  
     &                fluxu0(:pobs),          fluxl0(:pobs),              &  
     &                xobs  (:pobs),           zobs   (:pobs),            &  
     &                tfbons(:pnfeed),           tfbon  (:pnfeed),        &  
     &                tfbofs(:pnfeed),           tfbof  (:pnfeed),        &  
     &                fbfac (:pnfeed),           fbcon  (:pnfeed),        &  
     &                fbfac1(:pnfeed),           tzerfba(:pnfeed),        &  
     &                fbfacd(:pnfeed),           tperfba(:pnfeed),        &  
     &                fbfaci(:pnfeed),                                    &  
     &                finc4a

      CASE ("begc4b")
      read(lu) begc4b,                                                    &  
     &                gcur  (:ptpts,:pngroup),                            &  
     &                pcur  (:ptpts),           ppres (:ptpts),           &  
     &                tpro  (:ptpts),           tpros (:ptpts),           &  
     &                gzerov(:ptpts),           beamp (:ptpts),           &  
     &                rnorm (:ptpts),           fact  (:ptpts),           &  
     &                fbchia(:ptpts),           fbchiia(:ptpts),          &
     &                                          facd  (:ptpts),           &  
     &                vloopv(:ptpts),           tpest (:100)  ,           &  
     &                xmagz (:ptpts),           zmagz (:ptpts),           &  
     &                tevv0 (:ptpts),           ffac0 (:ptpts),           &  
     &                resgs (:pngroup),         zeffv(:ptpts),            &  
     &                rzerv(:ptpts),            azerv(:ptpts),            &  
     &                finc4ba

      CASE ("begc4ba")
      read(lu) begc4ba,                                                   &  
     &                ezerv(:ptpts),            dzerv(:ptpts),            &  
     &                plhamp(:ptpts),           betarv(:ptpts),           &  
     &                alpharv(:ptpts),          frcparv(:ptpts),          &  
     &                alhd(:ptpts),             dlhd(:ptpts),             &  
     &                a1lhd(:ptpts),            a2lhd(:ptpts),            &  
     &                aclhd(:ptpts),            dclhd(:ptpts),            &  
     &                a1clhd(:ptpts),           a2clhd(:ptpts),           &  
     &                afwd(:ptpts),             dfwd(:ptpts),             &  
     &                a1fwd(:ptpts),            a2fwd(:ptpts),            &  
     &                acfwd(:ptpts),            dcfwd(:ptpts),            &  
     &                a1cfwd(:ptpts),           a2cfwd(:ptpts),           &  
     &                thalov(:ptpts),           whalov(:ptpts),           &  
     &                picrh(:ptpts),            fwcd(:ptpts),            &
     &                heactv(:ptpts),           tpedv(:ptpts),            &
     &                qaddv(:ptpts),            fhmodeiv(:ptpts),         &
     &                pwidthcv(:ptpts),         chipedv(:ptpts),          &
     &                fraciv(:8,:ptpts),        nflagv(:ptpts),           &
     &                expn1v(:ptpts),           expn2v(:ptpts),           &
     &                firitbv(:ptpts),          secitbv(:ptpts),          &
     &                fracn0v(:ptpts),          newdenv(:ptpts), finc4b,  &
     &                pecrh(:ptpts),            eccd(:ptpts),             &
     &                aecd(:ptpts),             decd(:ptpts),             &
     &                a1ecd(:ptpts),            a2ecd(:ptpts)


      CASE ("begc4c")
 
      read(lu) begc4c,                                                    &  
     &          atnvw1(:ptpts,:pnwire),                                   &  
     &           atnvw2(:ptpts,:pnwire),                                  &  
     &          atnvw3(:ptpts,:pnwire),                                   &  
     &          atnvw4(:ptpts,:pnwire),                                   &  
     &         atnvc1(:ptpts,:pncoil),                                    &  
     &         atnvc2(:ptpts,:pncoil),                                    &  
     &          atnvc3(:ptpts,:pncoil),                                   &  
     &        atnvc4(:ptpts,:pncoil),                                     &  
     &                finc4c

      CASE ("begcm5")
      read(lu) begcm5,                                                    &  
     &         xcoil  (:pncoil), zcoil  (:pncoil), ccoil  (:pncoil),      &  
     &         dcoil  (:pncoil), ccoils (:pncoil),                        &  
     &         aturnsc(:pncoil),  rscoil (:pncoil),                       &  
     &         rscoils(:pncoil), aindc  (:pncoil), ccoil0 (:pncoil),      &  
     &         xwire  (:pnwire), zwire  (:pnwire),                        &  
     &         cwics  (:pnwire), rswire (:pnwire), rswires(:pnwire),      &  
     &         cwires (:pnwire), cwire0 (:pnwire), cwiradi(:pnwire),      &  
     &         ajwire (:pnwire),                                          &  
     &         aturnsw(:pnwire),                                          &  
     &         aindw  (:pnwire), resave(:pnwire),                         &  
     &         xlima  (:pnlim) ,   zlima(:pnlim) ,                        &  
     &         worksp(:pncoil) ,                                          &  
     &                fincm5

      CASE ("begc5a")
      read(lu) begc5a,                                                    &  
     &     as0 (:ppsi),as1 (:ppsi),as2 (:ppsi),as3 (:ppsi),as4 (:ppsi),   &  
     &     cs0 (:ppsi),cs1 (:ppsi),cs2 (:ppsi),cs3 (:ppsi),cs4 (:ppsi),   &  
     &     dsi0(:ppsi),dsi1(:ppsi),dsi2(:ppsi),dsi3(:ppsi),dsi4(:ppsi),   &  
     &     dse0(:ppsi),dse1(:ppsi),dse2(:ppsi),dse3(:ppsi),dse4(:ppsi),   &  
     &     gam (:ppsi),as  (:ppsi),cs  (:ppsi),dsi (:ppsi),dse (:ppsi),   &  
     &     vary(:ppsi),anue(:ppsi),vp  (:ppsi),vp2 (:ppsi),vpg (:ppsi),   &  
     &     xmja(:ppsi),adpo(:ppsi),adeo(:ppsi),adno(:ppsi),adio(:ppsi),   &  
     &     adi (:ppsi),adp (:ppsi),ade (:ppsi),adn (:ppsi),adpi(:ppsi),   &  
     &     qamin(:ppsi),bsqar(:ppsi),bsqiar(:ppsi),bpolar(:ppsi),         &  
     &     flux0(:ppsi),flux10(:ppsi),flux20(:ppsi),flux30(:ppsi),        &  
     &     chiitima(:ppsi),chietema(:ppsi),diffnema(:ppsi),               &  
     &     chiitimao(:ppsi),chietemao(:ppsi),diffnemao(:ppsi),            &  
     &     dchii_dtip(:ppsi),dchii_dtep(:ppsi),dchie_dtip(:ppsi),         &
     &     dchie_dtep(:ppsi),dti_dpsi(:ppsi),dte_dpsi(:ppsi),             &
     &     velthi11(:ppsi), velthi21(:ppsi), velthi31(:ppsi),             &  
     &     ratioi(:ppsi),ratioe(:ppsi),                                   &  
     &              finc5a

      CASE ("begc5b" )
      read(lu) begc5b,                                                    &  
     &      equila(:ppsi),       etpara(:ppsi),       deltaa(:ppsi),      &  
     &      gxmja (:ppsi),       gxmja2(:ppsi),        xmja2(:ppsi),      &  
     &      savei (:ppsi),        savee(:ppsi),        srave(:ppsi),      &  
     &      polflx(:ppsi),       torflx(:ppsi),        qprof(:ppsi),      &  
     &      temflx(:ppsi),       volflx(:ppsi),        tprof(:ppsi),      &  
     &      gprof (:ppsi),       gemflx(:ppsi),        xsv2 (:ppsi),      &  
     &      xsv   (:ppsi),        ftrap(:ppsi),       qprof2(:ppsi+1),    &  
     &      savia (:ppsi),        avhyp(:ppsi),     sraveb(:ppsi),        &  
     &      vlooph(:ppsi),      alphanu(:ppsi),      sraveedg(:ppsi),     &  
     &      dperpa(:ppsi),     aplosts(:ppsi),     aplostsr(:ppsi),       &  
     &      savea(:ppsi),       savebre(:ppsi),    savecyc(:ppsi),        &  
     &      saveimp(:ppsi),     savefw(:ppsi),     savifw(:ppsi),         &
     &                                             anhe(:ppsi),           &  
     &      savebm(:ppsi),      savibm(:ppsi),     savelh(:ppsi),         &  
     &      savilh(:ppsi),      rimp1(:ppsi),    chiesec(:ppsi),          &  
     &      alphade(:ppsi),     alphapr(:ppsi),   sighot(:ppsi),          &  
     &      chiisec(:ppsi),      anheb(:ppsi),    ajaveq(:ppsi),          &  
     &      densum(:ppsi),    gja2(:ppsi),                                &  
     &      rminora(:ppsi),   rmajora(:ppsi),    elonga(:ppsi),           &  
     &      ajave(:ppsi), ajbtsq(:ppsi), ajbsq(:ppsi),                    &  
     &                finc5b

      CASE ("begc5c")
      read(lu)  begc5c,                                                   &  
     &           asv   (:ppsi,:pneq),     bsv(:ppsi,:pneq),               &  
     &           csv   (:ppsi,:pneq),     dsv(:ppsi,:pneq),               &  
     &                finc5c

      CASE ("begcap")
      read(lu) begcap,                                                    &  
     &           vgain (:pncoil) , voltmax(:pncoil) , tramp (:pncoil)  ,  &  
     &           voltn (:pncoil) , diold  (:pncoil) , vgain2(:pncoil)  ,  &  
     &           rewire(:pncoil) , diol2  (:pncoil) , vgain3(:pncoil)  ,  &  
     &           thrumx(:pncoil) , thrtp  (:pncoil) , thrtm (:pncoil)  ,  &  
     &           thrur (:pncoil) , thrps  (:pncoil) , vfeed (:pncoil)  ,  &  
     &           thruc (:pncoil) , x2ave (:ppsi),                         &  
     &           save1 (:pncoil) , save2  (:pncoil) , save3 (:pncoil)  ,  &  
     &           save4 (:pncoil) , save5  (:pncoil) , save6 (:pncoil)  ,  &  
     &           cwirin (:999)   , cwold(:24)       , gcursv(:24)      ,  &  
     &           gmary(:24,:24)   , cdot(:24)        , dtcorr          ,  &  
     &            rctim  ,thrp   ,                                        &  
     &                xrplate(:pnplat)  , xlplate(:pnplat)  ,             &  
     &                zrplate(:pnplat)  , zlplate(:pnplat)  ,             &  
     &                strike(:pnplat,:2) , fplate(:pnplat,:2) ,           &  
     &                dsep(:pnplat,:2)   ,                                &  
     &                hplate(:pnplat,:pnseg) , dplate(:pnplat,:pnseg) ,   &  
     &                xsega(:pnplat,:pnseg+1), zsega(:pnplat,:pnseg+1),   &  
     &                  fincap

      CASE ("begc6a")
      read(lu) begc6a,                                                    &  
     &              rsceg  (:pngroup),                                    &  
!cj  &         ameg (:pngroup,:pngroup),   ameginv(:pngroup,:pngroup),    &  
     &         ((ameg(k,j),k=1,pngroup),j=1,pngroup),                     &  
     &         ((ameginv(k,j),k=1,pngroup),j=1,pngroup),                  &  
     &         ceg  (:pngroup),  cegp1  (:pngroup),  cegm1  (:pngroup),   &  
     &         ceg0 (:pngroup),  ceg0m1 (:pngroup),  aintceg(:pngroup),   &  
     &         flxeg(:pngroup),  flxegm1(:pngroup),  flxegm2(:pngroup),   &  
     &         veg  (:pngroup),  vegm1  (:pngroup),  aturneg(:pngroup),   &  
     &         voheg(:pngroup),  vohegm1(:pngroup),  ceg0m2(:pngroup),    &  
     &         vegplot(:pncoil), gvolt(:ptpts,:pngroup),                  &  
     &         teineg (:pngroup),tauseg (:pngroup),                       &  
     &         gainpeg(:pngroup),gaindeg(:pngroup),gainieg(:pngroup),     &  
     &         vmineg(:pngroup),vmaxeg(:pngroup),dtsave,                  &  
     &         dcur   (:pngroup),dcurm1 (:pngroup), dcurm2 (:pngroup),    &  
     &                 finc6a

      CASE ("begcot")
      read(lu) begcot,                                                    &  
     &         tlast ,                                                    &  
     &         tempc(:pncoil)  ,  dxcoil(:pncoil)  ,  dzcoil(:pncoil)  ,  &  
     &         fcu(:pncoil)    ,  fss(:pncoil)     ,  templ(:pngroup)  ,  &  
     &         ccics(:pncoil)  ,                                          &  
     &                    fincot

      CASE ("begcsa")
      read(lu) begcsa,                                                    &  
     &        ane(:ppsi), zeffa(:ppsi), avez(:ppsi), simpe(:ppsi),        &  
     &        te(:ppsi),ti(:ppsi),sradion(:ppsi),zeffa2(:ppsi),           &  
     &        sinzrec(:pnthe),ajavlh(:ppsi),ajpary(:ppsi),tpsir,tpsil,ajavfw(:ppsi),    &
     &        ajavec(:ppsi),ajavlh2(:ppsi),djdelh2(:ppsi),djdetsc2(:ppsi),              &  
     &        voltlp(:ppsi), anne(:ppsi), animp(:ppsi), sradiono(:ppsi),  &  
     & aneimp(:ppsi),avezimpa(:ppsi),amassimpa(:ppsi),amasshyda(:ppsi),   &  
     & aimassa(:ppsi), sumnia(:ppsi), sumnqa(:ppsi),                      &  
     &           tpsiri,tpsili,pden(:ppsi),ajavbs(:ppsi),ajavcd(:ppsi),   &  
     &  anhy(:ppsi), anox(:ppsi), anca(:ppsi), anfe(:ppsi), anbe(:ppsi),  &  
     &      powtsc(:ppsi),currtsc(:ppsi),djdetsc(:ppsi),djdelh(:ppsi),    &  
     &      powtsci(:pimp+2,ppsi), sravejet(:ppsi), sjane0, sumjet,       &  
     &                fincsa

      CASE ("begcim")
      read(lu) begcim,                                                    &  
     &                nby,nzspec(:3),fracim(:3),                          &  
     &                fincim

      CASE ("begcsp")
      read(lu) begcsp,                                                    &  
!cj  &         nq (:pchrgmx,:pimp,:ppsi),                                 &  
     &         (((nq(i,k,j),i=1,pchrgmx),k=1,pimp),j=1,ppsi),             &  
!cj  &         nqo(:pchrgmx,:pimp,:ppsi),                                 &  
     &         (((nqo(i,k,j),i=1,pchrgmx),k=1,pimp),j=1,ppsi),            &  
!cj  &         dperi(:pchrgmx,:ppsi)     ,                                &  
     &         ((dperi(k,j),k=1,pchrgmx),j=1,ppsi),                       &  
!cj  &         dpari(:pchrgmx,:ppsi),                                     &  
     &         ((dpari(k,j),k=1,pchrgmx),j=1,ppsi),                       &  
!cj  &         ainz(:pchrgmx,:ppsi)      ,                                &  
     &         ((ainz(k,j),k=1,pchrgmx),j=1,ppsi),                        &  
!cj  &         rec(:pchrgmx,:ppsi),                                       &  
     &         ((rec(k,j),k=1,pchrgmx),j=1,ppsi),                         &  
!cj  &         rad(:pimp,:pchrgmx,:ppsi),                                 &  
     &         (((rad(i,k,j),i=1,pimp),k=1,pchrgmx),j=1,ppsi),            &  
     &         radpel, vxpel, vzpel, xpel,                                &  
     &         zpel, dndtpel , tepel , totradp ,imppel, impbnd,           &  
     &                fincsp

      CASE ("begcrt")
      read(lu) begcrt,                                                    &  
     & alinzOx(:pcgOx,:pte),                                              &  
     & alradOx(:pcgOx,:pte,:pne),alrecOx(:pcgOx,:pte,:pne),               &  
     & alinzC (:pcgC ,:pte),                                              &  
     & alradC (:pcgC ,:pte,:pne),alrecC (:pcgC ,:pte,:pne),               &  
     & alinzFe(:pcgFe,:pte),                                              &  
     & alradFe(:pcgFe,:pte,:pne),alrecFe(:pcgFe,:pte,:pne),               &  
     & alinzBe(:pcgBe,:pte),                                              &  
     & alradBe(:pcgBe,:pte,:pne),alrecBe(:pcgBe,:pte,:pne),               &  
     & alinzNe(:pcgNe,:pte),                                              &  
     & alradNe(:pcgNe,:pte,:pne),alrecNe(:pcgNe,:pte,:pne),               &  
     & alinzKr(:pcgKr,:pte),                                              &  
     & alradKr(:pcgKr,:pte,:pne),alrecKr(:pcgKr,:pte,:pne),               &  
     & alinzAr(:pcgAr,:pte),                                              &  
     & alradAr(:pcgAr,:pte,:pne),alrecAr(:pcgAr,:pte,:pne),               &  
     & alinzW (:pcgW ,:pte),                                              &  
     & alradW (:pcgW ,:pte,:pne),alrecW (:pcgW ,:pte,:pne),               &  
     &    altei(:pte,:pimp), alnei(:pne,:pimp), nne(:pimp),               &  
     &    nte(:pimp),       nchrgsr(:pimp),                               &  
     &                  fincrt

      CASE ("plascur")
      read(lu) plascur, dipdt,                                            &  
     &                vvl, vvr, cvvpol, vvdum(:16), vvfin

      CASE ("begtcv")
      read(lu) begtcv,                                                    &  
     & apla,aplp,ax(:pnvvel,:pnvvel+1),aph(:pmmay),a7(:pnfob,:pnshap),    &  
     & a8(:pnfob,:pnvvel),a8l(:pnob,:pnvvel),berr(:pnprob),bee(:pnprob),  &  
     & bxob(:pnbob),bzob(:pnbob),coesv(:pnwire),coth(:pnprob),            &  
     & cerr(:picur),curr(:picur),caltime,cx(:pnvvel),cxold(:pnvvel),      &  
     &  cxave(:pnvvel),cxold2(:pnvvel),                                   &  
     & dzel,dex(:pnprob,:picur),dep(:pnprob,:pmmay),dsh(:pnx,:pnz),       &  
     & ddd(:pnshap),dddd(:pnshap),dissi,elss(:pnshap,:pnshap),            &  
     & elsv(:pnshap,:pnvvel),elvvi(:pnvvel,:pnvvel),eee(:pnob),           &  
     & epsza,epszo,epsqa,epsqo,elvv(:pnvvel,:pnvvel),                     &  
     & ferr(:pnluup),flxa(:pnluup),flxo(:pnluup),flxd(:pnluup),           &  
     & ffob(:pnob),gx(:pnluup,:picur),gp(:pnluup,:pmmay),hshap(:pnshap),  &  
     & hvvel(:pnvvel),icur,isvf(:picur),isvl(:picur),                     &  
     & iprtcv,ielch,mmax,nelz,nprob,nprtcv,nelch,nshac,mshap,             &  
     & nluup,nob,nfob,nbxob,nbzob,nshap,nvvel,nvsob,                      &  
     & rhoa(:pnob),rhoo(:pnob),                                           &  
     & rscurr(:picur),rshap(:pnshap),rvvel(:pnvvel),rhn(:pnvvel),         &  
     & rho(:pnob),sith(:pnprob),tmstp,timtcv2,thb(:pnprob),timesnw,       &  
     & timesld,time1,time2,volttcv(:pnwire),umeg(:pnob,:pnshap),uc,       &  
     & vssh2,vssh,vssha,vssh0,                                            &  
     & vvv(:pnprob),wshap(:pnshap),wvvel(:pnvvel),www(:pnluup),           &  
     & xf(:pnluup),xip,xop,xb(:pnprob),xfob(:pnfob),                      &  
     & xbxob(:pnbob),xbzob(:pnbob),                                       &  
     & xvsob,xniv(:pmmay,4),zniv(:pmmay),                                 &  
     & xhnn(:pmmay),zhnn(:pmmay),xhne(:pmmay),xhoo(:pmmay),zhoo(:pmmay),  &  
     & xhee(:pmmay),xhss(:pmmay),zhss(:pmmay),xhse(:pmmay),               &  
     & zf(:pnluup),zshap(:pnshap),zvvel(:pnvvel),zup,zlp,                 &  
     & zb(:pnprob),zfob(:pnfob),                                          &  
     & zbxob(:pnbob),zbzob(:pnbob),zvsob,                                 &  
     & fintcv

      CASE ("begufc")
      read(lu) begufc,                                                    &  
!    &          ufdata(:2*pnsave,:pnx,:30),                                 &
     &   (((ufdata(k,j,i),k=1,2*pnsave),j=1,pnx),i=1,30),                   &
     &                    finufc

      CASE ("begrun")
      read(lu) begrun                                                     &  
     &,         ajpre(:pnx,:pnz), ajprecc(:pnx,:pnz), ajpresf(:ppsi)      &  
     &,         anre(:ppsi), sresf(:ppsi), etafac(:ppsi)                  &  
     &,         adnre(:ppsi), recur, sumre, sreav                         &  
     &,                  finrun

      CASE ("begcwa")
      read(lu) begcwa,                                                    &  
!cj  &     xw (:pnthe,:ppsi),                                            &  
     &     ((xw(k,j),k=1,pnthe),j=1,ppsi),                               &  
!cj  &     zw (:pnthe,:ppsi)  ,                                          &  
     &     ((zw(k,j),k=1,pnthe),j=1,ppsi),                               &  
!cj  &     bmagfc(:pnthe,:ppsi),                                         &  
     &     ((bmagfc(k,j),k=1,pnthe),j=1,ppsi),                           &  
!cj  &     sqg(:pnthe ,:ppsi),                                           &  
     &     ((sqg(k,j),k=1,pnthe),j=1,ppsi),                              &  
!cj  &     g33(:pnthew,:ppsi),                                           &  
     &     ((g33(k,j),k=1,pnthew),j=1,ppsi),                             &  
!cj  &     g11(:pnthew,:ppsi),                                           &  
     &     ((g11(k,j),k=1,pnthew),j=1,ppsi),                             &  
!cj  &     g22(:pnthew,:ppsi),                                           &  
     &     ((g22(k,j),k=1,pnthew),j=1,ppsi),                             &  
!cj  &     g21(:pnthew,:ppsi),                                           &  
     &     ((g21(k,j),k=1,pnthew),j=1,ppsi),                             &  
!cj  &     g12(:pnthew,:ppsi),                                           &  
     &     ((g12(k,j),k=1,pnthew),j=1,ppsi),                             &  
!cj  &  xcentfc(:pnthe,:ppsi),                                           &  
     &     ((xcentfc(k,j),k=1,pnthe),j=1,ppsi),                          &  
!cj  &  apzone(:pnthe,:ppsi),                                            &  
     &     ((apzone(k,j),k=1,pnthe),j=1,ppsi),                           &  
!cj  &  aplost(:pnthe,:ppsi),                                            &  
     &     ((aplost(k,j),k=1,pnthe),j=1,ppsi),                           &  
!cj  &     zcentfc(:pnthe,:ppsi),                                        &  
     &     ((zcentfc(k,j),k=1,pnthe),j=1,ppsi),                          &  
!cj  &     ripcon(:pnthe,:ppsi),                                         &  
     &     ((ripcon(k,j),k=1,pnthe),j=1,ppsi),                           &  
!cj  &     zws(:pnthe,:ppsi),                                            &  
     &     ((zws(k,j),k=1,pnthe),j=1,ppsi),                              &  
!cj  &     aplostr(:pnthe,:ppsi),                                        &  
     &     ((aplostr(k,j),k=1,pnthe),j=1,ppsi),                          &
     &                fincwa


      CASE DEFAULT

      print *, "Error in buf_in: match not found", trim(idrss)

      END SELECT

      return
      END SUBROUTINE buf_in
