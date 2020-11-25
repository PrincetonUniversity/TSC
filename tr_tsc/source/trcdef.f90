!#include "f77_dcomplx.h"
      subroutine trcdef(i_newton)
!
!.....define as, cs, dsi, dse  coefficients for
!.....particular transport model - - surface centered
!
      USE trtsc, ONLY : use_user_rot,use_user_chie,use_user_pres
      USE CLINAM
      USE BALCLI
      USE RUNAWAY
      USE SAPROP
      USE SCR3
      USE NEWPLOT
      USE modmmm7_1

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,jpeak,npsihm,npsi9,lsaw,k,i,ii,jshoot,jmm,jmaxm
      INTEGER i_grad,i_delay,leigen,idengrad,nroot,iglf,irotstab
      INTEGER ibt_flag,npmax,nprint,matdim,nprout,lprint,lsuper
      INTEGER lreset,nerr,npoints,i_newton

!..Bateman, 15 March 2009, use iprint > 0 to control diagnostic output

      INTEGER :: iprint = 0

!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 slight,umu,aionm,chired,chicoef,chimax,elchar,sum2
      REAL*8 sum3,sum3i,sum4,sum5,sum6,sum7,sum8,sum9,cswitchre
      REAL*8 dvol,term,term2,term3,term4,enercn,ptot,wdot,temax
      REAL*8 formax,fp2,btmax,fcyc,fcyc2,fcyci,fcyci2,fpi2,y2x
      REAL*8 y2rat,accn,tfluxb,temid,anemid,timid,animid,timide
      REAL*8 avezm,allam,etperac,cr,si,alame,delt,rfac,etpera,w2lim
      REAL*8 zdiff,wratf,pnorm,psmid,esmid,aimh,rmid,emid,pmid
      REAL*8 temidd,timidd,face,faci,gps,densfac,gpsd,coef1,coef2
      REAL*8 coef3,coef4,coef5,chiesecj,chiisecj,diffsecj,formfs
      REAL*8 formis,an0,anel,vprime,tflux,chiemks,chiimks,floatj
      REAL*8 dperp,balfac,befo,torcur,chicmg1,chicmg2,chicmg3
      REAL*8 chicmg4,chicmg,pchinum,pchiden,pchirat,chi0,d0neo
      REAL*8 ailarr,allii,tauii,gval,gpval,gppval,bzero,collfii,aki
      REAL*8 chiemin,alarsq,animd1,tcoli,alfnc,tmsti,aki1,aki2,grs
      REAL*8 facs,zpte_in,zpti_in,zpne_in,zpni_in,bt_exp,arho_exp
      REAL*8 rmajor_exp,amassgas_exp,zimp_exp,amassimp_exp,alpha_e
      REAL*8 x_alpha,diffnem,chietem,chiitim,etaphim,etaparm
      REAL*8 etaperm,exchm,relfactor,relfactor2,rminormm,rmajormm
      REAL*8 elongmm,densemm,denshmm,densimpmm,densfemm,zeffmm
      REAL*8 tekevmm,tikevmm,qmm,gvalmm,gpvalmm,gppvalmm,btormm
      REAL*8 avezimpmm,amassimpmm,amasshydmm,aimassmm,wexbsmm
      REAL*8 befogrd,grdnemm,grdnimm,grdnhmm,grdnzmm,grdtemm
      REAL*8 grdtimm,grdqmm,thiig,thdig,theig,thzig,thirb,thdrb
      REAL*8 zcsnd0
      REAL*8 therb,thzrb,thikb,thdkb,thekb,thzkb,sqrgps,bsq,al0
      REAL*8 alstar,alinc,ay,coef3nc,chiinc,chienc,coef6
      REAL*8 sum, AREAL, ptoti, felec, fion
      REAL*8, DIMENSION(1) ::  zrminormm, zrmajormm, zelongmm &
        & , zdensemm, zdenshmm, zdensimpmm, zdensfemm &
        & , zzeffmm,  ztekevmm, ztikevmm, zqmm, zbtormm &
        & , zavezimpmm, zamassimpmm, zamasshydmm, zaimassmm &
        & , zwexbsmm, zbefogrd, zgrdnemm, zgrdnimm, zgrdnhmm &
        & , zgrdnzmm, zgrdtemm, zgrdtimm, zgrdqmm &
        & , zgvrin, zvtorin, zgvpin, zvpolin, zgvparin, zvparin &
        & , zeta,     zdkdeps,  zthtig,   zthttig &
        & , zthiig,   zthdig,   ztheig,   zthzig

!============ 
!     common /newplot/ chienca(ppsi),chiinca(ppsi),
!    1                 chiicopi(ppsi),chiecopi(ppsi),diffary(ppsi)
!
!     dimension formf(ppsi),aemid(ppsi),armid(ppsi),apmid(ppsi),
!    1           aface(ppsi),afaci(ppsi),achi0(ppsi),formi(ppsi),
!    2           chiineo(ppsi),ache0(ppsi)
!
!.....arrays needed for glf23
!     REAL*8 te_m(0:ppsi),ti_m(0:ppsi),ne_m(0:ppsi),ni_m(0:ppsi),
!    1       ns_m(0:ppsi),angrotp_exp(0:ppsi),egamma_exp(0:ppsi),
!    2       gamma_p_exp(0:ppsi),vphi_m(0:ppsi),vpar_m(0:ppsi),
!    3       vper_m(0:ppsi),zeff_exp(0:ppsi),rho(0:ppsi),
!    4       gradrho_exp(0:ppsi),gradrhosq_exp(0:ppsi),
!    5       q_exp(0:ppsi),shat_exp(0:ppsi),alpha_exp(0:ppsi),
!    6       elong_exp(0:ppsi),
!    7       diff_m(0:ppsi),chie_m(0:ppsi),chii_m(0:ppsi),
!    8       etaphi_m(0:ppsi),etapar_m(0:ppsi),etaper_m(0:ppsi),
!    9       exch_m(0:ppsi), egamma_m(0:ppsi), egamma_d(0:ppsi),
!    *       gamma_p_m(0:ppsi),anrate_m(0:ppsi),anrate2_m(0:ppsi),
!    1       anfreq_m(0:ppsi),anfreq2_m(0:ppsi),rmin_exp(0:ppsi),
!    2       rmaj_exp(0:ppsi)
      integer itport_pt(5)
!
!.....arrays needed for mmm95
      REAL*8 cswitch(23),fig(4),frb(4),fkb(4),difthi(5,5,1),velthi(5,1),  &  
     &       vflux(5,1),gamma(5,1),omega(5,1)
      integer lswitch(8)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: zchii
      REAL*8, ALLOCATABLE, DIMENSION(:) :: formf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aemid
      REAL*8, ALLOCATABLE, DIMENSION(:) :: armid
      REAL*8, ALLOCATABLE, DIMENSION(:) :: apmid
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aface
      REAL*8, ALLOCATABLE, DIMENSION(:) :: afaci
      REAL*8, ALLOCATABLE, DIMENSION(:) :: achi0
      REAL*8, ALLOCATABLE, DIMENSION(:) :: formi
      REAL*8, ALLOCATABLE, DIMENSION(:) :: chiineo
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ache0
      REAL*8, ALLOCATABLE, DIMENSION(:) :: te_m,ti_m,ne_m,ni_m,          &  
     &       ns_m,angrotp_exp,egamma_exp,                                &  
     &       gamma_p_exp,vphi_m,vpar_m,                                  &  
     &       vper_m,zeff_exp,rho,                                        &  
     &       gradrho_exp,gradrhosq_exp,                                  &  
     &       q_exp,shat_exp,alpha_exp,                                   &  
     &       elong_exp,                                                  &  
     &       diff_m,chie_m,chii_m,                                       &  
     &       etaphi_m,etapar_m,etaper_m,                                 &  
     &       exch_m, egamma_m, egamma_d,                                 &  
     &       gamma_p_m,anrate_m,anrate2_m,                               &  
     &       anfreq_m,anfreq2_m,rmin_exp,                                &  
     &       rmaj_exp
      REAL*8 :: temp1,temp2,temp21,temp22,temp31,temp32,temp41,temp42
      REAL*8, ALLOCATABLE, DIMENSION(:) :: te_m_sav,ti_m_sav
      REAL*8, ALLOCATABLE, DIMENSION(:) :: te_sav,ti_sav, ane_sav
      REAL*8 :: chimin,chimin_tmp,delta_t,t_j,t_j1
      INTEGER :: jstart
!============      
!     Jan 2012 - fmp - variables to read thermal conductivity from file
      INTEGER ntime_chie,nrad_chie
      LOGICAL :: ex_chie=.false.
      LOGICAL :: first_read_chie=.true.
      REAL*8 pchi1,pchi2,chie_read
      REAL*8, ALLOCATABLE, DIMENSION(:,:):: prof_chie
      REAL*8, ALLOCATABLE, DIMENSION(:) :: time_chie,rad_chie
!============      
      logical :: first_read_rot=.true.
      logical :: ex_rot=.false. 
      real*8 :: tint, rint, wrot1, wrot2, frot,tfluxd,x1,x2,            &
     &          pres1,pres2,presint
      real*8, dimension(:),allocatable :: time_rot, rad_rot
      real*8, dimension(:,:),allocatable :: wrot
      integer :: j1, j2, kk, l, lsav
      integer :: ntime_rot, nrad_rot
!============  
!     arrays needed for the CDBM model
      integer jj,icdbmodel,jmin,jmin_cdbm
      REAL*8 xchi,chiemksCT,chiimksCT,qmin_CT
      INTEGER isw,npsiCDBM,npsiCT,npsiCDBMmin
      REAL*8, DIMENSION(npsit) :: chicdbm,chicdbm0,chiicdbms,chiecdbms
!============  
      IF(.not.ALLOCATED(formf)) ALLOCATE( zchii(ppsi,5), STAT=istat)
      IF(.not.ALLOCATED(formf)) ALLOCATE( formf(ppsi), STAT=istat)
      IF(.not.ALLOCATED(aemid)) ALLOCATE( aemid(ppsi), STAT=istat)
      IF(.not.ALLOCATED(armid)) ALLOCATE( armid(ppsi), STAT=istat)
      IF(.not.ALLOCATED(apmid)) ALLOCATE( apmid(ppsi), STAT=istat)
      IF(.not.ALLOCATED(aface)) ALLOCATE( aface(ppsi), STAT=istat)
      IF(.not.ALLOCATED(afaci)) ALLOCATE( afaci(ppsi), STAT=istat)
      IF(.not.ALLOCATED(achi0)) ALLOCATE( achi0(ppsi), STAT=istat)
      IF(.not.ALLOCATED(formi)) ALLOCATE( formi(ppsi), STAT=istat)
      IF(.not.ALLOCATED(chiineo)) ALLOCATE( chiineo(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ache0)) ALLOCATE( ache0(ppsi), STAT=istat)
      IF(.not.ALLOCATED(te_m)) ALLOCATE( te_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(ti_m)) ALLOCATE( ti_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(ne_m)) ALLOCATE( ne_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(ni_m)) ALLOCATE( ni_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(ns_m)) ALLOCATE( ns_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(angrotp_exp)) ALLOCATE( angrotp_exp(0:ppsi),     &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(egamma_exp)) ALLOCATE( egamma_exp(0:ppsi),       &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(gamma_p_exp)) ALLOCATE( gamma_p_exp(0:ppsi),     &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(vphi_m)) ALLOCATE( vphi_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(vpar_m)) ALLOCATE( vpar_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(vper_m)) ALLOCATE( vper_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(zeff_exp)) ALLOCATE( zeff_exp(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(rho)) ALLOCATE( rho(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(gradrho_exp)) ALLOCATE( gradrho_exp(0:ppsi),     &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(gradrhosq_exp)) ALLOCATE( gradrhosq_exp(0:ppsi)  &  
     &                                   , STAT=istat)
      IF(.not.ALLOCATED(q_exp)) ALLOCATE( q_exp(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(shat_exp)) ALLOCATE( shat_exp(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(alpha_exp)) ALLOCATE( alpha_exp(0:ppsi),         &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(elong_exp)) ALLOCATE( elong_exp(0:ppsi),         &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(diff_m)) ALLOCATE( diff_m(0:ppsi),               &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(chie_m)) ALLOCATE( chie_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(chii_m)) ALLOCATE( chii_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(etaphi_m)) ALLOCATE( etaphi_m(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(etapar_m)) ALLOCATE( etapar_m(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(etaper_m)) ALLOCATE( etaper_m(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(exch_m)) ALLOCATE( exch_m(0:ppsi), STAT=istat)
      IF(.not.ALLOCATED(egamma_m)) ALLOCATE( egamma_m(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(egamma_d)) ALLOCATE( egamma_d(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(gamma_p_m)) ALLOCATE( gamma_p_m(0:ppsi),         &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(anrate_m)) ALLOCATE( anrate_m(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(anrate2_m)) ALLOCATE( anrate2_m(0:ppsi),         &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(anfreq_m)) ALLOCATE( anfreq_m(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(anfreq2_m)) ALLOCATE( anfreq2_m(0:ppsi),         &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(rmin_exp)) ALLOCATE( rmin_exp(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(rmaj_exp)) ALLOCATE( rmaj_exp(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(te_m_sav)) ALLOCATE( te_m_sav(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(ti_m_sav)) ALLOCATE( ti_m_sav(0:ppsi),           &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(ti_sav)) ALLOCATE( ti_sav(1:ppsi),               &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(te_sav)) ALLOCATE( te_sav(1:ppsi),               &  
     &                                   STAT=istat)
      IF(.not.ALLOCATED(ane_sav)) ALLOCATE( ane_sav(1:ppsi),             &  
     &                                   STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : trcdef  ' 
!============      
!
!
!.....constants used in the neoclassical ion transport model
      slight = 2.9979E10_R8
      umu = 1.66043E-24_R8
      aionm = (2.014_R8+3.016_R8)/2*umu
      chired = 1.0_R8
      chicoef = 2.25_R8
      chimax = 50.E6_R8/ max(apl,1._R8)
      chimin=1.0e-8
!     chimax=4.*chimax
!
      te_m = 0
      ti_m = 0
      te_m_sav = 0
      ti_m_sav = 0 
!
      elchar = 4.8032E-10_R8
      if(kcycle.le.1) then
      do j=1,ppsi
      chiitima(j) = 1.E-8_R8
      chietema(j) = 1.E-8_R8
      chiicopi(j) = 0.0_R8
      chiecopi(j) = 0._R8
      chienca(j) = 0._R8
      chiinca(j) = 0._R8
      enddo
      go to 309
      endif
      do j=1,ppsi
      chiitimao(j) = chiitima(j)
      chietemao(j) = chietema(j)
      diffnemao(j) = diffnema(j)
      enddo
      chiitimao(1) = chiitimao(2)
      chietemao(1) = chietemao(2)
      diffnemao(1) = diffnemao(2)
!
      zchii = 0.0_R8
!
!.....calculate global integrals and form factors
      sum = 0._R8
      sum2 = 0._R8
      sum3 = 0._R8
      sum3i= 0._R8
      sum4 = 0._R8
      sum5 = 0._R8
      sum6 = 0._R8
      sum7 = 0._R8
      sum8 = 0._R8
      sum9 = 0._R8
      cswitchre = 0._R8
      formf(1) = 0._R8
      formi(1) = 0._R8
      do 300 j=2,npsit
!
!
      dvol = vary(j)-vary(j-1)
!
      term = (.5_R8*(as(j-1)+as(j)) +vlooph(j) )/tpi                     &  
     &     *(gxmja2(j)-gxmja2(j-1))
      term2 = dvol*(savee(j) + savea(j)*ialpha                           &  
!    1 - sradion(j)*usdp/usdt - savebre(j) - savecyc(j) - saveimp(j)
!.....above line commented out 10/31/96....scj
     & + savefw(j) + savebm(j)  + savelh(j))
      term3 = dvol*(savei(j) + savia(j)*ialpha + savibm(j) + savilh(j)   &
     & +savifw(j) )
      term4 =                                                            &  
     &  dvol*(equila(j)*((1._R8+avez(j))*ade(j)-avez(j)*adp(j)))/vpg(j)
      if(term .le. 0) term=0._R8
!     if(term2.le. 0) term2=0.
!
      sum = sum + term
      sum2 = sum2 + 1.5_R8*dvol*adp(j)/vpg(j)
      sum3 = sum3 + term2  + term + term3
      sum3i = sum3i + term3 + term4
      sum4 = sum4 + te(j)*dvol
      sum5 = sum5 + ti(j)*dvol
      sum6 = sum6 + dvol
      sum7 = sum7 + dvol*(savea(j) + savia(j))
      sum8 = sum8 + dvol*(savei(j)+savee(j)+savifw(j)                    &  
     & + savefw(j) + savebm(j) + savibm(j) + savelh(j) + savilh(j))
      sum9 = sum9 + dvol*(                                               &  
     & - sradion(j)*usdp/usdt - savebre(j) - savecyc(j) - saveimp(j) )
!
  304 continue
      formf(j) = formf(j-1) + term + term2 + term3
      if(formf(j).lt.0.0_R8) formf(j) = 0._R8
      formi(j) = formi(j-1)  + term3 + term4
      if(formi(j).lt.0.0_R8) formi(j) = 0._R8
  300 continue
!cj      if(iecrh > 0) then
!cj         do l=1,ntpts-1
!cj            lsav = l
!cj            if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 410
!cj         enddo
!cj  410 continue
!cj         sum8 = sum8 + pecrh(lsav)
!cj      endif
      if(sum4.gt.0.0_R8.and.sum5.gt.0.0_R8) go to 305

!cj!cj...sepcial debug output
!cj      write(43,4343) sum4,sum5,sum6
!cj      do j=2,npsit
!cj        write(43,4344) j,vary(j),te(j),ti(j)
!cj      enddo
!cj 4343 format(" sum4, sum5, sum6", 1p3e12.4,/,"  j         vary te          ti")
!cj 4344 format(i3,1p3e12.4)
      ineg=43
      return
  305 pohmic = sum*udsp/udst
      if(pohmic.lt.0.0_R8) pohmic = 0.0_R8
      enercn = sum2*udsp
      teform = te(2)*sum6/sum4
      tiform = ti(2)*sum6/sum5
      palpha = sum7*ialpha*udsp/udst
      palphap = sum7*udsp/udst
      prad = sum9*udsp/udst
      paux = sum8*udsp/udst
!
!     write(*,*) "--- teform ---", kcycle, teform, te(2), tiform, ti(2) 
!
      if(sum3.le.1.E-12_R8) sum3 = 1.E-12_R8
      if(sum3i.le.1.E-12_R8) sum3i = 1.E-12_R8
      ptot = sum3*udsp/udst
      ptoti = sum3i*udsp/udst
!
!.....NOTE:  use Coppi/Tang during ohmic phase for itrmod .ge.8
      itrmode = itrmod
      if((paux .eq. 0.0_R8 .and. palpha .lt. 10.E6_R8)                   &  
     &                .and. itrmod.ge.8)       itrmode = 2
      if(itrmod.eq.12) itrmode = 2
      if(itrmod.eq.13) itrmode = 2
      if(itrmod.eq.14) itrmode = 2
      if(itrmod.eq.16) itrmode = 2
      if(itrmod.eq.17) itrmode = 2
!
!...==> special diagnostic
      if(iplt .ge. nskipl) then
      if(ialpha.gt.0) write(nterm,6661) paux,palpha,itrmod,itrmode
 6661 format("paux,palpha,itrmod,itrmode=",1p2e12.4,0p2i6)
      endif
!
!.....special fix for disruption modeling
      if(times .gt. acoef(95)) qsaw = acoef(96)
!
      wdot = wdotmw*1.E6_R8*usdp/usdt
      tauems = (sum2/(sum3-wdot))*udst*1000._R8
  309 continue
!
!.....increase chi inside max tempsurface
      temax = 0._R8
      do 1309 j=2,npsit-1
      if(te(j).le.temax) go to 1309
      temax = te(j)
      jpeak = j
 1309 continue
      if(jpeak.eq.2) go to 1409
      formax = formf(jpeak+1)
      do 1318 j=2,jpeak+1
 1318 formax = max(formax,formf(j))
      do 1319 j=2,jpeak
 1319 formf(j) = formax
 1409 continue
!
!
!
      do 3999 j=npsit+1,npsi
 3999 etpara(j) = etav
!.....If ilhcd is not zero, a power dependent contribution to the
!.....resistivity, resulting from the hot plasma population, is
!.....calculated.  sighot is the ratio of sigma_hot obtained from
!.....Eq. 11 of Fisch, Phys. Fluids 28, 245 (1985), to sigma_spitzer
!.....obtained from the Hirshman et al. paper refered to below.
!.....sigma_spitzer is the Hirshman result for conductivity when trapping
!.....is absent.
      IF (ilhcd .GT. 0) THEN
        fp2 = 806.2_R8* r0 / freqlh**2
        btmax = ABS(gs(imag,jmag))/xmag
        fcyc = 28.0_R8 * btmax / freqlh
        fcyc2 = fcyc * fcyc
        fcyci = fcyc * zion / (aion * 1836._R8)
        fcyci2 = fcyci * fcyci
        fpi2 = fp2 * zion**2 / (aion * 1836._R8)
        y2x = 1.0_R8/ (fcyci * fcyc)
        y2rat = y2x / (1.0_R8- y2x)
        IF ( y2x .GE. 1.0_R8.OR. fpi2 .LT. y2rat ) THEN
        accn = SQRT(fp2)/fcyc + SQRT(1._R8+ fp2/fcyc2 - fpi2)
        ELSE
        accn = SQRT( 1._R8/(1._R8- y2x) )
        ENDIF
      ENDIF
      tfluxb = (npsit-1)*dpsi
!
!----------------------------------------------------------------------------
!  read thermal diffusivity profiles from external file, if provided
!
      if(use_user_chie .and. first_read_chie) then
         inquire(file="user_chie_data", exist=ex_chie)
         if (ex_chie) then
            open(55,file="user_chie_data",form="formatted",status="old")
            read(55,*) ntime_chie, nrad_chie
            allocate(time_chie(ntime_chie), rad_chie(nrad_chie))
            allocate(prof_chie(nrad_chie,ntime_chie))
            read(55,*) time_chie(1:ntime_chie)
            read(55,*) rad_chie(1:nrad_chie)
            read(55,*) prof_chie(1:nrad_chie,1:ntime_chie)
            close (55)
         endif
         first_read_chie=.false.
      endif
!
!----------------------------------------------------------------------------
!	calculations for the CDBM model
      if(itrmod .eq. 16) call cdbm_mod(chicdbm)

!----------------------------------------------------------------------------
!.....start loop to define the arrays:
!     etpara, equila, csj , dsij , dsej     j=0,4
!     (note:  sawtooth corrections are added outside the loop)
!
      do 4000 j=1,npsit
!
!
!.....define surface averaged resistivity and equilibration
      temid = .5_R8*(te(j+1)+te(j))
      anemid = .5_R8*(ane(j+1)+ane(j))
      timid = .5_R8*(ti(j+1)+ti(j))
      if(timid.le.0.1_R8*temid) timid = 0.1_R8*temid
      animid = .5_R8*(anhy(j+1)+anhy(j))*1.E-6_R8
      timide = timid*1.6022E-12_R8
      avezm = .5_R8*(avez(j+1)+avez(j))
      allam = 24._R8-log(sqrt(anemid*1.E-6_R8)/temid)
      etperac= 1.03E-4_R8*allam*temid**(-1.5_R8)*usdr
!
!.....add neoclassical corrections
!
!......neoclassical conductivity from hirshman,hawryluk, nuclear fusion
!......                               17 (page 611) 1977
!
      cr = (0.56_R8/zeff)*(3.0_R8-zeff)/(3.0_R8+zeff)
      si = 0.58_R8+ 0.20_R8*zeff
      alame = (3.40_R8/zeff)*(1.13_R8+zeff)/(2.67_R8+zeff)
      delt = sqrt(vary(j)/(tpi*pi*xmag))/xmag
!
      anue(j) = 0._R8
      if(ftrap(j) .gt. 0.0_R8)                                           &  
     &  anue(j) = (xmag*qprof2(j)*anemid*allam)                          &  
     &          / (ftrap(j)*delt*temid**2*10.2E16_R8)
      rfac = alame*(1.0_R8-ftrap(j)/(1.0_R8+si*anue(j)))                 &  
     &     *       (1.0_R8-cr*ftrap(j)/(1.0_R8+si*anue(j)) )
!
!.....added 7/29/10 to be consistent with halo for disruption calculations
      if(acoef(109).gt.0.) rfac = 2.0/zeff
!     resred(j) = rfac
      etpera = etperac/rfac * acoef(810)
      sighot(j+1)=0._R8
!
      IF (ilhcd .GT. 0 .and. ifk.eq.0) THEN
!.......calculate w2lim based on lh accessibility condition.
!
!.......Normalize w2lim by v-thermal
        w2lim = 7.15E2_R8/(accn*SQRT(temid) )
!.......Since raytracing results show that w2lim become no larger than 10
!.......when w2lim is greater than 10 we replace it by 10.  Larger values
!.......of w2lim are likely to correspond to a runaway region.
        IF (w2lim .GT. 10.0_R8) w2lim = 10.0_R8
!.....Check that w2lim is no less than w1lim + 1.0e-3
        zdiff = ABS(w1lim)+ 1.0E-3_R8
        IF (w2lim .LT. zdiff) w2lim = zdiff
        IF (w1lim .LT. 0.0_R8) w2lim = -w2lim
        wratf = (w2lim**4 - w1lim**4) / (4.0_R8*log(w2lim/w1lim))
        pnorm = 260.2_R8* savelh(j+1)*(udsp/udst)*SQRT(temid)/anemid**2
        sighot(j+1)=( 8.0_R8/ (7.52_R8*alame *(zeff+3.0_R8)*allam*       &
     &  9.11E-28_R8) )*                                                  &  
     &              pnorm * wratf
      ENDIF
!
      etpera = min( etav , etpera / (1.0_R8+ sighot(j+1) ))
!
      if(temid.lt.30.0_R8.and. kcycle.gt.0)  then
        etpara(j) = 0.10_R8*(etpera - etpara(j)) + etpera
      else
        etpara(j) = etpera
      endif
      if(whalos.gt.0 .and. etpara(j).gt.etah) etpara(j)=etah
!
      equila(j) = anemid*3.1E-11_R8*etperac*zeffa2(j)*udsr/usdt
!
      if(kcycle.le.0) go to 4000
!
!     write(89,1088) kcycle,j,etperac,temid,ftrap(j),anue(j),etpera,etpara(j),etpara(j)*udsr
!1088 format(2i5,1p8e11.3)
!
      as0(j) = 0._R8
      as2(j) = 0._R8
      as1(j) = 0._R8
      as3(j) = 0._R8
!
!
      cs0(j) = 0._R8
      cs4(j) = 0._R8
      cs2(j) = 0._R8
      cs1(j) = 0._R8
      cs3(j) = 0._R8
!
!
      dsi0(j) = 0._R8
      dsi4(j) = 0._R8
      dsi2(j) = 0._R8
      dsi1(j) = 0._R8
      dsi3(j) = 0._R8
!
      dse0(j) = 0._R8
      dse4(j) = 0._R8
      dse2(j) = 0._R8
      dse1(j) = 0._R8
      dse3(j) = 0._R8
!
      psmid = .5_R8*(adp(j)+adp(j+1))
      esmid = .5_R8*(ade(j)+ade(j+1))
!
      aimh = .5_R8*(adi(j)+adi(j+1))
      rmid = .5_R8*(adn(j+1)/vp(j+1)+adn(j)/vp(j))
      emid = .5_R8*(ade(j+1)/vpg(j+1)+ade(j)/vpg(j))
      pmid = .5_R8*(adp(j+1)/vpg(j+1)+adp(j)/vpg(j))
      temidd = emid/rmid
      timidd = (pmid-emid)/rmid
      if(timidd.le.0.0_R8) timidd = 0.1_R8*temidd
!
!
!....note 5/2 changed to 3/2 on 11/04/99...scj
!....    /rmid removed on 05/25/00...scj
      face = emid*acoef(897)
      faci = (pmid-emid)*acoef(897)
!
!.....note:  gps is the surface averaged gradient of the toroidal flux squared
      gps = gja2(j)/vp2(j)*(tpi*qprof2(j))
      d2s(j) = gps
!
!     itrmod = 3 for the coppi mazzucato gruber transport model
!     itrmod = 4 for the Englade transport model
!     itrmod = 5 for marion Turner
!     itrmod = 6 for GLF2D + neoclassical
!     itrmod = 7 for         neoclassical only
!     itrmod = 8 for GLF2D + Coppi-Tang
!     itrmod = 9 for GLF2D + neo-Alcator
!     itrmod =10 for MMM95 + Coppi-Tang
!     itrmod =11 for MMM95 + neo-Alcator
!     itrmod =12 for Kessel-modified subroutine kcoppi
!     itrmod =15 for MMM_v7.1
!     itrmod =16 for CDBM model
!     itrmod =17 for Coppi-Tang+CDBM 
!....define flux surface where top of pedistal is to be applied
      npsihm = int(pwidthc * AREAL(npsit))

      IF ( 15 == itrmod ) GO TO 15
      go to(1,2,3,4,5,6,7,2,1,2,1),itrmode
!
    1 continue
!.......................................................................
!.....simple emperical model  ... goldstone-kaye model
!......................................................................
      densfac = adn(1)/vp(1)*udsd/1.E19_R8
      chiauxs = (acoef(107)*(densfac)*300._R8/apl)**2*ptot
!
      gpsd = (gps/usdt)*1.E19_R8/(udsd*rmid)
!.....particle diffusion
      coef1 = acoef(39)*gpsd/rmid
!.....electron thermal conductivity
      coef2 =-sqrt(acoef(35)**2+chiauxs)*gpsd
!.....ion thermal conductivity
      coef3 =-sqrt(acoef(37)**2+chiauxs)*gpsd
!.....off diagonal terms in electron and ion equations.
      coef4 =-sqrt(acoef(36)**2+chiauxs)*gpsd
      coef5 =-sqrt(acoef(38)**2+chiauxs)*gpsd
!
      if(gps.ne.0.0_R8) then
      chiesecj = -coef2/(udst*gps)
      chiisecj = -coef3/(udst*gps)
      diffsecj =  coef1*rmid/(udst*gps)
      endif
!
!......definitions for use with itrmode=9 and for plotting
      chiecopi(j) = chiesecj
      chiicopi(j) = chiisecj
      diffary(j) = diffsecj
!
!.....add GLF23 for itrmode=9 and MMM for itrmode=11
      if(itrmode.eq.9 ) go to 6
      if(itrmod.eq.10) go to 10
!
!.....particle diffusion
!
      cs1(j) = coef1
!
      dse1(j) = -coef2*emid/rmid - face*cs1(j)
      dse3(j) = coef2
!
!.....off-diagonal terms in the electron temperature equation
      dse1(j) = dse1(j) - coef4*(pmid-emid)/rmid
      dse2(j) = dse2(j) + coef4
      dse3(j) = dse3(j) - coef4
!
      dsi1(j) = -coef3*(pmid-emid)/rmid - faci*cs1(j)
      dsi2(j) = coef3
      dsi3(j) =-coef3
!
!.....off-diagonal terms in the ion temperature equation
      dsi1(j) = dsi1(j) - coef5*emid/rmid
      dsi3(j) = dsi3(j) + coef5
!
      go to 4000
    2 continue
!.......................................................................
!.....coppi/tang profile consistency model
!.......................................................................
!
      npsi9 = int(0.95_R8* AREAL(npsit) )
      q95pct = qprof2(npsi9)
      formfs = formf(j)/sum3
      formis = formi(j)/sum3i
      an0 = adn(2)/vp(2)*udsd
      anel = .5_R8*(adn(j+1)/vp(j+1)+adn(j)/vp(j))*udsd
      vprime = vp2(j)
      tflux = (j-1)*dpsi

!------------------------------------------------------------------
!
     
      SELECT CASE(itrmod)
         CASE(:11)      
           call coppi ( chiemks, chiimks, chiauxs, chiohms,             &
     &       an0,anel,xmag,tfluxb,tflux,ptot,formfs,formis,vprime,gps,  &
     &       q95pct,arad,gzero,zeff,alphar,acoef(121),acoef(122),qadd,  &
     &       acoef(124),acoef(126),acoef(3015),tfluxs,lsaw,dpsi,fhmodei,&
     &       pwidthc,firitb,secitb,chiped,npsihm,global(132))
         CASE(12)
            call kcoppi ( chiemks, chiimks, chiauxs, chiohms,           &
     &        an0,anel,xmag,tfluxb,tflux,ptot,formfs,formis,vprime,gps, &
     &        q95pct,arad,gzero,zeff,alphar,acoef(121),acoef(122),qadd, &
     &        acoef(124),acoef(126),tfluxs,lsaw,dpsi,fhmodei,pwidthc,   &
     &        firitb,secitb,chiped,npsihm,global(132))
         CASE(13,14)
!     NOTE:  subroutine d3dcoppi called for itrmod=13
!            need to supply it
!
!     call d3dcoppi ( chiemks, chiimks, chiauxs, chiohms,                &
!    & an0,anel,xmag,tfluxb,tflux,ptot,formfs,formis,vprime,gps,         &
!    & q95pct,arad,gzero,zeff,alphar,acoef(121),acoef(122),qadd,   &
!    & acoef(124),acoef(126),tfluxs,lsaw,dpsi,fhmodei,pwidthc,   &
!    & firitb,secitb,chiped,npsihm,global(132),times,     &
!    & felec,fion,ptoti)

         CASE(16)
           jmin_cdbm = 1
           qmin_CT = qprof2(3)
           do kk=3,npsit
              if (qprof2(kk) .le. qmin_CT) then
                 qmin_CT = qprof2(kk)
                 jmin_cdbm = kk
              endif
           enddo
           xchi = secitb
           if (secitb .gt. 0.) xchi = float(jmin_cdbm)/float(npsit)
           call coppi ( chiemks, chiimks, chiauxs, chiohms,             &
     &       an0,anel,xmag,tfluxb,tflux,ptot,formfs,formis,vprime,gps,  &
     &       q95pct,arad,gzero,zeff,alphar,acoef(121),acoef(122),qadd,  &
     &       acoef(124),acoef(126),acoef(3015),tfluxs,lsaw,dpsi,fhmodei,&
     &       pwidthc,firitb,secitb,chiped,npsihm,global(132))
!...neoclassical tansport
            jj = j
            if(j .eq. 1) jj = 2
            call geval(xsv2(jj),2,gval,gpval,gppval,imag,jmag)
            bsq = bsqar(jj)*1.e-8
            al0 = rmid*etpera*gxmja2(jj)/xmja2(jj)/bsq
            alstar = ftrap(jj)*gval**2*rmid*etpera/bsq
            alinc = 61.*(temidd/timidd)**1.5                            &
     &            *(al0*(1.+2.*qprof2(jj)**2)+0.46*alstar               &
     &            /(1.-0.54*ftrap(jj)))
            coef3nc = (tpi*qprof2(jj))**2*timidd*alinc
            chiinc = coef3nc/(udst*gps)
            if (j .eq. 1) chiinc = 0.5
            chienc = chiinc/(sqrt(2.*1835))
            chiimksCT = chiimks
            chiemksCT = chiemks
            npsiCDBM = int(acoef(3204)*AREAL(npsit))
            npsiCDBMmin = int(0.5*acoef(3204)*AREAL(npsit))
            npsiCT = int(acoef(3203)*AREAL(npsit))
            if (npsiCDBM .eq. npsiCT) npsiCDBM = npsiCDBM-3
!           discard CDBM calculations at the pedestal and use Coppi-Tang instead
            xchi = 0._R8
            if (j .le. npsiCDBMmin) then
               chiicdbms(j) = chicdbm(j)+chiinc
               chiecdbms(j) = chicdbm(j)+chiinc
              chiemks = chicdbm(j) + chiinc
              chiimks = chicdbm(j) + chiinc
               chiimks = chiicdbms(j)
               chiemks = chiimks
!           elseif (j .gt. npsiCDBM .and. j .le. npsiCT+2) then
            elseif (j .gt. npsiCDBMmin .and. j .le. npsihm) then
               if (acoef(3207) .gt. 0) then
                  jmin_cdbm = INT(0.5*(acoef(3203)+acoef(3204))*AREAL(npsit))
                  xchi = float(j-jmin_cdbm)/float(npsiCT-npsiCDBM)
                  xchi = 0.5*(TANH(xchi)+1.0)
               else
                  xchi = float(j-npsiCDBM)/float(npsiCT-npsiCDBM)
               endif
!              chiimks = (chicdbm(j)+chiinc)*(1._R8-xchi)               &
               chiicdbms(j) = (chicdbm(j)+chiinc)*(1._R8-xchi)          &
     &                 + chiimksCT*xchi
               chiecdbms(j) = (chicdbm(j)+chiinc)*(1._R8-xchi)          &
     &                 + chiemksCT*xchi
                chiimks = (chiicdbms(j-4)+chiicdbms(j-3)+chiicdbms(j-2)+&
     &                     chiicdbms(j-1)+chiicdbms(j))/5.0
                chiemks = (chiecdbms(j-4)+chiecdbms(j-3)+chiecdbms(j-2)+&
     &                     chiecdbms(j-1)+chiecdbms(j))/5.0
!               chiemks = chiimks
            elseif (j .gt. npsihm) then
                chiicdbms(j) = chiimksCT
                chiecdbms(j) = chiemksCT
                chiemks = chiemksCT
                chiimks = chiimksCT
            end if
         CASE(17)
            jmin_cdbm = 1
            qmin_CT = qprof2(3)
            do kk=3,npsit
               if (qprof2(kk) .le. qmin_CT) then
                  qmin_CT = qprof2(kk)
                  jmin_cdbm = kk
               endif
            enddo
            xchi = secitb
            if (secitb .gt. 0.) xchi = float(jmin_cdbm)/float(npsit)
            if (acoef(3205) .gt. 0._R8)                                 &
     &         chiped=acoef(3205)*qprof2(2)/qprof2(jmin_cdbm)
            call coppi ( chiemks, chiimks, chiauxs, chiohms,            &
     &       an0,anel,xmag,tfluxb,tflux,ptot,formfs,formis,vprime,gps,  &
     &       q95pct,arad,gzero,zeff,alphar,acoef(121),acoef(122),qadd,  &
     &       acoef(124),acoef(126),acoef(3015),tfluxs,lsaw,dpsi,fhmodei,&
     &       pwidthc,firitb,xchi,chiped,npsihm,global(132))

      END SELECT

!
      if(use_user_chie .and. ex_chie) then
         do  kk = 1,ntime_chie-1
            if(times .gt. time_chie(kk) .and. times .le. time_chie(kk+1)) then
              j1 = kk
              j2 = kk + 1
              tint = (times-time_chie(j1))/(time_chie(j2)-time_chie(j1))
              exit
            endif
         enddo
         if(times .le. time_chie(1)) then
           j1 = 1
           j2 = 2
           tint = 0.
         endif
         if(times .gt. time_chie(ntime_chie)) then
           j1 = ntime_chie-1
           j2 = ntime_chie
           tint = 1.
         endif
         tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
         do ii=1,nrad_chie-1
            if(tfluxd .gt. rad_chie(ii) .and. tfluxd .le. rad_chie(ii+1)) then
              rint = (tfluxd-rad_chie(ii))/(rad_chie(ii+1)-rad_chie(ii))
              pchi1 = (prof_chie(ii,j1)                                &
     &               + (prof_chie(ii+1,j1)-prof_chie(ii,j1))*rint)
              pchi2 = (prof_chie(ii,j2)                                &
     &               + (prof_chie(ii+1,j2)-prof_chie(ii,j2))*rint)
              chie_read = pchi1 + (pchi2-pchi1)*tint
              exit
            endif
         enddo             
         if(tfluxd .le. rad_chie(1)) then
           pchi1 = prof_chie(1,j1)
           pchi2 = prof_chie(1,j2)
           chie_read = pchi1 + (pchi2-pchi1)*tint
         endif
         if (acoef(3210) .eq. 0.0) then
            chiemks = chie_read
         else
            npsiCDBM = int(acoef(3204)*AREAL(npsit))
            npsiCDBMmin = int(0.5*acoef(3204)*AREAL(npsit))
            npsiCT = int(acoef(3203)*AREAL(npsit))
            if (npsiCDBM .eq. npsiCT) npsiCDBM = npsiCDBM-3
            xchi = 0._R8
            if (j .le. npsiCDBMmin) then
               chiemks = chie_read
            elseif (j .gt. npsiCDBMmin .and. j .le. npsihm) then
               if (acoef(3207) .gt. 0) then
                 jmin_cdbm = INT(0.5*(acoef(3203)+acoef(3204))*AREAL(npsit))
                 xchi = float(j-jmin_cdbm)/float(npsiCT-npsiCDBM)
                 xchi = 0.5*(TANH(xchi)+1.0)
               else
                 xchi = float(j-npsiCDBM)/float(npsiCT-npsiCDBM)
               endif
               chiemks = chie_read*(1._R8-xchi)+chiemksCT*xchi
            elseif (j .gt. npsihm) then
                chiemks = chiemksCT
            end if
         endif
         chiimks = chiemks*abs(acoef(126))
      end if
!
      go to 6776

6776  CONTINUE
!
!.....particle diffusion
      floatj = (real(j)-1._R8+1.E-6_R8)/(real(npsit)-1._R8)
      dperp = acoef(851) + acoef(875)*floatj**2
      if(floatj .ge. 0.75_R8) dperp = acoef(876)
      dperpa(j) = dperp
      coef1 = dperp*udst*fbchi*gps/rmid
!
      chiecopi(j) = chiemks*fbchi
      coef2 =      -chiemks*fbchi*udst*gps
!
      chiicopi(j) = chiimks*fbchi
      coef3 =      -chiimks*fbchi*udst*gps
      if(itrmod .eq. 13 .or. itrmod .eq. 14) then
      chiicopi(j) = chiimks*fbchii
      coef3 =      -chiimks*fbchii*udst*gps
      endif

!
!
!.....increase transport for balloon unstable surfaces if ibalsw=2
      if(ibalsw.eq.2) then
        balfac = 1._R8
        if(j+1.gt.3.and.idn(j-2).gt.0)     balfac = balfac+.040_R8*      &  
     & idn(j-2)
        if(j.gt.1.and.idn(j-1).gt.0)       balfac = balfac+.080_R8*      &  
     & idn(j-1)
        if(idn(j).gt.0)                    balfac = balfac+.240_R8*      &  
     & idn(j)
        if(idn(j+1  ).gt.0)                balfac = balfac+.080_R8*      &  
     & idn(j+1)
        if(j+1.lt.npsit.and.idn(j+2).gt.0) balfac = balfac+.040_R8*      &  
     & idn(j+2)
              coef2 = coef2*balfac
              coef3 = coef3*balfac
      endif
!
!.....add GLF23 for itrmode=8
      if(itrmode.eq.8 .or. itrmod .eq. 14) go to 6
!.....add MMM for itrmode=10
      if(itrmod.eq.10) go to 10
!
!
   20 continue
!
      cs1(j) =   coef1
      if(idens.eq.0) cs0(j) = coef1*rmid*acoef(852)/((npsit-1)*dpsi)
      if(idens.eq.3) cs0(j) = acoef(852)*udst*fbchi*sqrt(gps)
!
      dse0(j) = -face*cs0(j)
      dse1(j) = -coef2*emid/rmid - face*cs1(j)
      dse3(j) =  coef2
!
      dsi0(j) = -faci*cs0(j)
      dsi1(j) = -coef3*(pmid-emid)/rmid - faci*cs1(j)
      dsi2(j) =  coef3
      dsi3(j) = -coef3
      go to 4000
!
!.........................................................................
!.....coppi/mazzucato/gruber transport model (L.Sugiyama)
!.........................................................................
    3 befo = 1._R8
      torcur = gxmja2(j)/(usdi*tpi)
      chicmg1 = 1.5E13_R8*acoef(121)*torcur
      chicmg2 = ((.5_R8*(ane(j+1)+ane(j)))**0.8_R8)*te(j)
      chicmg3 = ((allam/2.4_R8)**0.4_R8)*((vary(npsit))**2)/(vp2(j))**2
      chicmg4 = 2.0_R8*sqrt(pi)*(((vary(npsit))/(tpi*xmag))**1.5_R8)
      chicmg = (chicmg1*chicmg3)/(chicmg2*chicmg4)
!
      pchinum = palpha+paux
      pchiden = palpha+paux+pohmic
      if (pchiden.eq.0._R8) then
      pchirat = 1
      else
      pchirat = pchinum/pchiden
      endif
      go to 3987
!.........................................................................
!.....coppi/mazzucato/gruber transport model  (R.Englade)
!.........................................................................
!
    4 befo = 1._R8
      torcur = gxmja2(j)/(usdi*tpi)
      chicmg1 = 2.296E12_R8*acoef(121)*torcur
      chicmg2 = ((.5_R8*(ane(j+1)+ane(j)))**0.8_R8)*te(j)
      chicmg3 = ((vary(npsit))**2)/(vp2(j))**2
      chicmg4 = (2.4_R8**0.4_R8)*(vary(npsit)/(2*(pi**2)*xmag))**1.5_R8
      chicmg = (chicmg1*chicmg3)/(chicmg2*chicmg4)
!
      pchinum = palpha+paux
      pchiden = palpha+paux+pohmic
      if (pchiden.eq.0._R8) then
      pchirat = 1
      else
      pchirat = (pchinum/pchiden)**2
      endif
!
 3987 continue
      chi0 = chicmg*(1._R8+(acoef(122)*betapol*pchirat))/usdt
!........ increase chielectron if Te>qadd
      if (te(1).gt.qadd) chi0=chi0*te(1)/qadd
      ache0(j+1) = chi0
!.....increase transport by 10 for balloon unstable surfaces if ibalsw=2
      if(ibalsw.eq.2 .and. idn(j).gt.0) befo = befo*10._R8
!
!.....enhance electron transport for sawtooth model
!
!.....neoclassical ion transport   model
      d0neo = qamin(j+1)/xmag
      ailarr = slight*sqrt(2*aionm*timide)/(elchar*bpolar(j+1))
      allii = 23.0_R8-log(sqrt(animid)/timid)
      tauii=3._R8*sqrt(aionm)*(timide**1.5_R8)/(4._R8*sqrt(pi)           &  
     &      *allii*animid*                                               &  
     &      elchar**4)
      call geval(psimin,2,gval,gpval,gppval,imag,jmag)
      bzero = gval/xmag*1.E4_R8
      collfii=bzero*(xmag*1.E2_R8)**1.5_R8/(tauii*sqrt(qamin(j+1)*       &  
     & 1.E2_R8)*                                                         &  
     &        bpolar(j+1)*sqrt(timide/aionm))
      aki=(.66_R8+1.88_R8*sqrt(d0neo)-1.54_R8*d0neo)/(1._R8+1.03_R8*     &  
     & sqrt(collfii)+.31_R8*                                             &  
     &collfii)*(bzero**2)*bsqiar(j+1)+(.58_R8*collfii*d0neo)/(1._R8      &  
     &   +.74_R8*                                                        &  
     &collfii*d0neo**1.5_R8)*((bzero**2)*bsqiar(j+1)-(bzero**2)/bsqar(j+  &  
     & 1))
      chiineo(j)=aki*zeff*ailarr**2*sqrt(d0neo)/tauii*1.E-4_R8
!.....electron thermal conductivity in m**2/sec
      if(gps.ne.0._R8) then
      chiesecj = chi0 *befo/(udst*gps)
      chiisecj = abs(acoef(126))*chiesecj+chiineo(j)
      if(cswitchre.eq.0._R8) then
      chiemin = chiesecj
      cswitchre = 1._R8
      endif
!     if (chiesecj.lt.chiemin) then
!     do 3998 k=2,j
!     ache0(k) = ache0(k)*chiesecj/chiesec(k)
!     chiesec(k) = chiesecj
!3998 continue
!     endif
      endif
!
!......no convection,  factor of abs(acoef(126)) less ion conduction
      coef1 = 0._R8
      coef2 = -befo*chi0
      coef3 = -chiisecj*udst*gps
!     coef3 = chicmg/usdt*(abs(acoef(126))+acoef(122)*betapol*
!    1 pchirat*(abs(acoef(126))+acoef(125))) *befo
!
!
      cs1(j) = coef1
!
      dse1(j) = -coef2*emid/rmid - face*cs1(j)
      dse3(j) = coef2
!
      dsi1(j) = -coef3*(pmid-emid)/rmid - faci*cs1(j)
      dsi2(j) = coef3
      dsi3(j) =-coef3
      aemid(j+1)=emid
      armid(j+1)=rmid
      aface(j+1)=face
!
      if (chiesecj.lt.chiemin) then
      chiemin = chiesecj
      do 3989 k=2,j
      coef2 = -befo*ache0(k)
      dse1(k-1) = -coef2*aemid(k)/armid(k)-aface(k)*cs1(k-1)
      dse3(k-1) = coef2
 3989 continue
      endif
      go to 4000
!
!
!
!.........................................................................
!.....Alcator Electron Transport + Neo-classical Ion Transport
!.....M.Turner, J.Connor.     April 1991
!.........................................................................
!
    5 continue
!
!.....enhance electron transport for sawtooth model
!
      gpsd = (gps/usdt)*1.E19_R8/(udsd*rmid)
!...
      coef1 = acoef(39)*gpsd/rmid
      coef2 =-acoef(35)*gpsd
!
!
!.....neoclassical ion transport   model
      d0neo = qamin(j+1)/xmag
      alarsq = 2.073_R8*amgas*timid/(bpolar(j+1)**2)
      allii = 23._R8-log(sqrt(animid)/timid)
      animd1 = animid*1.0E-11_R8
      tauii=6.57_R8*sqrt(amgas)*((timid*1.0E-3_R8)**1.5_R8)/             &  
     &(allii*animd1)
      call geval(psimin,2,gval,gpval,gppval,imag,jmag)
      bzero = gval/xmag*1.E4_R8
      tcoli=2.886_R8*(1.0E-4_R8*bpolar(j+1)/bzero)*d0neo*(timid**2)/     &  
     &(bzero*xmag*allii*animd1)
      alfnc = zeff - 1._R8
      tmsti = tcoli/(1.0_R8+ 1.54_R8*alfnc)
      aki1=(0.66_R8*(1.0_R8+1.54_R8*alfnc)+(1.88_R8*sqrt(d0neo)-1.54_R8*  &  
     & d0neo)*                                                           &  
     &(1.0_R8+1.375_R8*alfnc))*tmsti*(bzero**2)*bsqiar(j+1)/             &  
     &(tmsti+1.03_R8*sqrt(tmsti)+0.31_R8)
      aki2=(0.58_R8*d0neo)*(1.0_R8+1.33_R8*alfnc*(1.0_R8+0.6_R8*alfnc)/  &  
     &(1.0_R8+1.79_R8*alfnc))*                                           &  
     &(bzero**2)*(bsqiar(j+1)-1.0_R8/bsqar(j+1))/(tmsti+.74_R8*d0neo**   &  
     & 1.5_R8)
      chiineo(j)=(aki1+aki2)*zeff*alarsq*sqrt(d0neo)/tauii
!
!.....thermal conductivity in m**2/sec
      if(gps.ne.0) then
      chiesecj = -coef2/(udst*gps)
      chiisecj =  acoef(126)*chiesecj + chiineo(j)
      endif
      coef3 = -chiisecj*udst*gps
!
!
      cs1(j) = coef1
!
      dse1(j) = -coef2*emid/rmid - face*cs1(j)
      dse3(j) = coef2
!
      dsi1(j) = -coef3*(pmid-emid)/rmid - faci*cs1(j)
      dsi2(j) = coef3
      dsi3(j) =-coef3
      go to 4000
!
!
!.................................................................
!..... glf2D transport model
!.................................................................
!
    6 continue
!
!.....only call glf23 every nskipsf cycles to save time
      if(iskipsf.ne.1) go to 9
!
      if(j.ne.1) go to 106
!
      te_sav = te
      ti_sav = ti
!
!     smoothing functions
      if(1 .eq. 2)then
!     if(j .ge. npsihm-7) then
!     jstart=4
      jstart = npsihm-10
      do i=jstart, npsit-1
      te_sav(i)=(te(i-2)+te(i-1)+te(i)+te(i+1)+te(i+2))/5.0
      ti_sav(i)=(ti(i-2)+ti(i-1)+ti(i)+ti(i+1)+ti(i+2))/5.0
      enddo
      te_sav(npsit)=2.*te_sav(npsit-1)-te_sav(npsit-2)
      ti_sav(npsit)=2.*ti_sav(npsit-1)-ti_sav(npsit-2)
      if( jstart .eq. 4 ) then
      te_sav(3)=(te(2)+te(3)+te(4))/3.0
      ti_sav(3)=(ti(2)+ti(3)+ti(4))/3.0
      te_sav(2)=te(2)
      ti_sav(2)=ti(2)
!     te_sav(2)=(te(1)+te(2)+te(3))/3.0
!     ti_sav(2)=(ti(1)+ti(2)+ti(3))/3.0
!     te_sav(2)=2.*te_sav(3)-te_sav(4)
!     ti_sav(2)=2.*ti_sav(3)-ti_sav(4)
!     te_sav(2)=1.5*te_sav(3)-0.5*te_sav(4)
!     ti_sav(2)=1.5*ti_sav(3)-0.5*ti_sav(4)
!     te_sav(2)=3.*te_sav(3)-3.*te_sav(4)+te_sav(5)
!     ti_sav(2)=3.*ti_sav(3)-3.*ti_sav(4)+ti_sav(5)
      te_sav(1)=te_sav(2)
      ti_sav(1)=ti_sav(2)
      endif
      endif

      if(use_user_rot .and. first_read_rot) then
      inquire(file="user_rot_data", exist=ex_rot)
      if (ex_rot) then
      open(55,file="user_rot_data",form="formatted",status="old")
      read(55,*) ntime_rot, nrad_rot
      allocate(time_rot(ntime_rot), rad_rot(nrad_rot))
      allocate(wrot(nrad_rot,ntime_rot))
      read(55,*) time_rot(1:ntime_rot)
      read(55,*) rad_rot(1:nrad_rot)
      read(55,*) wrot(1:nrad_rot,1:ntime_rot)
      close (55)
      endif
      first_read_rot=.false.
      endif

      if(use_user_rot .and. ex_rot) then
      do  kk = 1,ntime_rot-1
      if(times .gt. time_rot(kk) .and. times .le. time_rot(kk+1)) then
      j1 = kk
      j2 = kk + 1
      tint = (times-time_rot(j1))/(time_rot(j2)-time_rot(j1))
      exit
      endif
      enddo
      if(times .le. time_rot(1)) then
      j1 = 1
      j2 = 2
      tint = 0.
      endif 
      if(times .gt. time_rot(ntime_rot)) then
      j1 = ntime_rot-1
      j2 = ntime_rot
      tint = 1.
      endif
      endif

      do i=0,npsit-1
      te_m(i) = 0.5_R8*(te(i+1)+te(i+2))*.001_R8
      ti_m(i) = 0.5_R8*(ti(i+1)+ti(i+2))*.001_R8
      if(ti_m(i) .le. 0.1_R8*te_m(i)) ti_m(i) = 0.1_R8*te_m(i)
      ne_m(i) = 0.5_R8*(ane(i+1)+ane(i+2))/1.E19_R8
      ni_m(i) = 0.5_R8*(ane(i+1)/avez(i+1)+ane(i+2)/avez(i+2))/1.E19_R8
      ns_m(i) = 0._R8
!
!.....toroidal angular velocity in m/s
!     angrotp_exp(i) = 0.0d3*float(npsit-i)/float(npsit)
      angrotp_exp(i) = acoef(4986)*float(npsit-i)/float(npsit)

      if(use_user_rot .and. ex_rot) then
      tfluxd = sqrt((float(i-1)*dpsi)/(float(npsit-1)*dpsi))
      do ii=1,nrad_rot
      if(tfluxd .gt. rad_rot(ii) .and. tfluxd .le. rad_rot(ii+1)) then
      rint = (tfluxd-rad_rot(ii))/(rad_rot(ii+1)-rad_rot(ii))
      wrot1 = wrot(ii,j1) + (wrot(ii+1,j1)-wrot(ii,j1))*rint
      wrot2 = wrot(ii,j2) + (wrot(ii+1,j2)-wrot(ii,j2))*rint
      frot = (wrot1 + (wrot2-wrot1)*tint)/tpi
      exit
      endif
      enddo
      angrotp_exp(i) = frot
      endif

      egamma_exp(i)  = 0._R8
      gamma_p_exp(i) = 0._R8
!
!.....toroidal velocity in m/s
      vphi_m(i) = 0._R8
      vpar_m(i) = 0._R8
      vper_m(i) = 0._R8
      zeff_exp(i) = 0.5_R8*(zeffa(i+1)+zeffa(i+2))
      if(i.ne.0) then
      rho(i) = sqrt(AREAL(i)/(npsit-1._R8))
      rmin_exp(i) =rminora(i+1)
      grs = gja2(i+1)*(tpi*qprof2(i+1))/vp2(i+1)
!
!....must verify these are order unity
      gradrho_exp(i) = .5_R8*sqrt(grs*xplas/(pi*gzero*i*dpsi))
      gradrhosq_exp(i) = gradrho_exp(i)*gradrho_exp(i)
      shat_exp(i) = (i+1)*(qprof2(i+2)-qprof2(i))/qprof2(i+1)
      else
      rho(0) = 0._R8
      rmin_exp(i) = 0._R8
      gradrho_exp(0) = 0._R8
      gradrhosq_exp(0) = 0._R8
      shat_exp(0) = 0._R8
      endif
      q_exp(i) = qprof2(i+1)
!
!.....so-called alpha effect not turned on
!     alpha_exp(i) = 0._R8
      alpha_exp(i) = acoef(4987)
      elong_exp(i) = elonga(i+1)
      rmaj_exp(i) = rmajora(i+1)
      enddo
      do ii=1,npsit-1
      i = npsit-ii
      if(te_m(i-1).eq.te_m(i))te_m(i-1)=1.01_R8*te_m(i)
      if(ti_m(i-1).eq.ti_m(i))ti_m(i-1)=1.01_R8*ti_m(i)
      enddo

      te_m_sav = te_m
      ti_m_sav = ti_m

  106 continue

!..Bateman, 15 March 2009, change to jshoot = 1
!  to avoid shift to zone boundaries made in callglf2d

!      jshoot = 0

      jshoot = 1

      jmm = j-1
      jmaxm = npsit
!     itport_pt(1) = 0
!     itport_pt(2) = 1
!     itport_pt(3) = 1
!     itport_pt(4) = 0
!     itport_pt(5) = 0
      itport_pt(1)=int(acoef(4981))
      itport_pt(2)=int(acoef(4982))
      itport_pt(3)=int(acoef(4983))
      itport_pt(4)=int(acoef(4984))
      itport_pt(5)=int(acoef(4985))
!
      i_grad = 1
!
!..Changes made by Bateman 12 March 2009 suggested by Long-Poe Ku

       facs = 2._R8*sqrt(AREAL(j-1)/AREAL(npsit-1))

       zpte_in =-2.*facs*float(npsit-1)*                                  &
      &              (te(j+1)-te(j))/(te(j+1)+te(j))
       zpti_in =-2.*facs*float(npsit-1)*                                  &
      &              (ti(j+1)-ti(j))/(ti(j+1)+ti(j))
       if(zpti_in .le. 0) zpti_in = 0.1_R8*zpte_in
       zpne_in =-2._R8*facs*(npsit-1)*(ane(j+1)-ane(j))/(ane(j+1)+ane(j))
       zpni_in =-2._R8*facs*(npsit-1)*(ane(j+1)/avez(j+1)-ane(j)/avez(j))/ &
      &  (ane(j+1)/avez(j+1)+ane(j)/avez(j))

! WRITE(*,*)' j. npsit zpti_in facs ti(j+1) ti(j)'
! WRITE(*,"(2I5,9ES12.4)") j, npsit, zpti_in &
!   , facs, ti(j+1), ti(j)

      bt_exp = gzero/xplas
      arho_exp = sqrt(tfluxb/(pi*bt_exp))
      rmajor_exp = xmag
      amassgas_exp = amgas
      zimp_exp = zimp
      amassimp_exp = 2*zimp

!  Bateman, 27 March 2009, restore alpha_e = 1.35_R8

      alpha_e = 1.35_R8

      x_alpha = 0._R8
      i_delay = 0
      leigen = 0
      idengrad = 2
      nroot = 8
      iglf = 1
!     irotstab = 0
      irotstab=int(acoef(4988))
      ibt_flag = 1
!
      npmax = npsit
      if(itrmode.eq.8 .or. itrmode.eq.9) npmax = npsihm
      if(itrmode.eq.6) npmax = npsihm
      chimin_tmp = chimin
      if(j .le. 5) chimin_tmp=1.*chimin

IF ( iprint > 0 ) WRITE(*,*) ' npmax = ', npmax

      if(j.gt.1 .and. j.le.npmax) then

!..Diagnostic printout by Bateman

IF ( iprint > 0 ) THEN

  WRITE(*,*) &
   'j, rmin_exp(j), ti(j+1), ti(j), ti_m, te(j+1), te(j), te_m'

  WRITE(*,"(' 1# ',I5,9ES12.4)") &
    j, rmin_exp(j), ti(j+1), ti(j), ti_m(j), te(j+1), te(j), te_m(j)

ENDIF
!

       call callglf2d(leigen,nroot,iglf,                                 &  
     & jshoot,jmm, jmaxm, itport_pt,irotstab,te_m,ti_m,ne_m,ni_m,ns_m,   &  
     & i_grad, idengrad,                                                 &  
     & zpte_in, zpti_in, zpne_in,zpni_in,angrotp_exp,egamma_exp,         &  
     & gamma_p_exp,vphi_m,vpar_m,vper_m,zeff_exp,bt_exp,ibt_flag,        &  
     & rho, arho_exp,                                                    &  
     & gradrho_exp,gradrhosq_exp,rmin_exp,rmaj_exp, rmajor_exp,          &  
     & zimp_exp, amassimp_exp, q_exp,                                    &  
     & shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, x_alpha,   &  
     & i_delay,                                                          &  
!    &                                                                   &  
!    &                !OUTPUTS                                           &
     & diffnem,chietem,chiitim,etaphim, etaparm,etaperm, exchm, diff_m,  &  
     & chie_m, chii_m, etaphi_m, etapar_m, etaper_m, exch_m, egamma_m,   &  
     & egamma_d, gamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m )
      nprint = 0.8_R8*npsit
      if(j.ge.nprint .and. (zpte_in .le. 0._R8 .or. zpti_in .le. 0._R8))  &  
     & write(nterm,6666)" diag6666", kcycle, j, zpte_in,zpti_in,chietem,  &  
!    &                                                                    &  
     &                  chiitim
 6666 format(a10,2i6,1p4e12.4)
!
!..Diagnostic printout by Bateman

      IF ( iprint > 0 ) THEN

        WRITE(*,"(' 1# ',I5,9ES12.4)") &
          j, chiitim, chietem, diffnem, rmin_exp(j-1)

      ENDIF
!
!--->end of diag
!
!.....put limits on chietem and chiitim
      if(chiitim .lt. chimin_tmp) chiitim = chimin_tmp
      if(chietem .lt. chimin_tmp) chietem = chimin_tmp
      if(diffnem .lt. chimin_tmp) diffnem = chimin_tmp
      if(chiitim .gt. chimax)  chiitim = chimax
      if(chietem .gt. chimax)  chietem = chimax
      if(diffnem .gt. chimax)  diffnem = chimax
!
      zchii(j,1) = chiitim
!
!.....apply relaxation factor to deal with stiff equation
      relfactor = acoef(3010)
      chiitima(j) = (1._R8-relfactor)*chiitima(j) + relfactor*1.0_R8*    &  
     & chiitim
      chietema(j) = (1._R8-relfactor)*chietema(j) + relfactor*1.0_R8*    &  
     & chietem
      diffnema(j) = (1._R8-relfactor)*diffnema(j) + relfactor*diffnem
!     cs1(j) = (diffnem/rmid)*gps*udst
!     [ note factor of 1.5 in going from ITER to TSC convention ]
!     [...removed 06/25/02]

!
!...linearization of Chi
!
!      delta_t = -1.0e-06
       delta_t=-1.0
       if(acoef(4958) .gt. 0.0 ) delta_t = acoef(4958)
       dchie_dtip(j) = 0.0
       dchii_dtip(j) = 0.0
       dchie_dtep(j) = 0.0
       dchii_dtep(j) = 0.0
       dti_dpsi(j)=0.0
       dte_dpsi(j)=0.0
       if(delta_t .gt. 0.0 .and. acoef(4957) .gt. 0.0) then
          te_m = te_m_sav
          t_j=te_sav(j)*(1+delta_t)
          t_j1=te_sav(j+1)*(1-delta_t)
          temp1=-2.*facs*float(npsit-1)*(t_j1-t_j)/(t_j1+t_j)

          te_m(jmm)=0.5*(t_j1+t_j)*0.001
          te_m(jmm-1)=0.5*(t_j+te_sav(j-1))*0.001
          te_m(jmm+1)=0.5*(t_j1+te_sav(j+2))*0.001

ii = 2

          call callglf2d(leigen,nroot,iglf,                                &
     &    jshoot,jmm, jmaxm, itport_pt,irotstab,te_m,ti_m,ne_m,ni_m,ns_m,  &
     &    i_grad, idengrad,                                                &
     &    temp1, zpti_in, zpne_in,zpni_in,angrotp_exp,egamma_exp,          &
     &    gamma_p_exp,vphi_m,vpar_m,vper_m,zeff_exp,bt_exp,ibt_flag,       &
     &    rho, arho_exp,                                                   &
     &    gradrho_exp,gradrhosq_exp,rmin_exp,rmaj_exp, rmajor_exp,         &
     &    zimp_exp, amassimp_exp, q_exp,                                   &
     &    shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, x_alpha,  &
     &    i_delay,                                                         &
!    &                 !OUTPUTS                                         &
     &    temp21,temp31,temp41,etaphim, etaparm,etaperm, exchm, diff_m,    &
     &    chie_m, chii_m, etaphi_m, etapar_m, etaper_m, exch_m, egamma_m,  &
     &    egamma_d, gamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m )

          te_m = te_m_sav
          t_j=te_sav(j)*(1-delta_t)
          t_j1=te_sav(j+1)*(1+delta_t)
          temp2=-2.*facs*float(npsit-1)*(t_j1-t_j)/(t_j1+t_j)

          te_m(jmm)=0.5*(t_j1+t_j)*0.001
          te_m(jmm-1)=0.5*(t_j+te_sav(j-1))*0.001
          te_m(jmm+1)=0.5*(t_j1+te_sav(j+2))*0.001

       call callglf2d(leigen,nroot,iglf,                                &
     &    jshoot,jmm, jmaxm, itport_pt,irotstab,te_m,ti_m,ne_m,ni_m,ns_m,  &
     &    i_grad, idengrad,                                                &
     &    temp2, zpti_in, zpne_in,zpni_in,angrotp_exp,egamma_exp,          &
     &    gamma_p_exp,vphi_m,vpar_m,vper_m,zeff_exp,bt_exp,ibt_flag,       &
     &    rho, arho_exp,                                                   &
     &    gradrho_exp,gradrhosq_exp,rmin_exp,rmaj_exp, rmajor_exp,         &
     &    zimp_exp, amassimp_exp, q_exp,                                   &
     &    shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, x_alpha,  &
     &    i_delay,                                                         &
!    &                 !OUTPUTS                                         &
     &    temp22,temp32,temp42,etaphim, etaparm,etaperm, exchm, diff_m,    &
     &    chie_m, chii_m, etaphi_m, etapar_m, etaper_m, exch_m, egamma_m,  &
     &    egamma_d, gamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m )   

         if(temp31 .lt. chimin_tmp) temp31 = chimin_tmp
         if(temp32 .lt. chimin_tmp) temp32 = chimin_tmp
         if(temp41 .lt. chimin_tmp) temp41 = chimin_tmp
         if(temp42 .lt. chimin_tmp) temp42 = chimin_tmp
         if(temp31 .gt. chimax)  temp31 = chimax
         if(temp32 .gt. chimax)  temp32 = chimax
         if(temp41 .gt. chimax)  temp41 = chimax
         if(temp42 .gt. chimax)  temp42 = chimax

         zchii(j,2) = temp41
         zchii(j,3) = temp42

         te_m = te_m_sav
         dchie_dtep(j) = (temp31-temp32)/(temp1-temp2)
         dchie_dtep(j) = -facs*2.0/(te_sav(j)+te_sav(j+1))*dchie_dtep(j)
         dchie_dtep(j) = dchie_dtep(j)*dpsi*(npsit-1)
         dchii_dtep(j) = (temp41-temp42)/(temp1-temp2)
         dchii_dtep(j) = -facs*2.0/(te_sav(j)+te_sav(j+1))*dchii_dtep(j)
         dchii_dtep(j) = dchii_dtep(j)*dpsi*(npsit-1)


         ti_m = ti_m_sav
         t_j=ti_sav(j)*(1+delta_t)
         t_j1=ti_sav(j+1)*(1-delta_t)
         temp1=-2.*facs*float(npsit-1)*(t_j1-t_j)/(t_j1+t_j)

         ti_m(jmm)=0.5*(t_j1+t_j)*0.001
         ti_m(jmm-1)=0.5*(t_j+ti_sav(j-1))*0.001
         ti_m(jmm+1)=0.5*(t_j1+ti_sav(j+2))*0.001

         call callglf2d(leigen,nroot,iglf,                                &
     &   jshoot,jmm, jmaxm, itport_pt,irotstab,te_m,ti_m,ne_m,ni_m,ns_m,  &
     &   i_grad, idengrad,                                                &
     &   zpte_in, temp1, zpne_in,zpni_in,angrotp_exp,egamma_exp,          &
     &   gamma_p_exp,vphi_m,vpar_m,vper_m,zeff_exp,bt_exp,ibt_flag,       &
     &   rho, arho_exp,                                                   &
     &   gradrho_exp,gradrhosq_exp,rmin_exp,rmaj_exp, rmajor_exp,         &
     &   zimp_exp, amassimp_exp, q_exp,                                   &
     &   shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, x_alpha,  &
     &   i_delay,                                                         &
!    &                 !OUTPUTS                                         &
     &   temp21,temp31,temp41,etaphim, etaparm,etaperm, exchm, diff_m,    &
     &   chie_m, chii_m, etaphi_m, etapar_m, etaper_m, exch_m, egamma_m,  &
     &   egamma_d, gamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m )   

         ti_m = ti_m_sav
         t_j=ti_sav(j)*(1-delta_t)
         t_j1=ti_sav(j+1)*(1+delta_t)
         temp2=-2.*facs*float(npsit-1)*(t_j1-t_j)/(t_j1+t_j)

         ti_m(jmm)=0.5*(t_j1+t_j)*0.001
         ti_m(jmm-1)=0.5*(t_j+ti_sav(j-1))*0.001
         ti_m(jmm+1)=0.5*(t_j1+ti_sav(j+2))*0.001
         call callglf2d(leigen,nroot,iglf,                                &
     &   jshoot,jmm, jmaxm, itport_pt,irotstab,te_m,ti_m,ne_m,ni_m,ns_m,  &
     &   i_grad, idengrad,                                                &
     &   zpte_in, temp2, zpne_in,zpni_in,angrotp_exp,egamma_exp,          &
     &   gamma_p_exp,vphi_m,vpar_m,vper_m,zeff_exp,bt_exp,ibt_flag,       &
     &   rho, arho_exp,                                                   &
     &   gradrho_exp,gradrhosq_exp,rmin_exp,rmaj_exp, rmajor_exp,         &
     &   zimp_exp, amassimp_exp, q_exp,                                   &
     &   shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, x_alpha,  &
     &   i_delay,                                                         &
!    &                 !OUTPUTS                                         &
     &   temp22,temp32,temp42,etaphim, etaparm,etaperm, exchm, diff_m,    &
     &   chie_m, chii_m, etaphi_m, etapar_m, etaper_m, exch_m, egamma_m,  &
     &   egamma_d, gamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m )

        if(temp31 .lt. chimin_tmp) temp31 = chimin_tmp
        if(temp32 .lt. chimin_tmp) temp32 = chimin_tmp
        if(temp41 .lt. chimin_tmp) temp41 = chimin_tmp
        if(temp42 .lt. chimin_tmp) temp42 = chimin_tmp
        if(temp31 .gt. chimax)  temp31 = chimax
        if(temp32 .gt. chimax)  temp32 = chimax
        if(temp41 .gt. chimax)  temp41 = chimax
        if(temp42 .gt. chimax)  temp42 = chimax

         zchii(j,4) = temp41
         zchii(j,5) = temp42

        ti_m = ti_m_sav
        dchie_dtip(j) = (temp31-temp32)/(temp1-temp2)
        dchie_dtip(j) = -facs*2.0/(ti_sav(j)+ti_sav(j+1))*dchie_dtip(j)
        dchie_dtip(j) = dchie_dtip(j)*dpsi*(npsit-1)
        dchii_dtip(j) = (temp41-temp42)/(temp1-temp2)
        dchii_dtip(j) = -facs*2.0/(ti_sav(j)+ti_sav(j+1))*dchii_dtip(j)
        dchii_dtip(j) = dchii_dtip(j)*dpsi*(npsit-1)
        dti_dpsi(j) = (ti_sav(j+1)-ti_sav(j))/dpsi
        dte_dpsi(j) = (te_sav(j+1)-te_sav(j))/dpsi

      endif

      else
      chiitima(j) = 1.E-8_R8
      chietema(j) = 1.E-8_R8
      diffnema(j) = 1.E-8_R8
      dchie_dtip(j) = 0.0
      dchii_dtip(j) = 0.0
      dchie_dtep(j) = 0.0
      dchii_dtep(j) = 0.0
      dti_dpsi(j)=0.0
      dte_dpsi(j)=0.0
      endif
!
!     
!...connect experimental chi fro d3dcoppi to GLF
!     
      tflux = float(j-1)*dpsi
      if(itrmod .eq. 14 .and. pwidthc .gt. 0.                       &
     & .and. tflux .ge. npsihm*dpsi) then
      chiitima(j) = chiicopi(j)
      chietema(j) = chiecopi(j)
      endif
!
!...pedestal modification to GLF23 if itrmode is 6 or 8
!
      tflux = float(j-1)*dpsi
      if(itrmod .eq. 6 .or. itrmod .eq. 8) then
      if(chiped .ne. 0. .and. pwidthc                          &
     & .ne. 0. .and. tflux .ge. npsihm*dpsi) then
      chiitima(j) = chiped*                                        &
     &(tflux/(float(npsihm-1)*dpsi))**4
      chietema(j) = chiped*                                        &
     &(tflux/(float(npsihm-1)*dpsi))**4
!     chiitima(j) = chiitima(j)*0.1 + chiitimao(j)*0.9
!     chietema(j) = chietema(j)*0.1 + chietemao(j)*0.9
       dchie_dtip(j) = 0.0
       dchii_dtip(j) = 0.0
       dchie_dtep(j) = 0.0
       dchii_dtep(j) = 0.0
      dti_dpsi(j)=0.0
      dte_dpsi(j)=0.0
      endif
      endif
!
!     if(j .ge. npsihm-7)then
!     chiitima(j) = chiitima(j)*0.9 + chiitim*0.1
!     chietema(j) = chietema(j)*0.9 + chietem*0.1
!     endif

    9 continue
      if(j .gt. 1) then
      relfactor2 = acoef(3009)
      chiitima(j) = (1._R8-relfactor2)*chiitima(j)                       &  
     &      + .5_R8*relfactor2*(chiitimao(j+1)+chiitimao(j-1))
      chietema(j) = (1._R8-relfactor2)*chietema(j)                       &  
     &      + .5_R8*relfactor2*(chietemao(j+1)+chietemao(j-1))
      diffnema(j) = (1._R8-relfactor2)*diffnema(j)                       &  
     &      + .5_R8*relfactor2*(diffnemao(j+1)+diffnemao(j-1))
      endif
 
      if(chiitima(j) .lt. 1.E-8_R8) chiitima(j) = 1.E-8_R8
      if(chietema(j) .lt. 1.E-8_R8) chietema(j) = 1.E-8_R8
      if(diffnema(j) .lt. 1.E-8_R8) diffnema(j) = 1.E-8_R8
      if(chiitima(j) .gt. chimax)  chiitima(j) = chimax
      if(chietema(j) .gt. chimax)  chietema(j) = chimax
      if(diffnema(j) .gt. chimax)  diffnema(j) = chimax
      if(itrmode.eq.8) go to 8
      go to 7
!
!.......................................................................
!     multimode model MMM95
!.......................................................................
 10   continue
!.....only call mmm95 every nskipsf cycles to save time
      if(iskipsf.ne.1) go to 19
!
      matdim = 5
      nprout = nout
      lprint = 0
      lsuper = 0
      lreset = 0
      nerr = 0
      npoints = 1
      rminormm = rminora(j)
      rmajormm = rmajora(j)
      elongmm  = elonga(j)
      densemm = .5_R8*(ane(j)+ane(j+1))
      denshmm = .5_R8*(anhy(j)+anhy(j+1))
      densimpmm = .5_R8*(animp(j)+animp(j+1))
      densfemm = 0._R8
      zeffmm = .5_R8*(zeffa(j)+zeffa(j+1))
      tekevmm= .5_R8*(te(j)+te(j+1))*1.E-3_R8
      tikevmm= .5_R8*(ti(j)+ti(j+1))*1.E-3_R8
      qmm = qprof2(j)
      call geval(xsv2(j),2,gvalmm,gpvalmm,gppvalmm,imag,jmag)
      btormm      = gvalmm/rmajormm
      avezimpmm   =  .5_R8*(avezimpa(j)+avezimpa(j+1))
      amassimpmm  = .5_R8*(amassimpa(j)+amassimpa(j+1))
      amasshydmm  =  .5_R8*(amasshyda(j)+amasshyda(j+1))
      aimassmm    = .5_R8*(aimassa(j)+aimassa(j+1))
      wexbsmm     = 0._R8
      befogrd     = -4._R8*rmajormm/(rminora(j+1)-rminora(j-1))
      grdnemm     = befogrd*(ane(j+1)-ane(j))/(ane(j)+ane(j+1))
      grdnimm  = befogrd*(sumnia(j+1)-sumnia(j))/(sumnia(j)+sumnia(j+1))    
      grdnhmm     = befogrd*(anhy(j+1)-anhy(j))/(anhy(j)+anhy(j+1))
      grdnzmm  = befogrd*(sumnqa(j+1)-sumnqa(j))/(sumnqa(j)+sumnqa(j+1))    
      grdtemm     = befogrd*(te(j+1)-te(j))/(te(j)+te(j+1))
      grdtimm     = befogrd*(ti(j+1)-ti(j))/(ti(j)+ti(j+1))
      grdqmm      = -0.25_R8*befogrd*(qprof2(j+1)-qprof2(j-1))/qprof2(j)    
!
      npmax = npsit
      if(itrmod.eq.10 .or. itrmod.eq.10) npmax = npsihm
      if(j.gt.1 .and. j.le.npmax) then
 
!
      call mmm95 (                                                       &  
     &   rminormm,  rmajormm,   elongmm                                  &  
     & , densemm,   denshmm,    densimpmm,  densfemm                     &  
     & , zeffmm,   tekevmm,    tikevmm,    qmm,       btormm             &  
     & , avezimpmm, amassimpmm, amasshydmm, aimassmm,  wexbsmm           &  
     & , grdnemm, grdnimm, grdnhmm, grdnzmm, grdtemm, grdtimm, grdqmm    &  
     & , thiig,   thdig,    theig,    thzig                              &  
     & , thirb,   thdrb,    therb,    thzrb                              &  
     & , thikb,   thdkb,    thekb,    thzkb                              &  
     & , gamma,   omega,    difthi,   velthi,  vflux                     &  
     & , matdim,  npoints,  nprout,   lprint,  nerr                      &  
     & , lsuper,  lreset,   lswitch,  cswitch, fig,    frb,     fkb)
      if(nerr .ne. 0) then
      ineg=56
      write(nout,6656) j,nerr
      write(nterm,6656) j,nerr
 6656 format(" j, nerr after call to mmm95", 2i5)
      return
      endif
!
!.....difthi(j1,j2,jz) and velthi(j1,j2,jz) matrices contain all of
!     the anomalous transport   (note factor of 3/2 in going from ITER to TSC convention)
      chiitim = 1.5_R8*difthi(1,1,1)
      chietem = 1.5_R8*difthi(3,3,1)
      diffnem =     difthi(2,2,1)
!     write(nterm,7771) kcycle, j, diffnem, thdig ,thdrb ,thdkb
!7771 format(2i5,1p4e12.4)
!
!      write(nterm,7771) kcycle, j, chiitim, chietem, rminormm,
!     1  rmajormm, elongmm, densemm, denshmm, densimpmm, densfemm,
!     2  zeffmm, tekevmm, tikevmm, qmm, btormm,
!     3  avezimpmm, amassimpmm, amasshydmm, aimassmm, wexbsmm,
!     4  grdnemm, grdnimm, grdnhmm, grdnzmm, grdtemm, grdtimm,grdqmm
! 7771 format(2i5, 1p5e12.5,/,10x,5e12.4,/,10x,5e12.4,/,10x,5e12.4,/,
!     5       10x,7e12.4)
!
!.....put limits on chietem and chiitim
      if(chiitim .lt. 1.E-8_R8) chiitim = 1.E-8_R8
      if(chietem .lt. 1.E-8_R8) chietem = 1.E-8_R8
      if(diffnem .lt. 1.E-8_R8) diffnem = 1.E-8_R8
      if(chiitim .gt. chimax)  chiitim = chimax
      if(chietem .gt. chimax)  chietem = chimax
      if(diffnem .gt. chimax)  diffnem = chimax
!
!
!.....apply relaxation factor to deal with stiff equation
      relfactor = acoef(3010)
      chiitima(j) = (1._R8-relfactor)*chiitima(j) + relfactor*chiitim
      chietema(j) = (1._R8-relfactor)*chietema(j) + relfactor*chietem
      diffnema(j) = (1._R8-relfactor)*diffnema(j) + relfactor*diffnem
      velthi11(j) = (1._R8-relfactor)*velthi11(j) + relfactor*velthi(1,  &  
     & 1)
      velthi21(j) = (1._R8-relfactor)*velthi21(j) + relfactor*velthi(2,  &  
     & 1)
      velthi31(j) = (1._R8-relfactor)*velthi31(j) + relfactor*velthi(3,  &  
     & 1)
      else
        chiitima(j) = 1.E-8_R8
        chietema(j) = 1.E-8_R8
        diffnema(j) = 1.E-8_R8
      endif
!
   19 continue
!
      sqrgps = sqrt(gps)
       dsi0(j) = 1.5_R8*(pmid-emid)*sqrgps*velthi11(j)*udst
       dse0(j) = 1.5_R8*(emid     )*sqrgps*velthi31(j)*udst
       cs0(j)  =                 sqrgps*velthi21(j)*udst
      if(chiitima(j) .lt. 1.E-8_R8) chiitima(j) = 1.E-8_R8
      if(chietema(j) .lt. 1.E-8_R8) chietema(j) = 1.E-8_R8
      if(diffnema(j) .lt. 1.E-8_R8) diffnema(J) = 1.E-8_R8
      if(chiitima(j) .gt. chimax)  chiitima(j) = chimax
      if(chietema(j) .gt. chimax)  chietema(j) = chimax
      if(diffnema(j) .gt. chimax)  diffnema(j) = chimax
      if(itrmod.eq.10) go to 8
!     write(nterm,7771) kcycle, j, chiitim, chietem, diffnem
!    1   ,chiitima(j),chietema(j),diffnema(j)
!7771 format(2i5,1p6e12.4)
      go to 7
!
!
!.......................................................................
!     multimode model MMM_v7.1
!.......................................................................
 15   continue
!.....only call mmm_v7.1 every nskipsf cycles to save time

    IF ( 1 == iskipsf ) THEN

      matdim = 5
      nprout = nout
      lprint = 0
      lsuper = 0
      lreset = 0
      nerr = 0
      npoints = 1
      zrminormm(1) = rminora(j)
      zrmajormm(1) = rmajora(j)
      zelongmm(1)  = elonga(j)
      zdensemm(1) = .5_R8*(ane(j)+ane(j+1))
      zdenshmm(1) = .5_R8*(anhy(j)+anhy(j+1))
      zdensimpmm(1) = .5_R8*(animp(j)+animp(j+1))
      zdensfemm(1) = 0._R8
      zzeffmm(1) = .5_R8*(zeffa(j)+zeffa(j+1))
      ztekevmm(1)= .5_R8*(te(j)+te(j+1))*1.E-3_R8
      ztikevmm(1)= .5_R8*(ti(j)+ti(j+1))*1.E-3_R8
      zqmm(1) = qprof2(j)
      call geval(xsv2(j),2,gvalmm,gpvalmm,gppvalmm,imag,jmag)
      zbtormm(1)      = gvalmm/rmajormm
      zavezimpmm(1)   =  .5_R8*(avezimpa(j)+avezimpa(j+1))
      zamassimpmm(1)  = .5_R8*(amassimpa(j)+amassimpa(j+1))
      zamasshydmm(1)  =  .5_R8*(amasshyda(j)+amasshyda(j+1))
      zaimassmm(1)    = .5_R8*(aimassa(j)+aimassa(j+1))
      zwexbsmm(1)     = 0._R8
      zbefogrd(1)     = -4._R8*rmajormm/(rminora(j+1)-rminora(j-1))
      zgrdnemm(1)     = befogrd*(ane(j+1)-ane(j))/(ane(j)+ane(j+1))
      zgrdnimm(1)  = befogrd*(sumnia(j+1)-sumnia(j))/(sumnia(j)+sumnia(j+1))
      zgrdnhmm(1)     = befogrd*(anhy(j+1)-anhy(j))/(anhy(j)+anhy(j+1))
      zgrdnzmm(1)  = befogrd*(sumnqa(j+1)-sumnqa(j))/(sumnqa(j)+sumnqa(j+1))
      zgrdtemm(1)     = befogrd*(te(j+1)-te(j))/(te(j)+te(j+1))
      zgrdtimm(1)     = befogrd*(ti(j+1)-ti(j))/(ti(j)+ti(j+1))
      zgrdqmm(1)      = -0.25_R8*befogrd*(qprof2(j+1)-qprof2(j-1))/qprof2(j)

      zgvrin(1) = 0.0_R8
      zvtorin(1) = 0.0_R8
      zgvpin(1) = 0.0_R8
      zvpolin(1) = 0.0_R8
      zgvparin(1) = 0.0_R8
      zvparin(1) = 0.0_R8
!
! neoclassical electrical resistivity

      zeta(1) = etpara(j)

! zcsnd0 = speed of sound at magnetic axis

      zcsnd0 = sqrt ( 1.602176487D-19 * te(2) / &
             & ( 1.672621638D-27 * aimassmm ) )

! dkdeps = derivative of elongation wrt r/R

      zdkdeps(1) = ( elonga(j+1) - elonga(j-1) ) / &
             & ( 1.0E-6_R8 + rminora(j+1) / rmajora(j+1) &
             &   - rminora(j-1) / rmajora(j-1) )

! zthtig  = toroidal momentum diffusivity
! zthttig = poloidal momentum diffusivity
!
      npmax = npsit
      if(itrmod.eq.10 .or. itrmod.eq.10) npmax = npsihm

    IF ( j.gt.1 .and. j.le.npmax ) THEN
 
!
      call mmm7_1 (                                                      &  
     &   zrminormm,  zrmajormm,   zelongmm                               &  
     & , zdensemm,   zdenshmm,    zdensimpmm,  zdensfemm                 &  
     & , zzeffmm,    ztekevmm,    ztikevmm,    zqmm,       zbtormm       &  
     & , zavezimpmm, zamassimpmm, zamasshydmm, zaimassmm,  zwexbsmm      &  
     & , zgrdnemm, zgrdnimm, zgrdnhmm, zgrdnzmm, zgrdtemm, zgrdtimm, zgrdqmm &
     & , zgvrin, zvtorin, zgvpin, zvpolin, zgvparin, zvparin &
     & , zeta,     zcsnd0,   zdkdeps  &
     & , zthiig,  zthdig,   ztheig, zthzig,  zthtig, zthttig   &  
     & , velthi = velthi,   vflux = vflux                    &  
     & , npoints = npoints, lprint = lprint,  nprout = nprout &
     & , nerr = nerr )

    IF ( nerr .ne. 0 ) THEN
      ineg=56
      write(nout,6657) j,nerr
      write(nterm,6657) j,nerr
 6657 format(" j, nerr after call to mmm_v7.1", 2i5)
      return
    ENDIF
!
!.....difthi(j1,j2,jz) and velthi(j1,j2,jz) matrices contain all of
!     the anomalous transport
!     (note factor of 3/2 in going from ITER to TSC convention)

      chiitim = 1.5_R8 * zthiig(1)
      chietem = 1.5_R8 * ztheig(1)
      diffnem =          zthdig(1)

IF ( iprint > 0 ) THEN
     write(nterm,'("#2",2i5,1p4e12.4)') &
 kcycle, j, diffnem, thdig ,thdrb ,thdkb
ENDIF

!7771 format(2i5,1p4e12.4)
!
!      write(nterm,7771) kcycle, j, chiitim, chietem, rminormm,
!     1  rmajormm, elongmm, densemm, denshmm, densimpmm, densfemm,
!     2  zeffmm, tekevmm, tikevmm, qmm, btormm,
!     3  avezimpmm, amassimpmm, amasshydmm, aimassmm, wexbsmm,
!     4  grdnemm, grdnimm, grdnhmm, grdnzmm, grdtemm, grdtimm,grdqmm
! 7771 format(2i5, 1p5e12.5,/,10x,5e12.4,/,10x,5e12.4,/,10x,5e12.4,/,
!     5       10x,7e12.4)
!
!.....put limits on chietem and chiitim
      if(chiitim .lt. 1.E-8_R8) chiitim = 1.E-8_R8
      if(chietem .lt. 1.E-8_R8) chietem = 1.E-8_R8
      if(diffnem .lt. 1.E-8_R8) diffnem = 1.E-8_R8
      if(chiitim .gt. chimax)  chiitim = chimax
      if(chietem .gt. chimax)  chietem = chimax
      if(diffnem .gt. chimax)  diffnem = chimax
!
!
!.....apply relaxation factor to deal with stiff equation

      relfactor = acoef(3010)
      chiitima(j) = (1._R8-relfactor)*chiitima(j) + relfactor*chiitim
      chietema(j) = (1._R8-relfactor)*chietema(j) + relfactor*chietem
      diffnema(j) = (1._R8-relfactor)*diffnema(j) + relfactor*diffnem
      velthi11(j) = (1._R8-relfactor)*velthi11(j) + relfactor*velthi(1,  &  
     & 1)
      velthi21(j) = (1._R8-relfactor)*velthi21(j) + relfactor*velthi(2,  &  
     & 1)
      velthi31(j) = (1._R8-relfactor)*velthi31(j) + relfactor*velthi(3,  &  
     & 1)

     ELSE
        chiitima(j) = 1.E-8_R8
        chietema(j) = 1.E-8_R8
        diffnema(j) = 1.E-8_R8
     ENDIF  ! end of IF ( j.gt.1 .and. j.le.npmax ) THEN
!
   ENDIF ! end of IF ( 1 == iskipsf ) THEN
!
      sqrgps = sqrt(gps)
       dsi0(j) = 1.5_R8*(pmid-emid)*sqrgps*velthi11(j)*udst
       dse0(j) = 1.5_R8*(emid     )*sqrgps*velthi31(j)*udst
       cs0(j)  =                 sqrgps*velthi21(j)*udst
      if(chiitima(j) .lt. 1.E-8_R8) chiitima(j) = 1.E-8_R8
      if(chietema(j) .lt. 1.E-8_R8) chietema(j) = 1.E-8_R8
      if(diffnema(j) .lt. 1.E-8_R8) diffnema(J) = 1.E-8_R8
      if(chiitima(j) .gt. chimax)  chiitima(j) = chimax
      if(chietema(j) .gt. chimax)  chietema(j) = chimax
      if(diffnema(j) .gt. chimax)  diffnema(j) = chimax
      if(itrmod.eq.10) go to 8
!     write(nterm,7771) kcycle, j, chiitim, chietem, diffnem
!    1   ,chiitima(j),chietema(j),diffnema(j)
!7771 format(2i5,1p6e12.4)
      go to 7
!
!
 
!
!.......................................................................
!     neoclassical transport model
!.......................................................................
!
 7    continue
!
!
!.....calculate neoclassical transport coefficients
!     REF: Hirshman & Jardin PF 22 (1979) p 741
!
      if(j.gt.1) then
      call geval (xsv2(j),2,gval,gpval,gppval,imag,jmag)
      bsq = bsqar(j)*1.E-8_R8
 
      al0 = rmid *etpera*gxmja2(j)/xmja2(j)/bsq
      alstar = ftrap(j)*gval**2*rmid*etpera/bsq
!     al11bp = alstar*(1.53 - 0.53*ftrap(j))
!     al12bp = alstar*(3.12 - 0.63*ftrap(j))
!     al13   = gval*ftrap(j)*(1.68 - 0.68*ftrap(j))
!     al23   = 1.25*gval*ftrap(j)*(1.-ftrap(j))
!     al11t  = al0*(1. + 2.*qprof2(j))**2
!     al12t  = 1.5*al11t
!     al22nc = 4.66*(al11t + alstar)
!     al33   = -1.26*ftrap(j)*(1. - 0.18*ftrap(j))
      alinc  = 61._R8*(temidd/timidd)**1.5_R8*(al0*(1._R8+ 2._R8*        &  
     & qprof2(j)**2)                                                     &  
     &       + 0.46_R8*alstar/(1._R8-0.54_R8*ftrap(j)) )
      ay = - 1.17_R8*(1._R8-ftrap(j))/(1._R8-0.54_R8*ftrap(j))
      cs2(j) = 0._R8
      cs3(j) = 0._R8
      cs4(j) = 0._R8
!     if(idens.eq.0) then
!     combo1=(tpi*qprof2(j))**2*(al11bp/rmid+al13**2 * etpara(j)/bsq)
!     combo2=(tpi*qprof2(j))**2*(al12bp/rmid+al13*al23*etpara(j)/bsq)
!     cs0(j) = 0.
!     cs1(j) = -combo1*ay*timidd + combo2*temidd
!     cs2(j) =  combo1*(1.+ay)
!     cs3(j) = -combo1*ay       - combo2
!     cs4(j) =  (tpi*qprof2(j))*al13*etpara(j)*gval**2
!    1         /(bsq*tpi**2*vp2(j))
!
!
      coef3nc = -(tpi*qprof2(j))**2*timidd*alinc
!
!     dsi0(j) = 0.
!     dsi1(j) = -timidd*coef3nc
!     dsi2(j) =        coef3nc
!     dsi3(j) =       -coef3nc
!     dsi4(j) = 0.
      chiinc = -coef3nc / (udst*gps)
!
!....NOTE:  Set electron thermal conductivity to ion value.
!     dse1(j) = -coef3nc*temidd
!     dse3(j) =  coef3nc
!
!     befo = (tpi*qprof2(j))**2*temidd
!     combo1 = al12bp + rmid*al23*al13*etpara(j)/bsq
!     combo2 = al22nc + rmid*al23**2 * etpara(j)/bsq
!     dse0(j) = 0
!     dse1(j) = befo*(-ay*timidd*combo1 + temidd*combo2)
!     dse2(j) = befo*(ali12t + combo1*(1.+ay))
!     dse3(j) =-befo*( ay*combo1 + combo2)
!     dse4(j) = (tpi*qprof2(j))*temidd*rmid*al23*etpara(j)*gval**2
!    1         /(bsq*tpi**2*vp2(j))
!     chienc = -coef3nc/(udst*gps)* acoef(73)
      chienc = -coef3nc*acoef(73)/(udst*gps)/sqrt(2._R8*1835._R8)
!
      chiisecj = chiinc
      chiinca(j) = chiinc
      chiesecj = chienc
      chienca(j) = chienc
!
!
      endif
!
!     add contribution from GLF23
!
      if(itrmode.eq.7) go to 4000
      if(itrmod .eq. 14) go to 88
!
!.....enter here for itrmode=8,9,10,11
    8 continue
      chiesecj = sqrt(chietema(j)**2 + chienca(j)**2                     &  
     &                 +chiecopi(j)**2)
      chiisecj = sqrt(chiitima(j)**2 + chiinca(j)**2                     &  
     &                 +chiicopi(j)**2)
      go to 888
88    continue
      tflux = float(j-1)*dpsi
      if(tflux .lt. npsihm*dpsi) then
      chiesecj = sqrt(chietema(j)**2 + chienca(j)**2)
      chiisecj = sqrt(chiitima(j)**2 + chiinca(j)**2)
      else
      chiesecj = chietema(j)
      chiisecj = chiitima(j)
      endif
 888  continue
!     coef3 = - chiisecj*udst*gps
!     coef2 = - chiesecj*udst*gps
!
!     set bounds
!
      if( chiitima(j)+dchii_dtip(j)*dti_dpsi(j) .lt. chiinca(j) ) then
      chiitima(j)=chimin_tmp
      dchii_dtep(j) = 0.0
      dchii_dtip(j) = 0.0
      endif
      if( chietema(j)+dchie_dtep(j)*dte_dpsi(j) .lt. chienca(j) ) then
      chietema(j)=chimin_tmp
      dchie_dtep(j) = 0.0
      dchie_dtip(j) = 0.0
      endif

      coef3 = -(chiisecj+chiitima(j)/chiisecj*dchii_dtip(j)*dti_dpsi(j))
      coef4 = -chiitima(j)/chiisecj*dchii_dtep(j)*dti_dpsi(j)
      coef5 = chiitima(j)/chiisecj*(                                    &
     &                    dchii_dtip(j)*dti_dpsi(j)*dti_dpsi(j)+        &
     &                    dchii_dtep(j)*dti_dpsi(j)*dte_dpsi(j))        
      coef1 = -chietema(j)/chiesecj*dchie_dtip(j)*dte_dpsi(j)
      coef2 = -(chiesecj+chietema(j)/chiesecj*dchie_dtep(j)*dte_dpsi(j))
      coef6 = chietema(j)/chiesecj*(                                    &
     &                    dchie_dtep(j)*dte_dpsi(j)*dte_dpsi(j)+        &
     &                    dchie_dtip(j)*dte_dpsi(j)*dti_dpsi(j))        
 
!
      coef3 = coef3*udst*gps
      coef4 = coef4*udst*gps
      coef5 = coef5*udst*gps*(temidd/temid)
!
      coef1 = coef1*udst*gps
      coef2 = coef2*udst*gps
      coef6 = coef6*udst*gps*(temidd/temid)

      dsi0(j) = coef5*rmid
      dsi1(j) = -(coef3*timidd+coef4*temidd)
      dsi2(j) = coef3
      dsi3(j) = -coef3 + coef4
      dse0(j) = coef6*rmid
      dse1(j) = -coef2*temidd-coef1*timidd
      dse2(j) = coef1
      dse3(j) = coef2 - coef1

!     diffnema(j)=1.0e-08

!     dsi1(j) =  - coef3*timidd
!     dsi2(j) =    coef3
!     dsi3(j) =  - coef3
!
!     dse1(j) =  - coef2*temidd
!     dse3(j) =    coef2
!
      cs1(j) = sqrt(diffnema(j)**2 + diffary(j)**2)*gps*udst/rmid
!
!.....add in convection term only for diffary part  (added 03/28/03)
      dse1(j) = dse1(j) - face*diffary(j)*gps*udst/rmid
      dsi1(j) = dsi1(j) - faci*diffary(j)*gps*udst/rmid
 4000 continue
!
!..Bateman, 11 Dec 2009, Diagnostic printout
!
IF ( iprint > 0 .AND. kcycle > 0 ) THEN

    WRITE(*,*)
    WRITE(*,*) ' In trcdef, kcycle = ', kcycle
    WRITE(*,*)
     WRITE(*,*) '   j   zchii(j,1)    zchii(j,2)   zchii(j,3)' &
      , '   zchii(j,4)   zchii(j,5)'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, zchii(j,1), zchii(j,2), zchii(j,3), zchii(j,4), zchii(j,5)
    ENDDO

    WRITE(*,*)
    WRITE(*,*) ' In trcdef, kcycle = ', kcycle

    WRITE(*,*)
    WRITE(*,*) '   j   rminora' &
      , '  ti_sav(j)   te_sav(j)'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, rminora(j) &
        , ti_sav(j), te_sav(j)
    ENDDO

    WRITE(*,*)
    WRITE(*,*) '   j  rminora(j)' &
      , '  chiitima(j)   chietema(j)   diffnema(j)'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, rminora(j) &
        , chiitima(j), chietema(j), diffnema(j)
    ENDDO

    WRITE(*,*)
    WRITE(*,*) '   j   rminora' &
      , '  dchii_dtip  dchii_dtep  dchie_dtep  dchie_dtip ' &
      , '  dte_dpsi    dti_dpsi'
    DO j=1,npsit
      WRITE(*,"(I5,9ES12.4)") &
        j, rminora(j) &
        , dchii_dtip(j), dchii_dtep(j), dchie_dtep(j), dchie_dtip(j) & 
        , dte_dpsi(j), dti_dpsi(j)
    ENDDO

ENDIF
!
!.....implement sawtooth model to modify:
!                   etpara
!                   dsej , dsij, csj     for j=0,4
!
      call sawtooth
!
!.....define energy confinement times in physical units
      do 4001 j=1,npsit
      gps = gja2(j)/vp2(j)*(tpi*qprof2(j))
      if(gps.eq.0) go to 4001
      chiesec(j) = -(dse3(j)+dse2(j))/(udst*gps)
      chiisec(j) = -dsi2(j)/(udst*gps)
 4001 continue
      chiesec(1)=chiesec(2)
      chiisec(1)=chiisec(2)
      chiecopi(1) = chiecopi(2)
      chiicopi(1) = chiicopi(2)
      chienca(1)  = chienca(2)
      chiinca(1)  = chiinca(2)
!     debug printout
!     do j=1,npsit
!     write(73+i_newton,4073) j,chiesec(j),chienca(j),chietema(j),chiisec(j),chiinca(j),chiitima(j)
!4073 format(i5,1p6e12.4)
!     enddo
!
      if(kcycle.lt.0) return
!
!
!     note:  <J.B> =  (gval**2/(tpi**2*vp2(j)))
!    2                *(adi(j)*gxmja(j)*xmja(j)
!    3                - adi(j-1)*gxmja(j-1)*xmja(j-1))*rdpsi
      do 4003 j=1,npsit
      aimh = .5_R8*(adi(j)+adi(j+1))
!
!.....define loop voltage due to hyperresistivity term
      vlooph(j+1) = -(avhyp(j+1)-avhyp(j))/dpsi
!
!.....subtract off bootstrap and current-drive terms
      as0(j) = -(tpi**2/(aimh*xmja2(j)))*etpara(j)*etafac(j)*            &  
     & tpi*(ajavbs(j)+ajavcd(j)+ajavfw(j)+ajavlh2(j)+ajavec(j))
!
!.....now define cs and ds arrays for convenience - - surface centered
!
      as(j) =                                                            &  
     &    as0(j)                                                         &  
     &  + as1(j)*(adn(j+1)/vp(j+1) - adn(j)/vp(j))*rdpsi                 &  
     &  + as2(j)*(adp(j+1)/vpg(j+1) - adp(j)/vpg(j))*rdpsi               &  
     &  + as3(j)*(ade(j+1)/vpg(j+1) - ade(j)/vpg(j))*rdpsi               &  
     &  + rdpsi*etpara(j)*etafac(j)*tpi**2/(xmja2(j)*aimh)**2            &  
     &  *(adi(j+1)*gxmja(j+1)*xmja(j+1)                                  &  
     &                            -adi(j)*gxmja(j)*xmja(j))
!
!
      cs(j) =                                                            &  
     &    cs0(j)                                                         &  
     &  + cs1(j)*(adn(j+1)/vp(j+1) - adn(j)/vp(j))*rdpsi                 &  
     &  + cs2(j)*(adp(j+1)/vpg(j+1) - adp(j)/vpg(j))*rdpsi               &  
     &  + cs3(j)*(ade(j+1)/vpg(j+1) - ade(j)/vpg(j))*rdpsi               &  
     &  + cs4(j)*(adi(j+1)*gxmja(j+1)*xmja(j+1)                          &  
     &            - adi(j)*gxmja(j)*xmja(j))*rdpsi
!
      dsi(j) =                                                           &  
     &     dsi0(j)                                                       &  
     &   + dsi1(j)*(adn(j+1)/vp(j+1)-adn(j)/vp(j))*rdpsi                 &  
     &   + dsi2(j)*(adp(j+1)/vpg(j+1)-adp(j)/vpg(j))*rdpsi               &  
     &   + dsi3(j)*(ade(j+1)/vpg(j+1)-ade(j)/vpg(j))*rdpsi               &  
     &   + dsi4(j)*(adi(j+1)*gxmja(j+1)*xmja(j+1)                        &  
     &              - adi(j)*gxmja(j)*xmja(j))*rdpsi
!
      dse(j) =                                                           &  
     &     dse0(j)                                                       &  
     &   + dse1(j)*(adn(j+1)/vp(j+1)-adn(j)/vp(j))*rdpsi                 &  
     &   + dse2(j)*(adp(j+1)/vpg(j+1)-adp(j)/vpg(j))*rdpsi               &  
     &   + dse3(j)*(ade(j+1)/vpg(j+1)-ade(j)/vpg(j))*rdpsi               &  
     &   + dse4(j)*(adi(j+1)*gxmja(j+1)*xmja(j+1)                        &  
     &              - adi(j)*gxmja(j)*xmja(j))*rdpsi
!
      emid = .5_R8*(ade(j+1)/vpg(j+1)+ade(j)/vpg(j))
      pmid = .5_R8*(adp(j+1)/vpg(j+1)+adp(j)/vpg(j))
      if(dse(j).ne.0) ratioe(j) = cs(j)*emid/dse(j)
      if(dsi(j).ne.0) ratioi(j) = cs(j)*(pmid-emid)/dsi(j)
 4003 continue
!
!.....define cs and ds at j=1 for each transport model
!
!
      cs(1) = 0._R8
      dsi(1) = 0._R8
      dse(1) = 0._R8
!
!.....adjust value of as(npsit) to be continuous
      as(npsit) = as(npsit-1)
      do 4005 j=npsit+1,npsi
 4005 as(j) = as(npsit)
!
!.....define as array at j=1
      as(1) = (2._R8*as(2) - as(3))
!
!.....second def of energy confinement time
      enerst2 = sum2*udsp
!
      return
!c
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
