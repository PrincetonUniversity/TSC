       subroutine callglf2d(                                             &  
                      !INPUTS
       leigen,                                                           &  
                      ! eigenvalue solver
                      ! 0 for cgg (default), 1 for tomsqz, 2 for zgeev
       nroot,                                                            &  
                      ! no. roots,8 for default, 12 for impurity dynamics
       iglf,                                                             &  
                      ! 0 for original GLF23, 1 for retuned version
       jshoot,                                                           &  
                      ! jshoot=0 time-dep code;jshoot=1 shooting code
       jmm,                                                              &  
                      ! grid number;jmm=0 does full grid jm=1 to jmaxm-1
       jmaxm,                                                            &  
                      ! profile grids 0 to jmaxm
       itport_pt,                                                        &  
                      ! 1:5 transport flags
       irotstab,                                                         &  
                      ! 0 to use egamma_exp; 1 use egamma_m
       te_m,                                                             &  
                      ! 0:jmaxm te Kev           itport_pt(2)=1 transport
       ti_m,                                                             &  
                      ! 0:jmaxm ti Kev           itport_pt(3)=1 transport
       ne_m,                                                             &  
                      ! 0:jmaxm ne 10**19 1/m**3
       ni_m,                                                             &  
                      ! 0:jmaxm ni 10**19 1/m**3 itport_pt(1)=1 transport
       ns_m,                                                             &  
                      ! 0:jmaxm ns 10**19 1/m**3
       i_grad,                                                           &  
                      ! default 0, for D-V method use i_grad=1 to input gradients
       idengrad,                                                         &  
                      ! default 2, for simple dilution
       zpte_in,                                                          &  
                      ! externally provided log gradient te w.r.t rho (i_grad=1)
       zpti_in,                                                          &  
                      ! externally provided log gradient ti w.r.t rho
       zpne_in,                                                          &  
                      ! externally provided log gradient ne w.r.t rho
       zpni_in,                                                          &  
                      ! externally provided log gradient ni w.r.t rho
       angrotp_exp,                                                      &  
                      ! 0:jmaxm exp plasma toroidal angular velocity 1/sec
                      ! if itport_pt(4)=0 itport_pt(5)=0
       egamma_exp,                                                       &  
                      ! 0:jmaxm exp exb shear rate in units of csda_exp
                      ! if itport_pt(4)=-1 itport_pt(5)=0
       gamma_p_exp,                                                      &  
                      ! 0:jmaxm exp par. vel. shear rate in units of csda_exp
                      ! if itport_pt(4)=-1 itport_pt(5)=0
       vphi_m,                                                           &  
                      ! 0:jmaxm toroidal velocity m/sec
                      ! if itport_pt(4)=1 itport_pt(5)=0 otherwise output
       vpar_m,                                                           &  
                      ! 0:jmaxm parallel velocity m/sec
                      ! if itport_pt(4)=1 itport_pt(5)=1 otherwise output
       vper_m,                                                           &  
                      ! 0:jmaxm perp. velocity m/sec
                      ! if itport_pt(4)=1 itport_pt(5)=1 otherwise output
       zeff_exp,                                                         &  
                      ! 0:jmaxm ne in 10**19 1/m**3
       bt_exp,                                                           &  
                      ! vaccuum axis toroidal field in tesla
       bt_flag,                                                          &  
                      ! switch for effective toroidal field use in rhosda
       rho,                                                              &  
                      ! 0:jmaxm 0 < rho < 1 normalized toroidal flux (rho=rho/rho(a))
       arho_exp,                                                         &  
                      ! rho(a), toroidal flux at last closed flux surface (LCFS)
                      !   toroidal flux= B0*rho_phys**2/2 (m)
                      !   B0=bt_exp, arho_exp=rho_phys_LCFS
       gradrho_exp,                                                      &  
                      ! 0:jmaxm dimensionless <|grad rho_phys |**2>
       gradrhosq_exp,                                                    &  
                      ! 0:jmaxm dimensionless <|grad rho_phys |>
                      !NOTE:can set arho_exp=1.,if gradrho_exp=<|grad rho |>
                      !                 and gradrhosq_exp = <|grad rho |**2>
       rmin_exp,                                                         &  
                      ! 0:jmaxm minor radius in meters
       rmaj_exp,                                                         &  
                      ! 0:jmaxm major radius in meters
       rmajor_exp,                                                       &  
                      ! axis major radius
       zimp_exp,                                                         &  
                      ! effective Z of impurity
       amassimp_exp,                                                     &  
                      ! effective A of impurity
       q_exp,                                                            &  
                      ! 0:jmaxm safety factor
       shat_exp,                                                         &  
                      ! 0:jmaxm magnetic shear, d (ln q_exp)/ d (ln rho)
       alpha_exp,                                                        &  
                      ! 0:jmaxm MHD alpha from experiment
       elong_exp,                                                        &  
                      ! 0:jmaxm elongation
       amassgas_exp,                                                     &  
                      !  atomic number working hydrogen gas
       alpha_e,                                                          &  
                      ! 1 full (0 no) no ExB shear stab
       x_alpha,                                                          &  
                      ! 1 full (0 no) alpha stabilization  with alpha_exp
                      !-1 full (0 no) self consistent alpha_m stab.
       i_delay,                                                          &  
                      !i_delay time delay for ExB shear should be non-zero only
                      ! once per step and is less than or equal 10
                      !OUTPUTS
       diffnem,                                                          &  
                      ! ion plasma diffusivity in m**2/sec
       chietem,                                                          &  
                      ! electron ENERGY diffuivity in m**2/sec
       chiitim,                                                          &  
                      ! ion      ENERGY diffuivity in m**2/sec
       etaphim,                                                          &  
                      ! toroidal velocity diffusivity in m**2/sec
       etaparm,                                                          &  
                      ! parallel velocity diffusivity in m**2/sec
       etaperm,                                                          &  
                      ! perpendicular velocity diffusivity in m**2/sec
       exchm,                                                            &  
                      ! turbulent electron to ion ENERGY exchange in MW/m**3
                      ! 0:jmaxm values
        diff_m,                                                          &  
        chie_m,                                                          &  
        chii_m,                                                          &  
        etaphi_m,                                                        &  
        etapar_m,                                                        &  
        etaper_m,                                                        &  
        exch_m,                                                          &  
        egamma_m,                                                        &  
                      !0:jmaxm exb shear rate in units of local csda_m
        egamma_d,                                                        &  
                      !0:jmaxm exb shear rate delayed by i_delay steps
        gamma_p_m,                                                       &  
                      !0:jmaxm par. vel. shear rate in units of local csda_m
        anrate_m,                                                        &  
                      !0:jmaxm leading mode rate in unints of local csda_m
        anrate2_m,                                                       &  
                      !0:jmaxm 2nd mode rate in units of local csda_m
        anfreq_m,                                                        &  
                      !0:jmaxm leading mode frequency
        anfreq2_m                                                        &  
                      !0:jmaxm 2nd mode frequency
       )
!@callglf2d.f
! 12-feb-03 Kinsey, v1.61 retuned GLF23 model
! 11-apr-01 Kinsey, updated for v1.50
! 29-aug-00 Kinsey, added ave_ve for 3pt smoothing of ve (0 for none)
! 05-nov-99 Kinsey, precall routine for glf2d.f
! added i_dengrad switch for dilution
!***********************************************************************
!************************************************************************
 
! see glf2d documentation for use of diffusivities in transport equations
! ITER definitions of diffusivity used.
! must add neoclassical diffusion
 
!.......begin common block ....................
!      only common used for communication with glf2d....designation by xxxxx_gf
 
 
      USE GLF
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
!.......end common block....................
 
!.......begin dimensions.................
 
      REAL*8 epsilon, zeps, zpi
      parameter ( epsilon = 1.E-34_R8, zeps = 1.E-6_R8)
 
! external arrays
 
      integer jmaxm, jshoot, jmm, i_grad, idengrad, itport_pt(1:5),      &  
     &  i_delay, j, jin, jout, jm, irotstab, iglf,                       &  
     &  jpt, jptf, jptt, jj, ii, ik, bt_flag, leigen, nroot
      REAL*8 alpha_e, x_alpha, xalpha,                                   &  
     &  diffnem, chietem, chiitim, etaphim,                              &  
     &  etaparm, etaperm, exchm,                                         &  
     &  rmajor_exp, zimp_exp, amassimp_exp,                              &  
     &  bt_exp, arho_exp, amassgas_exp,                                  &  
     &  cbetae, cxnu, relx, cmodel, drho, zeff_e,                        &  
     &  zpmte, zpmti, zpmne, zpmni, vstar_sign,                          &  
     &  egeo_local, pgeo_local, rdrho_local, rdrho_local_p1, fc,         &  
     &  akappa1, alpha_neo, alpha_neo_hold,                              &  
     &  zpmne_q, zpmni_q, zpmnimp, gfac
      REAL*8 te_m(0:jmaxm),ti_m(0:jmaxm),                                &  
     &  ne_m(0:jmaxm),ni_m(0:jmaxm), ns_m(0:jmaxm),                      &  
     &  vphi_m(0:jmaxm),angrotp_exp(0:jmaxm),                            &  
     &  egamma_exp(0:jmaxm),gamma_p_exp(0:jmaxm),                        &  
     &  vpar_m(0:jmaxm),vper_m(0:jmaxm),                                 &  
     &  rho(0:jmaxm),rmin_exp(0:jmaxm),rmaj_exp(0:jmaxm),                &  
     &  gradrho_exp(0:jmaxm),gradrhosq_exp(0:jmaxm),                     &  
     &  zeff_exp(0:jmaxm),q_exp(0:jmaxm),shat_exp(0:jmaxm),              &  
     &  bteff_exp(0:jmaxm)
      REAL*8 alpha_exp(0:jmaxm),elong_exp(0:jmaxm),                      &  
     &  diff_m(0:jmaxm),chie_m(0:jmaxm),chii_m(0:jmaxm),                 &  
     &  etaphi_m(0:jmaxm),                                               &  
     &  etapar_m(0:jmaxm),etaper_m(0:jmaxm), exch_m(0:jmaxm),            &  
     &  egamma_m(0:jmaxm),egamma_d(0:jmaxm,10),gamma_p_m(0:jmaxm),       &  
     &  anrate_m(0:jmaxm), anrate2_m(0:jmaxm),                           &  
     &  anfreq_m(0:jmaxm), anfreq2_m(0:jmaxm)
 
! internal arrays (which can be converted to externals)
 
      REAL*8 zpte_in, zpti_in, zpne_in, zpni_in,                         &  
     &  zpte_m(0:jmaxm),zpti_m(0:jmaxm),                                 &  
     &  zpne_m(0:jmaxm),zpni_m(0:jmaxm),                                 &  
     &  drhodr(0:jmaxm),drhodrrrho(0:jmaxm),geofac(0:jmaxm),             &  
     &  rhosda_m(0:jmaxm),csda_m(0:jmaxm),cgyrobohm_m(0:jmaxm),          &  
     &  betae_m(0:jmaxm),xnu_m(0:jmaxm),                                 &  
     &  alpha_m(0:jmaxm),vstarp_m(0:jmaxm)
 
! working arrays and variables
 
      REAL*8 ve(0:jmaxm),vpar(0:jmaxm)
!     real*8 vmode(0:jmaxm)
!     real*8 kevdsecpmw
 
! diagnostic arrays (if used)
 
!     real*8 vstar_m(0:jmaxm),vexb_m(0:jmaxm)
!     real*8 vmode_m(0:jmaxm)
!     real*8 gamma_mode_m(0:jmaxm)
!     real*8 gamma_k_j(20,0:jmaxm),freq_k_j(20,0:jmaxm)
!     real*8 chie_k_j(20,0:jmaxm),chii_k_j(20,0:jmaxm)
!     real*8 vnewstare_m(0:jmaxm),vnewstari_m(0:jmaxm)
!     real*8 ky_j(0:jmaxm)
!     real*8 gamma_j(0:jmaxm,1:4),freq_j(0:jmaxm,1:4)
!     real*8 phi_norm_j(0:jmaxm,1:4)
!     real*8 dnrate_m(0:jmaxm), dtnrate_m(0:jmaxm)
!     real*8 dnfreq_m(0:jmaxm)
 
! some internals
 
      REAL*8 tem,tim,nem,nim,nsm,zeffm, aiwt_jp1,                        &  
     &       xnimp_jp1, xnimp, vnewk3x
 
!.......end   dimensions.................
 
!    jm is local grid 0 < jin_m < j < jout_m < jmaxm
!    jm must be greater than 0 and less than jmaxm
!
!
!...constants
      zpi = atan2 ( 0.0_R8, -1.0_R8)
!
!dmc      write(6,7701) zpi
!dmc 7701 format(' zpi = ',1pe17.10)
!dmc      write(6,7702) jshoot,jmm,jmaxm,(itport_pt(j),j=1,5)
!dmc 7702 format(/' jshoot,jmm,jmaxm = ',3(1x,i6)/' itport_pt = ',5(1x,i3))
!dmc      call echo('te_m',te_m,jmaxm)
!dmc      call echo('ti_m',ti_m,jmaxm)
!dmc      call echo('ne_m',ne_m,jmaxm)
!dmc      call echo('ni_m',ni_m,jmaxm)
!
!dmc      write(6,7703) i_grad,zpte_in,zpti_in,zpne_in,zpni_in
!dmc 7703 format(/' igz: ',i1,1x,4(1pe17.10))
!
!dmc      call echo('angrotp_exp',angrotp_exp,jmaxm)
!dmc      call echo('egamma_exp',egamma_exp,jmaxm)
!dmc      call echo('gamma_p_exp',gamma_p_exp,jmaxm)
!dmc      call echo('vphi_m',vphi_m,jmaxm)
!dmc      call echo('vpar_m',vpar_m,jmaxm)
!dmc      call echo('vper_m',vper_m,jmaxm)
!dmc      call echo('zeff_exp',zeff_exp,jmaxm)
!
!dmc      write(6,7705) bt_exp
!dmc 7705 format(/' bt_exp = ',1pe17.10)
!
!dmc      call echo('rho',rho,jmaxm)
!
!dmc      write(6,7706) arho_exp
!dmc 7706 format(/' arho_exp = ',1pe17.10)
!
!dmc      call echo('gradrho_exp',gradrho_exp,jmaxm)
!dmc      call echo('gradrhosq_exp',gradrhosq_exp,jmaxm)
!dmc      call echo('rmin_exp',rmin_exp,jmaxm)
!dmc      call echo('rmaj_exp',rmaj_exp,jmaxm)
!
!dmc      write(6,7707) rmajor_exp
!dmc 7707 format(/' rmajor_exp = ',1pe17.10)
!
!dmc      call echo('q_exp',q_exp,jmaxm)
!dmc      call echo('shat_exp',shat_exp,jmaxm)
!dmc      call echo('alpha_exp',alpha_exp,jmaxm)
!dmc      call echo('elong_exp',elong_exp,jmaxm)
!
!dmc      write(6,7708) amassgas_exp
!dmc 7708 format(/' atomic no.:  ',1pe17.10)
!
!dmc      write(6,7709) alpha_e,x_alpha
!dmc 7709 format(/' alpha_e, x_alpha = ',2(1x,1pe17.10))
!
!dmc      write(6,7710) i_delay
!dmc 7710 format(/' i_delay = ',i5)
!
!...initialize variables
!
      do j=0,jmaxm
        zpte_m(j) = 0._R8
        zpti_m(j) = 0._R8
        zpne_m(j) = 0._R8
        zpni_m(j) = 0._R8
 
        betae_m(j)     = 0._R8
        xnu_m(j)       = 0._R8
        cgyrobohm_m(j) = 0._R8
        rhosda_m(j)    = 0._R8
        csda_m(j)      = 0._R8
 
        geofac(j)      = 0._R8
        drhodr(j)      = 0._R8
        drhodrrrho(j)  = 0._R8
 
        gamma_p_m(j)   = 0._R8
        egamma_m(j)    = 0._R8
        ve(j)          = 0._R8
        vper_m(j)      = 0._R8
        vpar(j)        = 0._R8
!        vphi_m(j)      = 0.D0
        vstarp_m(j)    = 0._R8
        alpha_m(j)     = 0._R8
 
        anrate_m(j)    = 0._R8
        anrate2_m(j)   = 0._R8
        anfreq_m(j)    = 0._R8
        anfreq2_m(j)   = 0._R8
 
        exch_m(j)      = 0._R8
        diff_m(j)      = 0._R8
        chie_m(j)      = 0._R8
        chii_m(j)      = 0._R8
        etaphi_m(j)    = 0._R8
        etapar_m(j)    = 0._R8
        etaper_m(j)    = 0._R8
      enddo
!
! diagnostic arrays (not used)
!
!     do j=0,jmaxm
!       vstar_m(j)     = 0.
!       vexb_m(j)      = 0.
!       dnrate_m(j)    = 0.
!       dtnrate_m(j)   = 0.
!       dnfreq_m(j)    = 0.
!      do k=1,4
!       gamma_j(j,k)    = 0.
!       freq_j(j,k)     = 0.
!       phi_norm_j(j,k) = 0.
!      enddo
!     enddo
!
!     do j=0,jmaxm
!      do ik=1,20
!       gamma_k_j(ik,j) = 0.D0
!       freq_k_j(ik,j)  = 0.D0
!       chie_k_j(ik,j)  = 0.D0
!       chii_k_j(ik,j)  = 0.D0
!      enddo
!      vnewstare_m(j) = 0.
!      vnewstari_m(j) = 0.
!     enddo
!
!
!***********************************************************************
!mnt   profiles of quantities derived from model profiles of te,ti,ne
!mnt   revised derivedmodlocal to center between jm-jptf and jm+ptf
!***********************************************************************
 
!.......begin switches and settings......
 
!**********************************************************************
!     glf23 parameters
!xx      settings same as function glf23_v4_1_10
 
!k      write(6,*)  'jmm=',jmm,'going in'
 
      eigen_gf = leigen
!      nroot_gf=8     ! 8 for pure plasma, 12 for full impurity dynamics
      nroot_gf=nroot
      iflagin_gf(1)=0
      iflagin_gf(2)=1
      iflagin_gf(3)=1
      iflagin_gf(4)=0
      iflagin_gf(5)=3
 
      xparam_gf(1)=0._R8
      xparam_gf(2)=0
      xparam_gf(3)=.7_R8
      xparam_gf(4)=0._R8
      xparam_gf(6)=0._R8
      xparam_gf(7)=1._R8
      xparam_gf(8)=0._R8
      xparam_gf(9)=1._R8
      xparam_gf(10)=0._R8
      xparam_gf(11)=0._R8
      xparam_gf(12)=0._R8
      xparam_gf(13)=0.2_R8
      xparam_gf(14)=1._R8
      xparam_gf(15)=-0.1_R8
      xparam_gf(16)=0._R8
      xparam_gf(17)=0.1_R8
      xparam_gf(18)=.0_R8
      xparam_gf(19)=0._R8
      xparam_gf(20)=0._R8
      xparam_gf(21)=0._R8
      xparam_gf(22)=0._R8
      xparam_gf(23)=1._R8
      xparam_gf(24)=0._R8
      xparam_gf(25)=0._R8
      xparam_gf(26)=0._R8
      xparam_gf(27)=0._R8
      xparam_gf(28)=0._R8
      xparam_gf(29)=0._R8
      xparam_gf(30)=0._R8
!
      xky0_gf= .2_R8
      rms_theta_gf=zpi/3._R8
      park_gf  =0.7_R8
      ghat_gf  =1._R8
      gchat_gf =1._R8
 
      adamp_gf=.50_R8
      alpha_star_gf  =0._R8
      alpha_mode_gf=0._R8
      gamma_e_gf  =0._R8
!temp
      gamma_e_gf  =-.000000000001_R8
      xkdamp_gf     =0._R8
 
      alpha_p_gf=0.50_R8
 
!   cbetae=1 is full electromagetic
      cbetae=1.E-6_R8
!      cbetae=1.D0
!   full collisionality
      cxnu=1._R8
 
      cnorm_gf=100._R8
 
      ikymax_gf=10
      xkymin_gf=.02_R8
      xkymax_gf=.5_R8
 
!# non glf23 paramerter
 
       cmodel=1._R8
       xalpha=x_alpha
 
!      ialphastab=1
!      ineo=-2
 
!      iexch_m=1
!      iexp_exch=-1
 
!      i_dengrad=2
!      iexp_imp=1
!      igeo_m=3
 
!      irotstab=1
!      irot1=1
!      irot2=1
!
!xx      endf
 
!......begin important optional settings
 
! turn on self-consistant alpha-stabilization
!       ialphastab=1
 
!       turn on EXB shear stabilization
!       alpha_e_gf=1. full on ExB shear
        alpha_e_gf=alpha_e
 
!       turn on self consistant EXB stabilization
!       irotstab=1
 
! itport_pt(1)=1 plasma transport on; itport_pt(2:3)=1; electron and ion heat on
! itport_pt(4)=1 ;itport_pt(5)=0  transport vphi with neoclassical determining vtheta
! itport_pt(4)=1 ;itport_pt(5)=1  transport vphi and vtheta with fast time scale
!       neoclassical drag built into vphi and vtheta transport equations...
!       consult G.M. Staebler
 
! if only vphi_exp is available itport_pt(4)=0
 
!      grid-centering in computing EXB shear  span jm-jptf to jm+jpt
 
!      turn on high-k eta-e modes
        xparam_gf(10)=1._R8
!
!    relaxation turned off relx=0.
!    relaxation can be turned on for one call per step
       relx=0._R8
!
! settings for retuned GLF23 model
!
       if (iglf.eq.1) then      ! retuned model
         cnorm_gf=50._R8  ! ITG normalization (via GYRO runs)
         xparam_gf(10)=12._R8  ! ETG normalization (cnorm*xparam(10))
         iflagin_gf(5)=5        ! rms theta fit formula
         xparam_gf(13)=0.15_R8  ! rms_theta q-dependence
         xparam_gf(16)=0.15_R8  ! rms_theta shat dependence
         xparam_gf(17)=0.25_R8  ! rms_theta shat dependence
         xparam_gf(19)=1.0_R8  ! rms_theta alpha dependence
         adamp_gf=.70_R8  ! radial mode damping exponent
         alpha_p_gf=0.35_R8  ! parallel velocity shear fit
         park_gf=0.8_R8  ! parallel ion motion fit
         bt_flag=1              ! use real geometry ExB shear
       endif
!
!.......end important optional settings
 
!.......end switches and settings......
 
!.......start setups.....................
!
 
!*********************************************************************************
! GEOMETRY FACTORS NEEDED FOR SETUP
! external
!      rho(jm)    :toroidal flux co-ordinate 0:50 grids 0->1
!      gradrhosq_exp(jm) : <|grad rho_phys |**2> toroidal flux= B0*rho_phys**2/2
!                 rho_phys=rho*arho_exp
!                 hence   gradrhosq_exp ->1 for a circle
!      gradrho_exp(jm)  : <|grad rho_phys |>
! internal
!      drhodr(jm)
!      drhodrrrho(jm)
!      geofac(jm)
!
!        geofac(j)=gradrho_exp(j)*(rho(j+1)-rho(j))*arho_exp
!     >   /(rmin_exp(j+1)-rmin_exp(j))/gradrhosq_exp(j)
!
!        drhodr(j)=(rho(j+1)-rho(j))*arho_exp/
!     >   (rmin_exp(j+1)-rmin_exp(j))
!
!        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/arho_exp/rho(j)
 
!  surface factor for power flow
!        sfactor(j)=2.*pi_m*arho_exp*rho(j)*h_exp(j)*2.*pi_m*rmajor_exp
!        h_exp(j-1)=hcap_d(j) hcap in ONETWO
 
!******************************************************************************
 
      if(jmm.gt.0) then
       jin=jmm
       jout=jmm
      endif
      if(jmm.eq.0) then
       jin=1
       jout=jmaxm-1
      endif
      do jm=jin,jout
 
! time dependent codes jshoot=0
! diffusion coefficients and gradients between jm+1 and jm
 
! outside to inside shooting codes jshoot=1 diffusion coefficient at jm
! gradient is forward jm to jm-1
! backward gradient is jm to jm+1
! shear is between forward and backward gradient.
! forward is implicit and backward is already updated
 
       tem=(te_m(jm+1-jshoot)+te_m(jm))/2._R8
       tim=(ti_m(jm+1-jshoot)+ti_m(jm))/2._R8
       nem=(ne_m(jm+1-jshoot)+ne_m(jm))/2._R8
       nim=(ni_m(jm+1-jshoot)+ni_m(jm))/2._R8
       nsm=(ns_m(jm+1-jshoot)+ns_m(jm))/2._R8
       zeffm=(zeff_exp(jm+1-jshoot)+zeff_exp(jm))/2._R8
 
       betae_m(jm) = 400._R8*nem*tem/(1.E5_R8*bt_exp**2)
!      betai_m(jm) = 400.*nim*tim/(1.e5*bt_exp**2)
 
 
!rew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
       vnewk3x=                                                          &  
     &   0.117_R8*nem*tem**(-1.5_R8)/(tim**0.5_R8)*arho_exp*             &  
     &   (amassgas_exp/2._R8)**0.5_R8
!rew   as used in gks multiply by 1/2 and takout any 1/2 factor in solfp
!rew          vnewk3x=vnewk3x/2.
       xnu_m(jm) =vnewk3x/(2._R8*tem/tim)**0.5_R8
!rew  10/25/95 fixed zeff+1 factor: zeff col with ions;1 col with elecs.
       zeff_e=0._R8
       xnu_m(jm) = xnu_m(jm)*(zeff_exp(jm)+zeff_e)
 
 
!      vnewstare_m(jm)=zeff_exp(jm) *2.91e-6*nem*1.e13*15./
!    >  (tem*1.e3)**2*rmaj_exp(jm)*100.*q_exp(jm)
!    >     /(rmin_exp(jm)/rmaj_exp(jm)+epsilon)**1.5/4.19e7
!
!      vnewstari_m(jm)=4.78e-8*nem*1.e13*15./
!    >  (tim*1.e3)**2*rmaj_exp(jm)*100.*q_exp(jm)
!    >     /(rmin_exp(jm)/rmaj_exp(jm)+epsilon)**1.5/9.79e5
 
!
      if(jm.eq.0.or.jm.eq.jmaxm) then
       write(6,*) 'can not call callglf2d for this jm'
      endif
 
      jptf=0
       if(jm.eq.1) jptf=0
      jpt=1
      jptt=jpt
       if(jptt.lt.1) jptt=1
       if(jm.eq.jmaxm-1) jptt=1
      jpt=jptt
 
!  note: dependence of shear's  on zpxx_m(jm-1) and zpxx_m(jm+1)
!  and alpha(jm) depends on zpxx_m(jm+1)
 200  format(i2,2x,0p1f4.2,0p5f10.5)
!
      do j=jm-jptf,jm+jpt
!
!... some geometric factors
!
        geofac(j)=gradrho_exp(j)*(rho(j)-rho(j-1))*arho_exp              &  
     &   /(rmin_exp(j)-rmin_exp(j-1)+epsilon)/gradrhosq_exp(j)
 
        drhodr(j)=(rho(j)-rho(j-1))*arho_exp/                            &  
     &   (rmin_exp(j)-rmin_exp(j-1)+epsilon)
 
        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/                             &  
     &   arho_exp/(rho(j)+epsilon)
!
!... local rate unit
!
       csda_m(j)=9.79E5_R8*(te_m(j)*1.E3_R8)**.5_R8/                     &  
     &    (arho_exp*100._R8)/amassgas_exp**0.5_R8
!
!... local rho_star
! Note: use effective B-field if bt_flag > 0
!
       if (bt_flag .gt. 0) then
         bteff_exp(j)=bt_exp*rho(j)*arho_exp/                            &  
     &         rmin_exp(j)*drhodr(j)
         rhosda_m(j)=((1.02E2_R8*(te_m(j)*1.E3_R8)**.5_R8)/bteff_exp(j)  &  
     &         /1.E4_R8)*amassgas_exp**.5_R8/(arho_exp*100._R8)
       else
         rhosda_m(j)=((1.02E2_R8*(te_m(j)*1.E3_R8)**.5_R8)/bt_exp/       &  
     & 1.E4_R8)                                                          &  
     &         *amassgas_exp**.5_R8/(arho_exp*100._R8)
       endif
      enddo
!
!   local gyrobohm unit of diffusion
!
       cgyrobohm_m(jm)=1.E-4_R8*                                         &  
     &  9.79E5_R8*(tem*1.E3_R8)**.5_R8/(arho_exp*100._R8)                &  
     &  *(1.02E2_R8*(tem*1.E3_R8)**.5_R8/bt_exp/1.E4_R8)**2*             &  
     & amassgas_exp**.5_R8
!
      do j=jm-jptf, jm+jpt
        drho=rho(j-1)-rho(j)+epsilon
        zpte_m(j)=-(log(te_m(j-1))-log(te_m(j)))/drho
        zpti_m(j)=-(log(ti_m(j-1))-log(ti_m(j)))/drho
        zpne_m(j)=-(log(ne_m(j-1))-log(ne_m(j)))/drho
        zpni_m(j)=-(log(ni_m(j-1))-log(ni_m(j)))/drho
      enddo
 
        zpmte=zpte_m(jm+1-jshoot)
        zpmti=zpti_m(jm+1-jshoot)
        zpmne=zpne_m(jm+1-jshoot)
        zpmni=zpni_m(jm+1-jshoot)
!
!... check on zero norm gradients:
!
        if (abs(zpmti).lt.zeps) zpmti=zeps
        if (abs(zpmte).lt.zeps) zpmte=zeps
        if (abs(zpmne).lt.zeps) zpmne=zeps
        if (abs(zpmni).lt.zeps) zpmni=zeps
!
        if(i_grad.eq.1) then
         zpmte=zpte_in
         zpmti=zpti_in
         zpmne=zpne_in
         zpmni=zpni_in
        endif
 
! MHD alpha parameter
 
         alpha_m(jm)=drhodr(jm)*                                         &  
     &    q_exp(jm)**2*rmaj_exp(jm)/arho_exp*                            &  
     &    betae_m(jm)*((tim/tem*nim/nem)*                                &  
     &   (zpmni+zpmti)                                                   &  
     &    +zpmne+zpmte)
 
! vstarp_m is diamagnetic part of egamma_m (doppler shear rate)
! vstar_sign is negative for negative vstar_i. Thus for co-injection or positive
! angrot toroidal rotation cancels the diamgnetic rotation
 
       vstar_sign=-1._R8
 
        j=jm
          rho(j-jptf)=rho(j-jptf)+epsilon
          rho(j+jpt)=rho(j+jpt)+epsilon
          egeo_local=1._R8
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
          rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
 
        vstarp_m(jm)=(                                                   &  
     & +egeo_local*vstar_sign*                                           &  
     & (rho(j+jpt)*rdrho_local_p1+rho(j-jptf)*rdrho_local)/2._R8*(       &  
     &(ti_m(j+jpt)/te_m(j+jpt))*csda_m(j+jpt)*                           &  
     & (zpti_m(j+jpt)+zpni_m(j+jpt))                                     &  
     &  *pgeo_local*rhosda_m(j+jpt)/rho(j+jpt)/rdrho_local_p1            &  
     &-(ti_m(j-jptf)/te_m(j-jptf))*csda_m(j-jptf)*                       &  
     & (zpti_m(j-jptf)+zpni_m(j-jptf))                                   &  
     &  *pgeo_local*rhosda_m(j-jptf)/rho(j-jptf)/rdrho_local             &  
     &  )/(rho(j+jpt)-rho(j-jptf)+epsilon)/csda_m(j)                     &  
     &  )
!
       do jj=1,2
 
        if(jj.eq.1) j=jm-jptf
 
        if(jj.eq.2) j=jm+jpt
 
! banana regime ie collisionless limit formulas
 
        fc=1-1.46_R8*(rmin_exp(j)/rmaj_exp(j))**0.5_R8+                  &  
     &      0.46_R8*(rmin_exp(j)/rmaj_exp(j))**1.5_R8
        akappa1=0.8839_R8*fc/(0.3477_R8+0.4058_R8*fc)
        alpha_neo=-akappa1+1._R8
        alpha_neo_hold=alpha_neo
 
!xc angrotp is plasma rotation
!xc angrot is impurity rotation from experiment
!xc if angrotp_exp is not suppied must insert
!x
!x        akappa2=(1.-fc)/(1.+1.1671*fc)
!xc trace impurity limit to go from impurity rotation angrot_exp to
!xc plasma rotation angrotp_exp
!x        angrotp_exp(j)=angrot_exp(j)+
!x     >     akappa2*3./2.*csda_exp(j)*zpti_exp(j)*(ti_exp(j)/te_exp(j))
!x     >     /rho(j)*q_exp(j)*rhosda_exp(j)*pgeo_local/rdrho_local
!x        if(angrot_exp(j).eq.0.) angrotp_exp(j)=0.
!x        angrotp_exp(j)=corot*angrotp_exp(j)
 
          egeo_local=1._R8
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j)/arho_exp/(rho(j)+epsilon)
 
        ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*         &  
     & (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local             &  
     & -rho(j)*rdrho_local*                                              &  
     &  arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)
 
        vpar(j)=rmajor_exp*angrotp_exp(j)-vstar_sign*                    &  
     & (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*      &  
     &((alpha_neo-1._R8)*zpti_m(j))*rho(j)*rdrho_local*                  &  
     & arho_exp/rmajor_exp/q_exp(j)
 
!        vmode(j)=anfreq_m(j)/(ky_j(j)+epsilon)*
!     >    csda_m(j)*arho_exp*rhosda_m(j)
 
       if(itport_pt(4).eq.0.and.itport_pt(5).eq.0) then
         vphi_m(j)=rmajor_exp*angrotp_exp(j)
 
         vpar_m(j)=vpar(j)
 
         vper_m(j)=ve(j)                                                 &  
     &  +(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*               &  
     &   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
       if(abs(itport_pt(4)).eq.1.and.itport_pt(5).eq.1) then
         vpar(j)=vpar_m(j)
       endif
 
       if(abs(itport_pt(4)).eq.1.and.itport_pt(5).eq.0) then
! this option vpar is vphi vexb from neo+vphi
         ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*        &  
     &  (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local            &  
     &  -rho(j)*rdrho_local*                                             &  
     &   arho_exp/rmajor_exp/q_exp(j)*vphi_m(j)
 
         vpar(j)=vphi_m(j)-vstar_sign*                                   &  
     &   (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*    &  
     &   ((alpha_neo-1._R8)*zpti_m(j))*rho(j)*rdrho_local*               &  
     &   arho_exp/rmajor_exp/q_exp(j)
 
         vpar_m(j)=vpar(j)
 
         vper_m(j)=ve(j)                                                 &  
     &   +(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*              &  
     &   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
       if(itport_pt(5).eq.1) then
! this option vexb from vper and vpar with neo dampng built into
! vpar and vper transport equations
         ve(j)=vper_m(j)                                                 &  
     &   -(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*              &  
     &   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
!       vexb_m(j)=ve(j)
!       vmode_m(j)=vmode(j)
!       vstar_m(j)=(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
!    > (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
 
      enddo
!
! compute shears from outside minus inside
!
        j=jm
        rho(j-jptf)=rho(j-jptf)+epsilon
        rho(j+jpt)=rho(j+jpt)+epsilon
        egeo_local=1._R8
        pgeo_local=drhodr(j)
        rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
        rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
!
        egamma_m(jm)=relx*egamma_m(jm)+(1._R8-relx)*(                    &  
     &  egeo_local*drhodrrrho(j)*                                        &  
     &  (rho(j-jptf)+rho(j+jpt))/(q_exp(j-jptf)+q_exp(j+jpt))*           &  
     &  (ve(j+jpt)*q_exp(j+jpt)/rho(j+jpt)/rdrho_local_p1-               &  
     &  ve(j-jptf)*q_exp(j-jptf)/rho(j-jptf)/rdrho_local)/               &  
     &  (rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)              &  
     &  )
!
!       write(*,*) jm, rho(jm), egamma_m(jm), ' egamma'
!
!       gamma_mode_m(jm)=relx*gamma_mode_m(jm)+(1.-relx)*(
!    >  egeo_local*drhodrrrho(j)*
!    >  (rho(j-jptf)+rho(j+jpt))/2.*
!    >  (vmode(j+jpt)/rho(j+jpt)/rdrho_local_p1-
!    >  vmode(j-jptf)/rho(j-jptf)/rdrho_local)/
!    >  (rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
!    >  )
 
        gamma_p_m(jm)=relx*gamma_p_m(jm)+(1._R8-relx)*(                  &  
     &   -drhodr(j)*                                                     &  
     &   (vpar(j+jpt)-vpar(j-jptf))                                      &  
     &   /(rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)            &  
     &  )
 
       if (jm.eq.1.and.gamma_p_m(jm).gt.10)                              &  
     &    gamma_p_m(jm)=10._R8
        alpha_neo=alpha_neo_hold
 
!.......end   setups...........
 
!***********************************************************************
!vv      subroutine model
!***********************************************************************
 
!   units:
!     diffusion (m**2/sec) note: computed in derived modprofiles
!     density (10**13 cm**-3 or 10**19 m**-3)
!     arho_exp and rmajor_exp (m)
!     power (MW)
!     flow (MW/kev=kA)
!
!   kev/sec per MW
!   kevdsecpmw=1.6022e-19*1.0e3*1.e-6
!
!       cgyrobohm_m(jm)=1.e-4*
!    >  9.79e5*(tem*1.e3)**.5/(arho_exp*100.)
!    >  *(1.02e2*(tem*1.e3)**.5/bt_exp/1.e4)**2*(amassgas_exp)**.5
!
!   sfactor(j)=
!     >      2.*pi_m*arho_exp*rho(j)*h_exp(j)*2.*pi_m*rmajor_exp
!   units of meter**2
!
!   9000 format (1i6,7e10.3)
 
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!mnt                         the model
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!mnt
!mnt supply model for chietem,chietim,chienem
!mnt                  chiitem,chiitim,chiinem
!mnt                  difftem,difftim,diffnem
!mnt
!mnt    chi_s and diff_s must be in meters**2/sec units
!mnt     and recall chi_s refer to total energy flow
!mnt
!mnt    if the model chi_s refer to "heat conduction" flow
!mnt    then a convection term xconv*3./2.*t_m*flow_exp is added.
!mnt    normally input xconv=0. otherwise xconv=1. or 5./3.
!mnt
!mnt    it is also possible to build convection into the model
!mnt    with "aconv".  aconv and xconv should not be double counted.
!mnt
!mnt note: can use diagonal forms with off diagonal dependence on
!mnt zpmte,zpmti,zpmne intrinsic to the diagonal elements as in sample
!mnt normall models are written in diagonal for with dependence on off
!mnt diagonal gradients implicit
!mnt
!mnt note: when flow is large anomalous e-i exchange should be added
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!xx      if (imodel.eq.8) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 2DGLF quasilinear model  GLF23
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!  test effect of canonical density gradients
! note pzmn_sol only in zpmne_q and zpmni_q which go into
! glf2d drivers rlne_gf and rlni_gf
 
!       zpmne_q=(1.+xparam_gf(20))*zpmne*pzmn_sol(jm)-xparam_gf(19)*
!     > (alog(q_exp(jm-1))-alog(q_exp(jm)))/(rho(jm-1)-rho(jm))
!       zpmni_q=(1.+xparam_gf(20))*zpmni*pzmn_sol(jm)-xparam_gf(19)*
!     > (alog(q_exp(jm-1))-alog(q_exp(jm)))/(rho(jm-1)-rho(jm))
 
       zpmne_q= zpmne
       zpmni_q= zpmni
 
       rmaj_gf=rmaj_exp(jm)/arho_exp
       rmin_gf=rmin_exp(jm)/arho_exp
       q_gf=q_exp(jm)
       betae_gf=max(cbetae*betae_m(jm), 1.E-6_R8)
       shat_gf=shat_exp(jm)*drhodrrrho(jm)
 
       if(xalpha.lt.0._R8) alpha_gf=-xalpha*alpha_m(jm)
       if(xalpha.gt.0._R8) alpha_gf=xalpha*alpha_exp(jm)
       if(alpha_gf.gt.4._R8) alpha_gf=4._R8
 
       elong_gf=elong_exp(jm)
 
       xnu_gf=cxnu*xnu_m(jm)
       taui_gf=tim/tem
       amassgas_gf=amassgas_exp
 
       apwt_gf=1._R8
! impurity dynamics not turned on by default
! and simple dilution included (idengrad=2, dil_gf=1-nim/nem)
! to turn on impurity dynamics need to change number of roots
! supply zimp_exp, amassimp_exp, and fractional density weights
! apwt_gf and aiwt_gf
       dil_gf=0._R8
       aiwt_gf=0._R8
!      zimp_gf=6.D0
!      amassimp_gf=12.D0
       zimp_gf=zimp_exp
       amassimp_gf=amassimp_exp
       rlnimp_gf=1._R8
       zpmnimp=1._R8
       if (idengrad.eq.2) dil_gf=1._R8-nim/nem
       if (idengrad.eq.3) then
         apwt_gf=nim/nem
         aiwt_jp1=(zeffm*nem-ni_m(jm+1)                                  &  
     &            -ns_m(jm+1))/(zimp_gf**2*nem)
         xnimp_jp1=aiwt_jp1*ne_m(jm+1)
         aiwt_gf=(zeffm*nem                                              &  
     &           -nim-ns_m(jm))/(zimp_gf**2*ne_m(jm))
         xnimp=aiwt_gf*ne_m(jm)
         zpmnimp=-(log(xnimp_jp1)-log(xnimp))/                           &  
     &           (rho(jm+1)-rho(jm))
         rlnimp_gf=zpmnimp*elong_exp(jm)**0.5_R8
       endif
!       write(*,66) rho(jm),nem,nim,rlnimp_gf,zpmnimp,amassimp_gf,
!     >             xnimp,zimp_gf,dil_gf
! 66    format(2x,0p1f4.2,0p9f10.5)
 
         rlte_gf=zpmte*drhodr(jm)
         rlti_gf=zpmti*drhodr(jm)
         rlne_gf=zpmne_q*drhodr(jm)
         rlni_gf=zpmni_q*drhodr(jm)
         rlnimp_gf=zpmnimp*drhodr(jm)
!        write(*,200) jm, rho(jm), rlti_gf, rlte_gf, rlne_gf, rlni_gf
 
        gamma_star_gf=vstarp_m(jm)
        gamma_e_gf=egamma_m(jm)
        gamma_p_gf=gamma_p_m(jm)
        if(itport_pt(4).eq.-1) then
          gamma_e_gf=egamma_exp(jm)
          gamma_p_gf=gamma_p_exp(jm)
        endif
!..jek 8/15/00
        if (irotstab.eq.0) then
          gamma_e_gf=egamma_exp(jm)
          gamma_p_gf=gamma_p_exp(jm)
          egamma_m(jm)=egamma_exp(jm)
        endif
!       gamma_mode_gf=gamma_mode_m(jm)
        gamma_mode_gf=0.0_R8
 
        if(i_delay.ne.0) then
!   i_delay should be negative on any intermediate step
 
         if(i_delay.gt.1) then
          do ii=1,i_delay-1
           egamma_d(jm,ii)=egamma_d(jm,ii+1)
          enddo
 
          egamma_d(jm,i_delay)=egamma_m(jm)
         endif
          gamma_e_gf=egamma_d(jm,1)
        endif
 
 
!.......THE  CALL TO GLF23
 
        call glf2d(iglf)
 
!.......POST CALL TO GLF23
 
!...diagnostic arrays (not presently used)
!
!      ky_j(jm)=xkyf_gf
!      gamma_j(jm,1)=gamma_gf(1)
!      gamma_j(jm,2)=gamma_gf(2)
!      gamma_j(jm,3)=gamma_gf(3)
!      gamma_j(jm,4)=gamma_gf(4)
!
!      freq_j(jm,1)=freq_gf(1)
!      freq_j(jm,2)=freq_gf(2)
!      freq_j(jm,3)=freq_gf(3)
!      freq_j(jm,4)=freq_gf(4)
!
!      phi_norm_j(jm,1)=phi_norm_gf(1)
!      phi_norm_j(jm,2)=phi_norm_gf(2)
!      phi_norm_j(jm,3)=phi_norm_gf(3)
!      phi_norm_j(jm,4)=phi_norm_gf(4)
!
!      do ik=1,ikymax_gf
!       gamma_k_j(ik,jm)= gamma_k_gf(1,ik)
!       freq_k_j(ik,jm) = freq_k_gf(1,ik)
!       chie_k_j(ik,jm) = chie_k_gf(ik)
!       chii_k_j(ik,jm) = chii_k_gf(ik)
!      enddo
 
       anrate_m(jm)=gamma_gf(1)
       anrate2_m(jm)=gamma_gf(2)
!      dnrate_m(jm)=0.
!      dtnrate_m(jm)=0.
       anfreq_m(jm)=freq_gf(1)
       anfreq2_m(j)=freq_gf(2)
!      dnfreq_m(jm)=0.
!       xkymax_m(jm)=xky_gf(1)
!       xkymax2_m(jm)=xky_gf(2)
 
       gfac=geofac(jm)
!
! exch_m in MW/m**3
!   kev/sec per MW
!k       kevdsecpmw=1.6022e-19*1.0e3*1.e-6
 
!k      exch_m(jm)=1.e19*
!k     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2*exch_gf*cmodel
      exch_m(jm)=1.E19_R8*1.6022E-19_R8*1.0E3_R8*1.E-6_R8*               &  
     & nem*tem*csda_m(jm)*rhosda_m(jm)**2*exch_gf*cmodel
 
      exchm=exch_m(jm)
 
! exch_m is directly related to the flow
! for a single mode branch exch_gf=-(-freq_gf(1)/xkyf_gf)*diff_gf*rln_gf.
! we can  not expect to find exch_m without knowing flow_exp as input.
! and solving self consistant flow eq. flown=flow_exp for density
! density profile.
!
! however, knowing freq_gf(1) from the gf model we can compute exch_exp
! from flow_exp using
!       flowm=kevdsecpmw*1.*nem*1.e19/arho_exp*gradrhosq_exp(jm)*
!     >       sfactor(jm)*(difftem*zpmte+difftim*zpmti+diffnem*zpmne)
! we have:
 
!       diffgb_local=flow_exp(jm)/
!     > (kevdsecpmw*1.*nem*1.e19/arho_exp*gradrhosq_exp(jm)*sfactor(jm)*
!     > zpmne_q)/cgyrobohm_m(jm)
 
!       exchgb_local=-(-freq_gf(1)/xkyf_gf)*diffgb_local*rlni_gf
 
!       exch_exp(jm)=1.e19*
!     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2*exchgb_local
 
!       exch_exp(jm)=flow_exp(jm)*tem*(-1.)*(-freq_gf(1)/xkyf_gf)*
!     > sqrt(elong_exp(jm))*arho_exp/gradrhosq_exp(jm)/sfactor(jm)
!     >/arho_exp(jm)**2
 
 
!   note electron(ion) wave freq > 0(<0) cool(heat) electrons
! (-1) denotes electron to ion
 
! to emphasize, we can not know exch_exp better than we know flow_exp
 
! chietem, chiitim, diffen in m**2/sec
! ITER definition of "chi" assumes will be proceeded with gradrhosq factor
 
      chietem=cmodel*gfac*chie_gf*cgyrobohm_m(jm)
      chiitim=cmodel*gfac*chii_gf*cgyrobohm_m(jm)
      diffnem=cmodel*gfac*diff_gf*cgyrobohm_m(jm)
 
      etaphim=cmodel*gfac*eta_phi_gf*cgyrobohm_m(jm)
      etaparm=cmodel*gfac*eta_par_gf*cgyrobohm_m(jm)
      etaperm=cmodel*gfac*eta_per_gf*cgyrobohm_m(jm)
 
!xx endif
!
      if ( itport_pt(1) .eq. 0 ) diffnem = 0._R8
      if ( itport_pt(2) .eq. 0 ) chietem = 0._R8
      if ( itport_pt(3) .eq. 0 ) chiitim = 0._R8
      if ( (itport_pt(4) .eq. 0) .and. (itport_pt(5) .eq. 0) ) then
         etaphim=0._R8
         etaparm=0._R8
         etaperm=0._R8
      endif
!
      diff_m(jm)=diffnem
      chie_m(jm)=chietem
      chii_m(jm)=chiitim
 
      etaphi_m(jm)=etaphim
      etapar_m(jm)=etaparm
      etaper_m(jm)=etaperm
!
!     write(6,*) jm,rho(jm),zpmte, zpmti, zpmne, zpmni
 
        enddo
 
!      if(jmm.eq.0)  write(6,*) 'jmm=', jmm,'going out'
!
!dmc        write(6,7801) diffnem,chietem,chiitim,etaphim,etaparm,etaperm,
!dmc     >     exchm
!dmc 7801   format(//' OUTPUTS...'/
!dmc     >'  diffnem = ',1pe17.10,' chietem = ',1pe17.10/
!dmc     >'  chiitim = ',1pe17.10,' etaphim = ',1pe17.10/
!dmc     >'  etaparm = ',1pe17.10,' etaperm = ',1pe17.10/
!dmc     >'  exchm = ',1pe17.10)
!
!dmc        call echo('diff_m output',diff_m,jmaxm)
!dmc        call echo('chie_m output',chie_m,jmaxm)
!dmc        call echo('chii_m output',chii_m,jmaxm)
!dmc        call echo('etaphi_m output',etaphi_m,jmaxm)
!dmc        call echo('etapar_m output',etapar_m,jmaxm)
!dmc        call echo('etaper_m output',etaper_m,jmaxm)
!dmc        call echo('exch_m output',exch_m,jmaxm)
!
!dmc        call echo('egamma_m output',egamma_m,jmaxm)
!dmc        call echo('egamma_d output',egamma_d,jmaxm)
!dmc        call echo('gamma_p_m output',gamma_p_m,jmaxm)
!dmc        call echo('anrate_m output',anrate_m,jmaxm)
!dmc        call echo('anrate2_m output',anrate2_m,jmaxm)
!dmc        call echo('anfreq_m output',anfreq_m,jmaxm)
!dmc        call echo('anfreq2_m output',anfreq2_m,jmaxm)
!
       return
       end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
