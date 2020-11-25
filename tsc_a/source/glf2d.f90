!#include "f77_dcomplx.h"
           subroutine glf2d(iglf)
 
! 12mar2003 jek added rms-theta for retuned version
! 02feb2001 jek fixed section w/ xparam(24) for power law ExB
! 08may2000 jek merged ITG and ETG loops over ky for better MPI optimization
! 23nov1999 jek added zgeev eigenvalue solver (eigen_gf=2)
! 29mar1999 fgtok -s cgg.table "dmc:  rename eispack routines"
! 29mar1999 fgtok -s rr.table "dmc:  rename intrinsic REAL -> REAL"
!
! dmc -- Cray/workstation portable real*8<-->complex*16 conversion routines
 
!#include "f77_dcomplx.h"
 
!
!glf2d.f 12-mar-03 Kinsey
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
 
!
!***********************************************************************
! questions  should be addressed to
!  Ron Waltz 619-455-4584  or email: waltz@gav.gat.com
!***********************************************************************
 
 
! 2D GLF equations with massless isothermal passing electrons from
!  Waltz et al, Phys. of Plasmas 6(1995)2408
 
! In Eq. 24  p_u_par is replaced by n_u/(1-reps) +t_u and
! the isothermal conditions t_u= (betae/2)*w_s/k_par*(rlt)*a_par is
! used. Thus n_u = (1-reps)*ph (adiabatic passing electrons) for betae=0
!
! In Eq. 23 (p_u_par+p_u_per) is replaced by n_u+t_u
! using the isothermal condition the mHD beta limit is too low by
! 1/(1+ reps(-0.25+0.75rlte)/(taui*(rln+rlti)+(rln+rlte))
! It is possible to patch this up by replacing with x*n_u+y*t_u
! then solving for x and y to obtain
! universal MHD betae_crit=2*k_par**2/(w_d*(taui*(rln+rlti)+w_d*(rln+rlti)))
! beta_crit=(1+taui)*betae_crit=(1/2)(1/q**2)*(L_p/rmajor)
! 1/2 is replaced by s_hat in models with shear
 
! EVERYTHING else including normalizing units follows the paper.
 
!  unit of microlength is rho_s, unit of macrolength is a
!   a is typically rho value at separatrix
!  unit of time is a/c_s; unit of diffusion is (c_s/a)*rho_s**2
!  c_s=sqrt(Te/M_i), omega=eB/cM_i, rho_s=c_s/omega
 
! example balance equations to clarify meaning  of diffusivities
!
!       chie_hat is effective energy diffusivity
!
!  (3/2) d ne Te/ dt =
! -1/V(rho)_prime d /d rho V(rho)_prime
!             |grad rho|**2 ne chie_hat (c_s rho_s**2/a)(-d Te/ d rho)
!      -exch_hat (ne Te) c_s/a (rho_s/a)**2 +heating density
!
! and similarly for ion equation
! note that no convective part is added, ie "convection" is included
! inside chie_hat
! note no impurity energy flow is computed although it could be easily done
 
!        d_hat is effective plasma diffusivity for ions
!
!        d ni / dt =
! -1/V(rho)_prime d /d rho V(rho)_prime
!               |grad rho|**2 ne d_hat (c_s rho_s**2/a) (-d ni/ d rho)
!        + plasma source density
 
!        d_im_hat is the effective plasma diffusivity for impurities
 
!        eta_phi is the toroidal viscosity or momentum diffusivity
!
!       M ni d v_phi/ dt =
! error found 5/12/98  should be d (M ni v_phi)/ dt =
! -1/V(rho)_prime d /d rho V(rho)_prime
!      |grad rho|**2 ( M ni eta_phi_hat (c_s rho_s**2/a) (-d v_phi/ d rho)
!                    +    M v_phi ne d_hat (c_s rho_s**2/a) (-d ne/ d rho))
!        + toroidal momentum source density
!
! note that a convective part was added
!
!  eta_par_hat and eta_per_hat are diagnostic. See CAUTION on eta_phi_hat
!  at large gamma_p=  (-d v_phi/ d rho) /(c_s/a)
!
!  chie_e_gf is the eta_e mode electron transport which is te <-> ti
!  and mi <-> me isomorphic to eta_i (ITG) ion transport
!  with adiabatic electrons.
!  these mode obtain at high-k where the ions are adiabatic from
!  the gyro cut-off.
!  their wave numbers are sqrt(mi/me) larger than ITG modes and
!  since their frequencies are sqrt(mi/me) larger, they are not
!  rotationally shaer satbilized.
!  when xparam_gf(10).eq.0 xparam_gf(10)*chie_e_gf is added to
!  chie_gf and chie_e_gf is a diagnostic.
 
! input
 
!  eigen_gf = 0 use cgg eigenvalue solver (default)
!           = 1 use generalized tomsqz eigenvalue solver
!           = 2 use zgeev eigenvalue solver
!  nroot number of equations
!  iflagin(1:20) control flags
!   iflagin(1) 0 use ky=ky0; 1 use landau damping point
!   iflagin(2) 0. local w_d and k_par "2d"; 1 fit to trial function "3d"
!   iflagin(3) 0,1,and 2 fix up park low high beta and beta lim elong factor
!   iflagin(4) 0 trapped electron Waltz EoS 1 weiland EoS
!   iflagin(5) rms_theta 0:fixed; 1 inverse to q/2 ; 2 inverse to root q/2
!                        3: inverse to xparam(13)*(q/2-1)+1.
!              5 for retuned rms-theta
!  xparam(1:20) control parameters
!   xparam(1:2): idelta=xi*xparam(1)+xparam(2) nonadiabatic electron response
!   xparam(3) multiplier park_gf(high betae)/ park_gf(low betae) -1
!   xparam(6)+1. is enhancement of xnueff
!   xparam(7) coef of resistivity
!   xparam(8) cut off on rotational stabilization
!   xparam(9)+1. is shape (triangularity) enhancement to beta_crit
!   xparam(10) is high k electron mode enhancement
!   xparam(11:12) lamda parameters
!   xparam(13) rms_theta q-dependence
!   xparam(14)  adjustment to gamma_p avoiding negative viscosity
!   xparam(15)   (1+xparam(15)*reps trapped electron fraction
!   xparam(16) rms_theta shat dependence
!   xparam(17) ""
!   xparam(18) rms_theta betae dependence
!   xparam(19:20)  extra
!   xparam(21) 1 add impurity energy diffusivity to ion energy diffusivity
!   xparam(22) >0 keeps gamma_e from changeing spectrum
!   xparam(23) 1. kills kx**2 in k_m**2
!   xparam(24) exb damping model
!  ky0=k_theta*rho_s; k_theta= nq/r; normally 0.3
!  rms_theta width of phi**2 mode function for best fit near pi/3
!  rlti=a/L_Ti   a/L_f= sqrt(kappa) a d ln f / d rho
!  rlte=a/L_Te
!  rlne= a/L_ne
!  rlni= a/L_ni
!  rlnimp= a/L_nim
!  dil=1.-ni_0/ne_0  dilution
!  apwt = ni_0/ne_0
!  aiwt = nim_0/ne_0
!  taui=Ti/Te
!  rmin=r/a
!  rmaj=Rmaj/a
!  xnu=nu_ei/(c_s/a)
!  betae=neTe/(B**2/(8pi))  0 is electrostatic
!  shat= dlnr/drho used only for parallel dynamics part
!  alpha local shear parameter or MHD pressure grad (s-alpha diagram)
!  elong= local elongation or kappa
!  xwell amount of magnetic well xwell*min(alpha,1)
!  park=1  (0) is a control parameter to turn on (off) parallel motion
!       0.405 best at zero beta and 2.5x larger at high beta..see iflagin(3)
!  ghat=1  (0) is a control parameter to turn on (off) curvature drift
!  gchat=1 (0) is a control parameter to turn on (off) div EXB motion
!  adamp= radial mode damping exponent  1/4 < adamp < 3/4
!       0.25 from direct fit of simulations varying radial mode damping
!   but 0.75 is better fit to rlti dependence
!  alpha_star O(1-3)  gyyrobohm breaking coef for diamg. rot. shear
!  gamma_star ion diamagnetic rot shear rate in units of c_s/a
!  alpha_e O(1-3)   doppler rot shear coef
!  gamma_e    doppler rot shear rate in units of c_s/a
!  alpha_p 1.5  fit for parallel velocity shear effect at rmaj=3 and q=2
!  gamma_p    parallel velocity shear rate (-d v_phi/ drho) in units of c_s/a
!  kdamp model damping normally 0.
 
! output
 
!  yparam(20) output diagnostics
! kyf  value of ky used
! gamma   leading mode growth rate in c_s/a
! freq    leading mode freq rate in c_s/a
! ph_m    (e phi /T_e)/(rho_s/a)  saturation value
! d_hat    plasma diffusivity for ions
! d_im_hat    plasma diffusivity for impurities
! chii_hat ion energy diffusivity
! chie_hat electron energy diffusivity
! exch_hat anomalous e to i energy exchange
! eta_par_hat parallel component of toroidal momentum diffusivity
! eta_per_hat perpendicular    ""
! eta_phi_hat toroidal momentun diffusivity
 
! internal definitions
! nroot = number of equations,
!   nroot=12 full impurity dynamics
!   nroot=9 exb convective impurity dynamics
!   nroot=8 full pure plasma, nrout=6 (betae=0), nrout=5 (betae=0 and park=0)
! v(i)  12 solution vector
!   v(1)=n_i,v(2)=p_par,v(3)=p_per,v(4)=n_t,v(5)=p_t
!   v(6)=u_par, v(7)=n_u, v(8)=a_par
!   v(9)=n_im, v(10)=p_im_par,v(11)=p_im_per,v(12)=u_im_par
! -i*omega v(i)= sum_j amat(i,j) v(j) where omega=freq+xi*gamma
! quasineitrality is
!  (-idelta+(1-dil)*(1/taui)*(1-g1))*ph=f0*ph=(1-dil)*n_i-n_t-n_u
!  or (1-idelta-reps+(1-dil)*(1/taui)*(1-g1))*ph=f1*ph=(1-dil)*n_i-n_t
!  for betae=0
! k_m  inverse mixing length
 
! numerical controls and flags
!
!
!...Dimensions
!
! neq maximum number of equations
! ieq actual number used
 
!***********************************************************************
!***********************************************************************
!
!
      USE GLF
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
! Glf is common block, which must contain all the _gf inputs and outputs
!
      character*1 jobvr, jobvl
      integer neq, iflagin(30), ilhmax, ilh, ikymaxtot,                  &  
     &  lprint, ieq, j1, j2, j, i, jmax,                                 &  
     &  iar, ifail, jroot(4), itheta, iky, iky0, iroot, iglf
      REAL*8 epsilon
      REAL*8 AREAL
      parameter ( neq = 12, epsilon = 1.E-34_R8)
!
      REAL*8 pi, xparam(30),yparam(2*nmode),                             &  
     &  nroot,ky0,rms_theta,rlti,rlte,rlne,rlni,dil,taui,                &  
     &  rmin,rmaj,q,rlnimp,amassimp,zimp,mimp,                           &  
     &  aikymax, aiky, apwt, aiwt,                                       &  
     &  alpha_mode, gamma_mode, alpha_p, gamma_p,                        &  
     &  amassgas, chi_par_1, chi_per_1, x43,                             &  
     &  anorm, ave_g, ave_g0,                                            &  
     &  ave_cos, ave_theta2, ave_kxdky2,                                 &  
     &  dtheta, theta, phi2, ave_k_par, chk,                             &  
     &  alpha_n, yk, byk, fnn, fnp, fnf, fpn, fpp, fpf,                  &  
     &  xnu_col, amass_e, chkf, xadiabat,                                &  
     &  eta_par_hat, eta_per_hat, anorm_k, del_k,                        &  
     &  gamma_k_max, gamma_gross_net,                                    &  
     &  xnu,betae,shat,alpha,elong,xwell,                                &  
     &  park,ghat,gchat,kdamp,                                           &  
     &  adamp,alpha_star,gamma_star,alpha_e,gamma_e,                     &  
     &  kyf,gamma,freq,ph_m,d_hat,d_im_hat,                              &  
     &  chii_hat,chie_hat,exch_hat
      COMPLEX*16 xi, idelta,                                             &  
     &  v(1:12), amat(1:12,1:12),                                        &  
     &  n_i,p_par,p_per,n_t,p_t,u_par,n_u,a_par,ph,t_u,n_e,              &  
     &  n_im,p_im_par,p_im_per
!     complex u_im_par
      REAL*8 b0,g0,g1,g2,g3,g12,g23,                                     &  
     &  b0i,g0i,g1i,g2i,g3i,g12i,g23i
      COMPLEX*16 f0,f1
      REAL*8 k_par,ky,kx,k_per,k_m,                                      &  
     &  w_s, w_d, w_d0, w_cd,                                            &  
     &  reps,xnueff,betae0,k_par0
      COMPLEX*16 xmu,lamda_d,                                            &  
     &  xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per
      REAL*8 gam_par,gam_per,x_par,x_per,xt_mhd,yt_mhd,                  &  
     &  th,tc,fh,fc,                                                     &  
     &  phi_norm,gamma_r
      COMPLEX*16 chknu,chknt,chknt2
      REAL*8 phi_renorm,gamma_net
!
!...Declarations for eigenvaluesolver
!
      REAL*8 zgamax
!
!... solver varaibles
!
      parameter ( iar=neq )
!     integer iai, ivr, ivi, intger(neq) ! if NAG solver f02ake used
!     parameter ( iai=neq, ivr=neq, ivi=neq )
      REAL*8 ar(iar,neq), ai(iar,neq), rr(neq), ri(neq)                  &  
     &  , vr(iar,neq), vi(iar,neq)
      REAL*8 br(iar,neq), bi(iar,neq), beta_tom(neq), ztemp1
 
      integer matz
      REAL*8 fv1(neq),fv2(neq),fv3(neq)
!
! amat(i,j) = complex matrix A
! zevec(j) = complex eigenvector
!
      integer lwork
      parameter ( lwork=198 )
      COMPLEX*16 mata(iar,neq),cvr(iar,neq),cvl(iar,neq),w(neq)
      COMPLEX*16 work(lwork)
      REAL*8 rwork(2*neq)
      COMPLEX*16 zevec(neq,neq), zomega(neq)
      REAL*8 gammaroot(4),freqroot(4),phi_normroot(4)
!
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
!
!   return zero if no unstable modes
!
      if(ipert_gf.eq.1.and.ngrow_k_gf(0).eq.0)go to 888
!
!...initialize variables
!
! inputs.........................................................
 
      do i=1,30
       iflagin(i)=iflagin_gf(i)
       xparam(i)=xparam_gf(i)
      enddo
 
      ilhmax=1
      ikymaxtot=ikymax_gf
!     if (xparam_gf(10).gt.0.) ilhmax=2
!
! If ETG modes included, then double ky spectrum
! Inside ky loop, ilh=1 low k ion modes and ilh=2 high k electron modes
! For ETG mode, use complete te <-> ti and mi <-> me isomorphism
! with adiabatic electron ITG mode then chii_hat will be electron
! transport in gyrobohm electron units with T_i.
! chie_e_gf converted back to c_s*rho_s**2/a units and added
! to chie_gf after ky loop
!
      if (xparam_gf(10).gt.0._R8) then
        ilhmax=2
        ikymaxtot=2*ikymax_gf
      endif
!
      nroot=nroot_gf
      ky0=xky0_gf
      rms_theta=rms_theta_gf
      rlti=rlti_gf
      rlte=rlte_gf
      rlne=rlne_gf
      rlni=rlni_gf
      rlnimp=rlnimp_gf
      dil=dil_gf
      apwt=apwt_gf
      aiwt=aiwt_gf
      taui=taui_gf
      rmin=rmin_gf
      rmaj=rmaj_gf
      q=q_gf
      xnu=xnu_gf
      betae=betae_gf
      shat=shat_gf
      alpha=alpha_gf
      elong=elong_gf
      xwell=xwell_gf
      park=park_gf
      ghat=ghat_gf
      gchat=gchat_gf
      adamp=adamp_gf
      alpha_star=alpha_star_gf
      gamma_star=gamma_star_gf
      alpha_e=alpha_e_gf
      gamma_e=gamma_e_gf
      alpha_mode=alpha_mode_gf
      gamma_mode=gamma_mode_gf
      alpha_p=alpha_p_gf
      gamma_p=gamma_p_gf
      kdamp=xkdamp_gf
      lprint=lprint_gf
      amassgas=amassgas_gf
      amassimp=amassimp_gf
      zimp=zimp_gf
      if(ipert_gf.eq.0)then
       do j=0,nmode
        ngrow_k_gf(j) = 0
       enddo
      endif
!
      idelta=0._R8
!     if(ilh.eq.1) idelta=xi*xparam(1)+xparam(2)
 
!.................................................................
!
      if (lprint.gt.0) open(1)
      ieq  = nroot
!
      if (lprint.eq.99) then
      write(1,*) 'ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q: ',  &  
          ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q
      write(1,*)'xnu,beta,shat,alpha,elong,xwell: ',                     &  
          xnu,betae,shat,alpha,elong,xwell
      write(1,*)'park, ghat, gchat: ',                                   &  
          park, ghat, gchat
      write(1,*)'adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp: ',   &  
          adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp
 
      endif
      if (lprint.eq.98) then
      write(2,*) 'ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q: ',  &  
          ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q
        write(2,*)'xnu,betae,shat,alpha,elong,xwell: ',                  &  
          xnu,betae,shat,alpha,elong,xwell
        write(2,*)'park, ghat, gchat: ',                                 &  
          park, ghat, gchat
        write(2,*)'adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp: ',  &  
          adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp
 
      endif
!
      xi=(0._R8,1._R8)
      pi=atan2 ( 0.0_R8, -1.0_R8)
 
! GLF model coefficients
 
      chi_par_1=2._R8*sqrt(2._R8/pi)
      chi_per_1=sqrt(2._R8/pi)
 
      gam_par=3._R8
      gam_per=1._R8
      x_par=2._R8
      x_per=3._R8/2._R8
 
      xmu=(0.80_R8+.57_R8*xi)
 
      xnu_par_par=(1._R8+xi)
      xnu_par_per=0._R8
      xnu_per_par=0._R8
      xnu_per_per=(1._R8+xi)
 
      lamda_d=(-0.7_R8-0.80_R8*xi)
      if(xparam(11).ne.0._R8.or.xparam(12).ne.0._R8)                     &  
     &  lamda_d=xparam(11)-xi*xparam(12)
      x43=1._R8
      if (iflagin(4).eq.1) then
       lamda_d=5._R8/3._R8
       x43=4._R8/3._R8
      endif
!
! 3d trial wave function analysis
!
      if(iflagin(2).ge.1) then
       if(rms_theta.eq.0._R8) rms_theta=pi/3._R8
       if(iflagin(5).eq.1) rms_theta=rms_theta_gf*(2._R8/q_gf)
       if(iflagin(5).eq.2) rms_theta=rms_theta_gf*(2._R8/q_gf)**0.5_R8
       if(iflagin(5).eq.3) rms_theta=rms_theta_gf/                       &  
     &      (xparam(13)*(q_gf/2._R8-1._R8)+1._R8)                        &  
     &  /sqrt(1._R8+xparam(16)*(shat_gf**2-1._R8)+xparam(17)*            &  
     &    (shat_gf-1._R8)**2)                                            &  
     &  /(1._R8+xparam(18)*sqrt(betae/.006_R8))
       if(iflagin(5).eq.4) rms_theta=rms_theta_gf/                       &  
     &      (xparam(13)*(q_gf/2._R8-1._R8)+1._R8)                        &  
     &  /sqrt(1._R8+xparam(16)*(shat_gf**2-1.0_R8)+xparam(17)*           &  
     &    (shat_gf-1._R8)**2+xparam(19)*(alpha-0.5_R8)**2._R8)           &  
     &  /(1._R8+xparam(18)*sqrt(betae/.006_R8))
       if(iflagin(5).eq.5) rms_theta=rms_theta_gf/                       &  
     &      (xparam(13)*((q_gf/2._R8)-1._R8)+1._R8)                      &  
     &  /sqrt(1._R8+xparam(16)*((shat_gf-                                &  
     &  xparam(19)*alpha)**2-0.5_R8)+xparam(17)*                         &  
     &  (shat_gf-xparam(19)*alpha-0.5_R8)**2)/                           &  
     &  taui_gf**0.25_R8
! along the filed line physics with wave function
! phi=exp(-theta**2/(4.*rms_theta**2))=W_even
! ave_F= [int(0 to inf) F phi**2 d_theta]/[int(0 to inf)  phi**2 d_theta]
! ave_theta**2=(rms_theta)**2
 
! phi, densities and pressures are even functions of theta
! the odd functions like u_par can be represented by
! W_odd= W*i*theta/rms_theta*W_even
! then for W=-1, the k_par operator = i*k_par=1/(rmaj*q) d/dtheta
! becomes 1/(rmaj*q)/(2.*rms_theta)*park  (park=1) in every equation
! ie ave_k_par=1/(2.*rms_theta)
! park is tuned to best fit.
!
! parallel velocity shear gamma_p breaks parity so wave functions
! become mixed W_even->(1-xi*gamma_p*alpha_n*theta/rms_theta)*W_even
!
! gamma_p*alpha_n mustbe dimensionless and independent of norm a
! hence alpha_n=(rmaj/3)*alpha_p since gamma_p is in units
! of c_s/a  and rmaj=rmajor/a. Since in the slab limit where
! where we can have the parallel shear drive, rmaj enters only with
! parameter rmaj*q, we further assume
! alpha_n=(rmaj/3.)*(q/2)*alpha_p as the appropriate scaling
! q=2 and rmaj=3 are the norm points for alpha_p=1.5
! For the extreme toroidal limit q-> infinity where rmaj and
! q are not product associated, we will lose the instability.
!
! to first order in gamma_p*alpha_n
! this leads to a weighting  factor [gamma_p*alpha_n] in the
! xi*ky*ph1*gamma_p linear drive and in the
! eta_phi_hat=conjg(u_par)*(-xi*ky*ph)/gamma_p toroidal vocosity.
!
! the correct dependence gamma-gamma0 going like gamma_p**2 is found
! but QLT eta_phi_hat goes like gamma_p**2 also
! CAUTION: this is worked out only to first order in gamma_p*alpha_n
! small.  It seems likely that there is a higher order saturation factor
! something like 1/(1+(gamma_p*alpha_n)**2) in  eta_phi_hat
!
! doppler (EXB) rotational shear also breaks parity
! thus there should a term egamma*alpha_n_e added to gamma_p*alpha_n
! but it is unclear how to weight alpha_n_e compared to alpha_n
!
! see R.R.Dominguez and G.M Staebler Phys. Fluids B5 (1993) 3876
! for a discussion of QLT theory of anomalous momentum transport in
! slab geometry
 
! compute weight factors
! fix later so these are computed only once per j grid point
 
!      ave_theta2=rms_theta**2
 
      anorm=0._R8
      ave_g=0._R8
      ave_g0=0._R8
      ave_cos=0._R8
      ave_theta2=0._R8
      ave_kxdky2=0._R8
 
      dtheta=4._R8*rms_theta/100._R8
      theta=0._R8
!
      do itheta=1,100
       theta=theta+dtheta
       phi2=exp(-theta**2/(2._R8*rms_theta**2))
       anorm=anorm+phi2*dtheta
       ave_theta2=ave_theta2+                                            &  
     &  theta**2*phi2*dtheta
       ave_g=ave_g +                                                     &  
     & (-xwell*min(1._R8,alpha)+cos(theta)+                              &  
     &    (shat*theta-alpha*sin(theta))*sin(theta))*phi2*dtheta
       ave_g0=ave_g0 + phi2*dtheta
       ave_kxdky2=ave_kxdky2+                                            &  
     &  (abs(shat*theta-alpha*sin(theta)))**2*phi2*dtheta
       ave_cos=ave_cos +                                                 &  
     &  cos(theta)*phi2*dtheta
      enddo
!
      ave_theta2=ave_theta2/anorm
      ave_g=ave_g/anorm
      ave_g0=ave_g0/anorm
      ave_kxdky2=ave_kxdky2/anorm
      ave_cos=ave_cos/anorm
!
      ave_k_par=1/(2._R8*rms_theta)
 
      chk=abs(ave_theta2-rms_theta**2)/rms_theta**2
      if (chk.gt..02_R8) write (6,*) 'chk:', chk
 
      alpha_n=(rmaj/3._R8)*(q/2._R8)*alpha_p
 
      if(lprint.eq.2) then
       write(6,*) 'rms_theta,chk :', rms_theta, chk
       write(6,*) 'ave_theta2,ave_g,ave_k_par,ave_cos:',                 &  
     &   ave_theta2,ave_g,ave_k_par,ave_cos
      endif
      endif
      if(iflagin(2).eq.0) then
       shat=1._R8
       ave_theta2=1._R8
       ave_g=1._R8
       ave_g0=1._R8
       ave_kxdky2=1._R8
       ave_k_par=1._R8
       ave_cos=1._R8
      endif
!
! start ky loop
! first half ITG, second half high-k ETG ... each with ikymax_gf modes
! ilh=1 low k ion modes  ilh=2 high k electron modes
 
      do iky0=1,ikymaxtot
!
!gms      iky=iky0
      iky = ikymax_gf+1-iky0
      ilh=1
!
! offset iky if in high-k range and set ilh=2
!
      if (iky0.gt.ikymax_gf) then
!gms         iky=iky0-ikymax_gf
         iky = ikymaxtot+1-iky0
         ilh=2
      endif
!
      if (ilh.eq.2) then
       nroot=6
       ieq=nroot
       xnu=0._R8
       betae=1.E-6_R8
       rlte=rlti_gf
       rlti=rlte_gf
       rlne=rlni_gf
       rlni=rlne_gf
       rlnimp=epsilon
       dil=1._R8-1._R8/(1._R8-dil_gf)
       apwt=1._R8
       aiwt=0._R8
       taui=1._R8/taui_gf
       rmin=epsilon
       xparam(7)=0._R8
       xparam(6)=-1._R8
       alpha_star=0._R8
       alpha_e=0._R8
       alpha_p=0._R8
       alpha_n=0._R8
       alpha_mode=0._R8
! check this for current driven mode
      endif
!
      idelta=0._R8
      if (ilh.eq.1) idelta=xi*xparam(1)+xparam(2)
!
! logarithmic ky grid
!
      if(ikymax_gf.gt.1) then
       aikymax=ikymax_gf
       aiky=iky
       yk=aiky/aikymax
 
       byk=log(xkymax_gf/xkymin_gf)/(1._R8-1._R8/aikymax)
       ky=xkymax_gf*exp(byk*(yk-1._R8))
      endif
      if(ikymax_gf.eq.1) then
!     ky=sqrt(2.*taui)/rlti/(rmaj*q)
!     from w_star_ti=v_i_th*k_par
!  possible physics basis of q (ie current) scaling ..to be determined
       if(iflagin(1).eq.0) ky=ky0
       if(iflagin(1).eq.1) ky=ky0*sqrt(taui)*(3._R8/rlti)*               &  
     &    (3._R8*2._R8/rmaj/q)
       if(iflagin(1).eq.2) ky=ky0*(2._R8/q)
       if(ky0.eq.0._R8) ky=0.3_R8
      endif
 
 
      kyf=ky
!
 
      kx=ky*sqrt(ave_kxdky2)
      k_per=sqrt(ky**2+kx**2)
      k_m=sqrt(ky**2+(1._R8-xparam(23))*kx**2) ! inverse mixing length model
 
       do iroot=1,4
        gammaroot(iroot)=0._R8
        freqroot(iroot)=0._R8
        phi_normroot(iroot)=0._R8
       enddo
        d_hat=0._R8
        d_im_hat=0._R8
        chie_hat=0._R8
        chii_hat=0._R8
        exch_hat=0._R8
        eta_par_hat=0._R8
        eta_per_hat=0._R8
        jroot(1)=0
        jroot(2)=0
        jroot(3)=0
!
! skip this k for perturbation if last call was stable
!
      if(ipert_gf.eq.1.and.ngrow_k_gf(iky0).eq.0)go to 777
! skip this k if the previous k was stable and 4 k's have been done
      if(iky.lt.ikymax_gf-4.and.ngrow_k_gf(iky0-1).eq.0)go to 777
 
! primary ions
 
      b0=taui*k_per**2
 
!     Pade aproximates...may use gamma functions later
      g0=1._R8
      g1=1._R8/(1+b0)
      g2=1._R8/(1+b0)*g1
      g3=1._R8/(1+b0)*g2
 
      g12=(g1+g2)/2._R8
      g23=(g2+g3)/2._R8
 
! impurity ions
 
      b0i=taui*k_per**2*amassimp/amassgas/zimp**2
 
!     Pade aproximates...may use gamma functions later
      g0i=1._R8
      g1i=1._R8/(1+b0i)
      g2i=1._R8/(1+b0i)*g1i
      g3i=1._R8/(1+b0i)*g2i
 
      g12i=(g1i+g2i)/2._R8
      g23i=(g2i+g3i)/2._R8
 
      mimp=amassimp/amassgas
 
 
      w_s=ky
      w_d=(ghat*2._R8/rmaj)*ky*ave_g
      w_d0=(ghat*2._R8/rmaj)*ky*ave_g0
      w_cd=(gchat*2._R8/rmaj)*ky*ave_g
 
      k_par=park/(rmaj*q)*ave_k_par*sqrt((1._R8+elong**2)/2._R8)
 
!     sqrt((1.+elong**2)/2.) to get higher beta_crit prop to k_par**2
!     roughly same as betae-> betae/((1.+elong**2)/2.)
!     physically like shortening the connection length to good curv.
 
      if (iflagin(3).eq.2) then
       betae=betae_gf/(1._R8+xparam(3))**2/(1._R8+xparam(9))
      endif
 
 
      if (iflagin(3).eq.1) then
       k_par=park/(rmaj*q)*ave_k_par
       betae=betae_gf/(1._R8+xparam(3))**2/(1._R8+xparam(9))             &  
     &      /((1._R8+elong**2)/2._R8)
      endif
 
!     we put the park enhancement directy into betae
!     ie park=.5 best at low beta and 2.5x.5=1.25 at high beta
 
!     option iglagin(3)=1 puts beta_crit elongation enhancement
!     directly into betae
 
!     option iflagin(3)=2 puts beta_crit elongation factor into
!     the connection length
 
!     an extra shape  factor 2 (triangularity) enhancement
!     is optained by (1.+xparam(9))=2.
!     if w_d negative flip sign of dissipative parts
      if(w_d.lt.0._R8) then
! error 12/21       lamda_d=-conjg(lamda_d)
       lamda_d=conjg(lamda_d)
 
       xmu=-conjg(xmu)
 
       xnu_par_par=-conjg(xnu_par_par)
       xnu_par_per=-conjg(xnu_par_per)
       xnu_per_par=-conjg(xnu_per_par)
       xnu_per_per=-conjg(xnu_par_per)
      endif
 
      reps=(1._R8+xparam(15))*                                           &  
     &   sqrt((rmin/rmaj)*(1._R8+ave_cos)/(1._R8+(rmin/rmaj)*ave_cos))
 
      if(nroot.le.3) reps=0._R8
 
! fix trapped eletron MHD limit
! 3/4*reps*(1+rlte) + xt_mhd*(1-reps)+yt_mhd*rlte=1+rlte
! solve for xt_mhd and yt_mhd
! 3/4*reps+yt_mhd=1; 3/4*reps+xt_mhd*(1-reps)=1
 
      yt_mhd=(1-x43*(3._R8/4._R8)*reps)
      xt_mhd=(1._R8-x43*(3._R8/4._R8)*reps)/(1._R8-reps)
 
! collision detrapping retrapping model
 
      xnueff=(1._R8+xparam(6))*xnu/(reps**2+1.E-6_R8)
 
! very difficult get xnueff correct hince add enhancement factor
! and fit to Beer or GKS
 
      th=4.08_R8
      tc=0.918_R8
      fh=0.184_R8
      fc=0.816_R8
 
      fnn=xnueff*((th/tc**(3._R8/2._R8))-(tc/th**(3._R8/2._R8)))/(th-tc)    
      fnp=xnueff*(3._R8/2._R8)*                                          &  
     &   ((1._R8/th)**(3._R8/2._R8)-(1._R8/tc)**(3._R8/2._R8))/(th-tc)
      fnf=xnueff*((fh/th**(3._R8/2._R8))+(fc/tc**(3._R8/2._R8)))
 
      fpn=xnueff*(2._R8/3._R8)*                                          &  
     &   ((th/tc**(1._R8/2._R8))-(tc/th**(1._R8/2._R8)))/(th-tc)
      fpp=xnueff*((1._R8/th)**(1._R8/2._R8)-(1._R8/tc)**(1._R8/2._R8))/  &  
     & (th-tc)
      fpf=xnueff*(2._R8/3._R8)*((fh/th**(1._R8/2._R8))+(fc/tc**(1._R8/   &  
     & 2._R8)))
 
!  collisional modes added with xnu_col
!  must fix for atomic mass dependence other than deuterium
      xnu_col=xparam(7)
      amass_e=2.7E-4_R8*(2._R8/amassgas)
 
! check adiabatic property that chkf should be 1.0  (finite xnu)
 
      chkf=(fnn*fpp-fpn*fnp)/((fnf*fpp-fpf*fnp)+epsilon)
      if (lprint.eq.2) write(6,*) 'chkf:', chkf
 
      if(neq.le.3) reps=0._R8
 
      f0=-idelta+(1._R8-dil)*apwt*(1/taui)*(1._R8-g1)                    &  
     &           +zimp**2*aiwt*(1/taui)*(1._R8-g1i)
      f1=1._R8-reps + f0
 
      xadiabat=0._R8
      if(nroot.le.6) then
        betae=0._R8
        f0=f1
        xadiabat=1._R8
      endif
      if(nroot.le.5) k_par=0._R8
 
      betae0=betae+epsilon
      k_par0=k_par+epsilon
 
      if (lprint.eq.98) then
        write(2,*) 'ky,g1,g2,g3,g12,g23,w_s,w_d,w_cd: ',                 &  
     &    ky,g1,g2,g12,g23,w_s,w_d,w_cd
        write(2,*) 'f0,taui,k_par,reps: ',                               &  
     &    f0,taui,k_par,k_per,reps
        write(2,*) 'chi_par_1,chi_per_1,gam_par,gam_per:',               &  
     &    chi_par_1,chi_per_1,gam_par,gam_per
        write(2,*) 'x_par,x_per,xmu:',                                   &  
     &    x_par,x_per,xmu
        write(2,*) 'xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per:',   &  
     &    xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per
        write(2,*) 'lamda_d,betae,xadiabat:',                            &  
     &    lamda_d,betae,xadiabat
        write(2,*) 'yt_mhd,xt_mhd:',                                     &  
     &    yt_mhd,xt_mhd
 
      endif
!
! matrix in order
! note ph=(n_i-n_t-n_u)/f0 results in (i,1)-(i,4)-(i,7) parts
!
! n_i equ #1
 
      amat(1,1)= (1._R8-dil)*apwt*                                       &  
     &  (-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,2)=                                                         &  
     &  +xi*w_d*taui*0.5_R8
 
      amat(1,3)=                                                         &  
     &  +xi*w_d*taui*0.5_R8
 
      amat(1,4)= -(-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,5)= 0._R8
 
      amat(1,6)=                                                         &  
     &  -xi*k_par
 
      amat(1,7)= -(-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,8)= 0._R8
 
      amat(1,9)= aiwt*zimp*                                              &  
     &  (-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,10)=0._R8
 
      amat(1,11)=0._R8
 
      amat(1,12)=0._R8
 
! p_par equ #2
 
      amat(2,1)= (1._R8-dil)*apwt*                                       &  
     &  (-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0                 &  
     &  +k_par*chi_par_1                                                 &  
     &  -(xi*w_d*taui*3._R8/2._R8-w_d*taui*xnu_par_par)                  &  
     &  -(xi*w_d*taui*1._R8/2._R8-w_d*taui*xnu_par_per)
 
      amat(2,2)=                                                         &  
     &  -k_par*chi_par_1                                                 &  
     &  +xi*w_d*taui*x_par +                                             &  
     &  (xi*w_d*taui*3._R8/2._R8-w_d*taui*xnu_par_par)
 
      amat(2,3)=                                                         &  
     &  (xi*w_d*taui*1._R8/2._R8-w_d*taui*xnu_par_per)
 
      amat(2,4)= -(-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,5)=0._R8
 
      amat(2,6)=                                                         &  
     &  -xi*gam_par*k_par
 
      amat(2,7)= -(-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,8)=0._R8
 
      amat(2,9)= aiwt*zimp*                                              &  
     &  (-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,10)=0._R8
 
      amat(2,11)=0._R8
 
      amat(2,12)=0._R8
 
! p_per equ #3
 
      amat(3,1)= (1._R8-dil)*apwt*                                       &  
     &  (-xi*w_s*((rlni-rlti)*g2+2._R8*rlti*g3)+xi*x_per*w_cd*g23)/f0    &  
     &  +k_par*chi_per_1                                                 &  
     &  -(xi*w_d*taui-w_d*taui*xnu_per_per)                              &  
     &  -(xi*w_d*taui*1._R8/2._R8-w_d*taui*xnu_per_par)
 
 
      amat(3,2)=                                                         &  
     &  +(xi*w_d*taui*1/2-w_d*taui*xnu_per_par)
 
      amat(3,3)=                                                         &  
     &  -k_par*chi_per_1                                                 &  
     &  +xi*w_d*taui*x_per   +(xi*w_d*taui-w_d*taui*xnu_per_per)
 
      amat(3,4)=                                                         &  
     & -(-xi*w_s*((rlni-rlti)*g2+2._R8*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,5)=0._R8
 
      amat(3,6)=                                                         &  
     &  -xi*gam_per*k_par
 
      amat(3,7)=                                                         &  
     &  -(-xi*w_s*((rlni-rlti)*g2+2._R8*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,8)=0._R8
 
      amat(3,9)= aiwt*zimp*                                              &  
     &   (-xi*w_s*((rlni-rlti)*g2+2._R8*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,10)=0._R8
 
      amat(3,11)=0._R8
 
      amat(3,12)=0._R8
 
! n_t equ #4
 
      amat(4,1)=(1._R8-dil)*apwt*                                        &  
     &  (-xi*w_s*rlne*reps*g0+xi*x43*3._R8/4*w_cd*reps*g0)/f0            &  
     &  -(1._R8-dil)*apwt*(-reps*fnf*(1._R8-reps)*g0/f0*xadiabat)
 
      amat(4,2)=0._R8
 
      amat(4,3)=0._R8
 
      amat(4,4)=-(-xi*w_s*rlne*reps*g0+xi*x43*3._R8/4*w_cd*reps*g0)/f0   &  
     &  -((1._R8-reps)*fnn)                                              &  
     &  -(-(-reps*fnf*(1._R8-reps)*g0/f0*xadiabat))
 
      amat(4,5)=                                                         &  
     &  -xi*w_d*x43*3._R8/4._R8                                          &  
     &  -((1._R8-reps)*fnp)
 
      amat(4,6)=0._R8
 
      amat(4,7)=-(-xi*w_s*rlne*reps*g0+xi*x43*3._R8/4*w_cd*reps*g0)/f0   &  
     &   -(-reps*fnf)
 
      amat(4,8)=0._R8
 
      amat(4,9)=aiwt*zimp*                                               &  
     &   (-xi*w_s*rlne*reps*g0+xi*x43*3._R8/4*w_cd*reps*g0)/f0           &  
     &   -aiwt*zimp*(-reps*fnf*(1._R8-reps)*g0/f0*xadiabat)
 
      amat(4,10)=0._R8
 
      amat(4,11)=0._R8
 
      amat(4,12)=0._R8
 
! p_t equ #5
 
      amat(5,1)= (1._R8-dil)*apwt*                                       &  
     &  (-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5._R8/4*w_cd*reps*g0)/f0     &  
     &  -(1._R8-dil)*apwt*(-reps*fpf*(1._R8-reps)*g0/f0*xadiabat)
 
      amat(5,2)=0._R8
 
      amat(5,3)=0._R8
 
      amat(5,4)=                                                         &  
     &  -(-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5._R8/4*w_cd*reps*g0)/f0    &  
     &            +xi*w_d*lamda_d                                        &  
     &  -((1._R8-reps)*fpn)                                              &  
     &  -(-(-reps*fpf*(1._R8-reps)*g0/f0*xadiabat))
 
      amat(5,5)=                                                         &  
     &  -xi*w_d*x43*5._R8/4._R8-xi*w_d*lamda_d                           &  
     &  -((1._R8-reps)*fpp)
 
      amat(5,6)=0._R8
 
      amat(5,7)=                                                         &  
     &  -(-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5._R8/4*w_cd*reps*g0)/f0    &  
     &  -(-reps*fpf)
 
      amat(5,8)=0._R8
 
      amat(5,9)= aiwt*zimp*                                              &  
     &  (-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5._R8/4*w_cd*reps*g0)/f0     &  
     &  -aiwt*zimp*(-reps*fpf*(1._R8-reps)*g0/f0*xadiabat)
 
      amat(5,10)=0._R8
 
      amat(5,11)=0._R8
 
      amat(5,12)=0._R8
 
! u_par equ #6
 
      amat(6,1)=(1._R8-dil)*apwt*                                        &  
     &   (-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0         &  
     &   -(betae/2._R8)*(-xi*k_par*(2._R8/betae0)*g0)/f0)
 
      amat(6,2)=-xi*k_par*taui
 
      amat(6,3)=0._R8
 
      amat(6,4)=                                                         &  
     &  -(-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0)        &  
     &  -(-(betae/2._R8)*(-xi*k_par*(2._R8/betae0)*g0)/f0)
 
      amat(6,5)=0._R8
 
      amat(6,6)=                                                         &  
     & +xi*w_d*(gam_par+gam_per)/2._R8-w_d*xmu
 
      amat(6,7)=                                                         &  
     &  -(-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0)        &  
     &  -(-(betae/2._R8)*(-xi*k_par*(2._R8/betae0)*g0)/f0)               &  
     &  -(betae/2._R8)*xi*k_par*(2._R8/betae0)/(1._R8-reps)
 
      amat(6,8)=                                                         &  
     &  -(betae/2._R8)*(-xi*w_s*(rlni*g1+rlti*g2))                       &  
     &  -(betae/2._R8)*(-xi*w_s*rlne)                                    &  
     &  +amass_e*xnu*xnu_col*k_per**2                                    &  
     &  -(betae/2._R8)*(-2._R8/betae0*amass_e*xnu*xnu_col*k_per**2)
! note there is no double counting in last two terms
 
      amat(6,9)=aiwt*zimp*                                               &  
     &  (-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0          &  
     &  -(betae/2._R8)*(-xi*k_par*(2._R8/betae0)*g0)/f0)
 
      amat(6,10)=0._R8
 
      amat(6,11)=0._R8
 
      amat(6,12)=0._R8
 
! n_u equ #7
 
      amat(7,1)=(1._R8-dil)*apwt*                                        &  
     &  (-xi*w_s*rlne*(1._R8-reps)*g0+xi*w_cd*                           &  
     &  (1._R8-x43*(3._R8/4._R8)*reps)*g0)/f0                            &  
     &  +(1._R8-dil)*apwt*(-reps*fnf*(1._R8-reps)*g0/f0*xadiabat)
 
      amat(7,2)=0._R8
 
      amat(7,3)=0._R8
 
      amat(7,4)=                                                         &  
     &  -(-xi*w_s*rlne*(1._R8-reps)*g0+xi*w_cd*                          &  
     &  (1._R8-x43*(3._R8/4._R8)*reps)*g0)/f0                            &  
     &  +((1._R8-reps)*fnn)                                              &  
     &  +(-(-reps*fnf*(1._R8-reps)*g0/f0*xadiabat))
 
      amat(7,5)=0._R8                                                    &  
     &  +((1._R8-reps)*fnp)
 
      amat(7,6)=-xi*k_par
 
      amat(7,7)=                                                         &  
     &  -(-xi*w_s*rlne*(1._R8-reps)*g0+xi*w_cd*                          &  
     &  (1._R8-x43*(3._R8/4._R8)*reps)*g0)/f0                            &  
     &  -xi*w_d*xt_mhd                                                   &  
     &  +(-reps*fnf)
 
      amat(7,8)=                                                         &  
     &  -xi*k_par*(-k_per**2)                                            &  
     &  -xi*w_d*yt_mhd*(w_s*(betae/2._R8)/k_par0*rlte)
 
      amat(7,9)=aiwt*zimp*                                               &  
     &  (-xi*w_s*rlne*(1._R8-reps)*g0+xi*w_cd*                           &  
     &  (1._R8-x43*(3._R8/4._R8)*reps)*g0)/f0                            &  
     &  +aiwt*zimp*(-reps*fnf*(1._R8-reps)*g0/f0*xadiabat)
 
      amat(7,10)=0._R8
 
      amat(7,11)=0._R8
 
      amat(7,12)=0._R8
 
! a_par equ #8
 
      amat(8,1)=(1._R8-dil)*apwt*(-xi*k_par*(2._R8/betae0)*g0/f0)
 
      amat(8,2)=0._R8
 
      amat(8,3)=0._R8
 
      amat(8,4)=-(-xi*k_par*(2._R8/betae0)*g0/f0)
 
      amat(8,5)=0._R8
 
      amat(8,6)=0._R8
 
      amat(8,7)=-(-xi*k_par*(2._R8/betae0)*g0/f0)                        &  
     &  +xi*k_par*(2._R8/betae0)/(1._R8-reps)
 
      amat(8,8)=-xi*w_s*rlne                                             &  
     &  -(2._R8/betae0)*amass_e*xnu*xnu_col*(k_per**2)
 
      amat(8,9)=aiwt*zimp*(-xi*k_par*(2._R8/betae0)*g0/f0)
 
      amat(8,10)=0._R8
 
      amat(8,11)=0._R8
 
      amat(8,12)=0._R8
 
! n_im equ #9
 
      amat(9,1)= (1._R8-dil)*apwt*                                       &  
     &  (-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,2)=0._R8
 
      amat(9,3)=0._R8
 
      amat(9,4)= -(-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/   &  
     & f0
 
      amat(9,5)= 0._R8
 
      amat(9,6)= 0._R8
 
      amat(9,7)= -(-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/   &  
     & f0
 
      amat(9,8)= 0._R8
 
      amat(9,9) = aiwt*zimp*                                             &  
     &  (-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,10)=                                                        &  
     &  +xi*w_d*taui*0.5_R8/zimp
 
      amat(9,11)=                                                        &  
     &  +xi*w_d*taui*0.5_R8/zimp
 
      amat(9,12)=                                                        &  
     &  -xi*k_par
 
! pim_par equ #10
 
      amat(10,1)= (1._R8-dil)*apwt*                                      &  
     &  (-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
 
      amat(10,2)=0._R8
 
      amat(10,3)=0._R8
 
      amat(10,4)= -(-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/   &  
     & f0
 
      amat(10,5)=0._R8
 
      amat(10,6)=0._R8
 
      amat(10,7)= -(-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/   &  
     & f0
 
      amat(10,8)=0._R8
 
      amat(10,9)= aiwt*zimp*                                             &  
     &   (-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0           &  
     &   +k_par*chi_par_1/sqrt(mimp)                                     &  
     &   -(xi*w_d*taui*3._R8/2._R8-w_d*taui*xnu_par_par)/zimp            &  
     &   -(xi*w_d*taui*1._R8/2._R8-w_d*taui*xnu_par_per)/zimp
 
      amat(10,10)=                                                       &  
     &  -k_par*chi_par_1/sqrt(mimp)                                      &  
     &  +xi*w_d*taui*x_par/zimp +(xi*w_d*taui*3._R8/2._R8                &  
     &  -w_d*taui*xnu_par_par)/zimp
 
      amat(10,11)=                                                       &  
     &  (xi*w_d*taui*1._R8/2._R8-w_d*taui*xnu_par_per)/zimp
 
      amat(10,12)=                                                       &  
     &  -xi*gam_par*k_par
 
! pim_per equ #11
 
      amat(11,1)= (1._R8-dil)*apwt*                                      &  
     &  (-xi*w_s*((rlnimp-rlti)*g2i+2._R8*rlti*g3i)+                     &  
     &  xi*x_per*w_cd*g23i)/f0
 
      amat(11,2)= 0._R8
 
      amat(11,3)= 0._R8
 
      amat(11,4)=                                                        &  
     &  -(-xi*w_s*((rlnimp-rlti)*g2i+2._R8*rlti*g3i)+                    &  
     &  xi*x_per*w_cd*g23i)/f0
 
      amat(11,5)=0._R8
 
      amat(11,6)=0._R8
 
      amat(11,7)=                                                        &  
     &  -(-xi*w_s*((rlnimp-rlti)*g2i+2._R8*rlti*g3i)+                    &  
     &  xi*x_per*w_cd*g23i)/f0
 
      amat(11,8)=0._R8
 
      amat(11,9)= aiwt*zimp*                                             &  
     &  (-xi*w_s*((rlnimp-rlti)*g2i+2._R8*rlti*g3i)+                     &  
     &  xi*x_per*w_cd*g23i)/f0                                           &  
     &  +k_par*chi_per_1/sqrt(mimp)                                      &  
     &  -(xi*w_d*taui-w_d*taui*xnu_per_per)/zimp                         &  
     &  -(xi*w_d*taui*1._R8/2._R8-w_d*taui*xnu_per_par)/zimp
 
      amat(11,10)=                                                       &  
     &  +(xi*w_d*taui*1/2-w_d*taui*xnu_per_par)/zimp
 
      amat(11,11)=                                                       &  
     &  -k_par*chi_per_1/sqrt(mimp)                                      &  
     &  +xi*w_d*taui*x_per/zimp                                          &  
     &  +(xi*w_d*taui-w_d*taui*xnu_per_per)/zimp
 
      amat(11,12)=                                                       &  
     &  -xi*gam_per*k_par
 
! uim_par equ #12
!gms 5/21/99 added gamma_p to amat(12,1),amat(12,4),amat(12,7),amat(12,9)
!    added xnu_col term to amat(12,8)
!    fixed mimp factor in amat(12,12)
 
      amat(12,1)=(1._R8/mimp)*(1._R8-dil)*apwt*                          &  
     &  ((-xi*k_par*g1i/f0)*zimp                                         &  
     &  -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0                         &  
     &  -(betae/2._R8)*zimp*(-xi*k_par*(2._R8/betae0)*g0i)/f0)
 
      amat(12,2)=0._R8
 
      amat(12,3)=0._R8
 
      amat(12,4)=                                                        &  
     & -(1._R8/mimp)*(-xi*k_par*g1i/f0)*zimp                             &  
     & -(1._R8/mimp)*(-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0)          &  
     & -(1._R8/mimp)*(-(betae/2._R8)*(-xi*k_par                          &  
     & *(2._R8/betae0)*g0i)/f0)*zimp
 
      amat(12,5)=0._R8
 
      amat(12,6)=0._R8
 
       amat(12,7)=                                                       &  
     & -(1._R8/mimp)*(-xi*k_par*g1i/f0)*zimp                             &  
     & -(1._R8/mimp)*(-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0)          &  
     & -(1._R8/mimp)*(-(betae/2._R8)*                                    &  
     &  (-xi*k_par*(2._R8/betae0)*g0i)/f0)*zimp                          &  
     & -(1._R8/mimp)*(betae/2._R8)*xi*                                   &  
     &  k_par*(2._R8/betae0)/(1._R8-reps)*zimp
 
      amat(12,8)=                                                        &  
     &  -(1._R8/mimp)*(betae/2._R8)*(-xi*w_s*(rlnimp*g1i+rlti*g2i))      &  
     &  -(1._R8/mimp)*(betae/2._R8)*(-xi*w_s*rlne)*zimp                  &  
     &  +(1._R8/mimp)*zimp*amass_e*xnu*xnu_col*k_per**2                  &  
     &  -(1._R8/mimp)*(betae/2._R8)*                                     &  
     &   (-2._R8/betae0*amass_e*xnu*xnu_col*k_per**2)*zimp
 
      amat(12,9)=(1._R8/mimp)*aiwt*zimp*                                 &  
     &   ((-xi*k_par*g1i/f0)*zimp                                        &  
     &  -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0                         &  
     &  -(betae/2._R8)*zimp*(-xi*k_par*(2._R8/betae0)*g0i)/f0)
 
      amat(12,10)=-(1._R8/mimp)*xi*k_par*taui
 
      amat(12,11)=0._R8
 
      amat(12,12)= (1._R8/mimp)*(xi*w_d*(gam_par+gam_per)                &  
     &  /2._R8/zimp -w_d*xmu/zimp)
!
! put in rot shear stabilization and possible source of gyrobohm breaking
! and model damping kdamp
!
!***********************************************************************
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
! solve 12x12 complex
! -xi*omega*v(i)=sum_j amat(i,j)*v(j)  omega=freq+xi*gamma
! upto nroot
! order with max gamma and find eigenvector v(i) with ant fixed norm.
!
!...Fill matricies for eigenvalue equation
!
      do j1=1,neq
        rr(j1) = 0.0_R8
        ri(j1) = 0.0_R8
        do j2=1,neq
          ar(j1,j2) = REAL(  amat(j1,j2) )
          ai(j1,j2) = AIMAG( amat(j1,j2) )
!...test tmp
!         ai(j1,j2) = 0.0
!         ar(j1,j2) = 0.0
!         if (j1.eq.j2) ar(j1,j2)=j1
!
          vr(j1,j2) = 0.0_R8
          vi(j1,j2) = 0.0_R8
        enddo
      enddo
!
!...diagnostic output
!
      if ( lprint .gt. 6 ) then
        write (1,*)
        write (1,*) ' ar(j1,j2)  j2 ->'
        do j1=1,neq
          write (1,192) (ar(j1,j2),j2=1,neq)
        enddo
!
        write (1,*)
        write (1,*) ' ai(j1,j2)  j2->'
        do j1=1,neq
          write (1,192) (ai(j1,j2),j2=1,neq)
        enddo
 192    format (1p8e10.2)
 193    format (1p8e12.4)
      endif
!
!..find the eigenvalues and eigenvectors
!
!.. eigen_gf = 0 use cgg solver (default)
!..          = 1 use tomsqz solver
!..          = 2 use zgeev solver
!.. not longer used:
!        call f02ake( ar,iar,ai,iai,ieq,rr,ri,vr,ivr,vi,ivi,
!     >               intger, ifail )
!
        ifail = 0
!
        if (eigen_gf .eq. 2 ) then
!
        jobvl = 'N'
        jobvr = 'V'
        do j1=1,neq
         do j2=1,ieq
           mata(j1,j2) = CMPLX(ar(j1,j2),ai(j1,j2),R8)
         enddo
        enddo
!
!        call zgeev(jobvl,jobvr,ieq,mata,neq,w,cvl,neq,cvr,
!     &             neq,work,lwork,rwork,ifail)
        do j1=1,neq
         rr(j1) = REAL(w(j1))
         ri(j1) = AIMAG(w(j1))
         do j2=1,ieq
           vr(j1,j2) = REAL(cvr(j1,j2))
           vi(j1,j2) = AIMAG(cvr(j1,j2))
         enddo
        enddo
!
        elseif (eigen_gf .eq. 1 ) then
!
        do j2=1,neq
           do j1=1,neq
              bi(j1,j2)=0.0_R8
              if(j1.eq.j2) then
                 br(j1,j2)=1.0_R8
              else
                 br(j1,j2)=0.0_R8
              endif
           enddo
        enddo
!
        call r8tomsqz(neq,ieq,ar,ai,br,bi, rr,ri,beta_tom, vr,vi, ifail)    
!
        do j1=1,ieq
           ztemp1 = beta_tom(j1)
           if ( abs(beta_tom(j1)) .lt. epsilon ) ztemp1 = epsilon
           rr(j1)=rr(j1) / ztemp1
           ri(j1)=ri(j1) / ztemp1
        enddo
!
        else
!
        matz=1
!       write(*,*) 'neq = ',neq
!       write(*,*) 'ieq = ',ieq
!       write(*,*) 'matz = ',matz
!       write (*,*) ' ar(j1,j2)  j2 ->'
!       do j1=1,neq
!         write (*,193) (ar(j1,j2),j2=1,neq)
!       enddo
!       write (*,*) ' ai(j1,j2)  j2 ->'
!       do j1=1,neq
!         write (*,193) (ai(j1,j2),j2=1,neq)
!       enddo
!
        call cgg_glf(neq,ieq,ar,ai,rr,ri,matz,vr,vi,fv1,fv2,fv3,ifail)
!
!       write (*,*) ' wr(j1) and wi(j1)'
!       do j1=1,neq
!         write (*,193) rr(j1), ri(j1)
!       enddo
!       write (*,*) ' zr(j1,j2)  j2 ->'
!       do j1=1,neq
!         write (*,193) (vr(j1,j2),j2=1,neq)
!       enddo
!       write (*,*) ' zi(j1,j2)  j2 ->'
!       do j1=1,neq
!         write (*,193) (vi(j1,j2),j2=1,neq)
!       enddo
 
        endif
!
        if ( lprint .gt. 1 ) then
          write (1,*) ifail,' = ifail routine '
        endif
!
!..print eigenvalues
!
        if ( lprint .gt. 6 ) then
          write (1,121)
          do j=1,ieq
            write (1,122) rr(j), ri(j)
          enddo
 121      format (/' Solution of the eigenvalue equations'               &  
     &     /t4,'real   ',t18,'imag   ')
 122      format (1p2e14.5)
        endif
!
!...Store the complex eigenvectors and eigenvalues
!...Note the routines here solve A.v = lambda v
!...but that the equation solved is A.v = -i omega v
!...The i-th column of v is the i-th eigenvector
!
        do j1=1,ieq
          zomega(j1) = xi*(rr(j1)+xi*ri(j1))
          do j2=1,ieq
            zevec(j2,j1) = vr(j2,j1) + xi*vi(j2,j1)
          enddo
        enddo
!
        if ( lprint .gt. 6 ) then
          write (6,123)
          do j=1,ieq
            write (6,122) REAL(zomega(j)), AIMAG(zomega(j))
          enddo
 123      format (/' Multiplied by i: '                                  &  
     &     /t4,'zomegar',t18,'zomegai')
        endif
        do iroot=1,4
!
!
!..save growth rates and frequencies in real variables
!
        zgamax = 0.0_R8
!temp
        zgamax = -1.E10_R8
        jmax=0
        gamma=0._R8
        do j=1,ieq
         if(j.ne.jroot(1).and.j.ne.jroot(2).and.j.ne.jroot(3)) then
          if (AIMAG(zomega(j)).gt. zgamax) then
            zgamax = AIMAG(zomega(j))
            jmax=j
          endif
         endif
        enddo
!
        freq = REAL( zomega(jmax) )
        gamma = AIMAG( zomega(jmax) )
!
! skip stable modes
!        if(gamma.lt.0.D0)go to 775
 
         jroot(iroot)=jmax
 
!temp        if(zgamax.lt.zepsqrt) gamma=0.
 
        if (jmax.ne.0) then
         gammaroot(iroot)=gamma
         freqroot(iroot)=freq
 
 
        do j=1,12
         v(j)=0._R8
        enddo
        do j=1,ieq
          v(j) = zevec(j,jmax)
        enddo
!
!***********************************************************************
 
      n_i=0._R8
      p_par=0._R8
      p_per=0._R8
      n_t=0._R8
      p_t=0._R8
      u_par=0._R8
      n_u=0._R8
      a_par=0._R8
      n_im=0._R8
      p_im_par=0._R8
      p_im_per=0._R8
!     u_im_par=0.
 
      t_u=0._R8
 
 
      n_i=v(1)
      p_par=v(2)
      p_per=v(3)
      if(ieq.ge.5) n_t=v(4)
      if(ieq.ge.5) p_t=v(5)
      if(ieq.ge.6) u_par=v(6)
      if(ieq.ge.8) n_u=v(7)
      if(ieq.ge.8) a_par=v(8)
      if(ieq.ge.9) n_im=v(9)
      if(ieq.eq.12) p_im_par=v(10)
      if(ieq.eq.12) p_im_per=v(11)
!     if(ieq.eq.12) u_im_par=v(12)
 
      if (ieq.ge.8) ph=((1._R8-dil)*apwt*n_i+aiwt*zimp*n_im-n_t-n_u)/f0
      if (ieq.lt.8) ph= ((1._R8-dil)*apwt*n_i-n_t)/f0
      if (ieq.le.3) ph= ((1._R8-dil)*apwt*n_i)/f0
      t_u=(betae/2._R8)*w_s/k_par0*rlte*a_par
 
      n_e=(1._R8-dil)*apwt*(n_i-(g0-g1)/taui*ph)                         &  
     &     +aiwt*zimp*(n_im-zimp*(g0i-g1i)/taui*ph)
 
! impurity trace convective limit
      if (aiwt.lt.-epsilon) then
       n_im=0._R8
       do j=1,8
        n_im=n_im+amat(9,j)*v(j)/(-xi*freq+gamma)
       enddo
      endif
 
! idelta=xi*yparam(1)+yparam(2)   for trapped electrons
 
      yparam(1)=AIMAG(-(n_t-reps*ph)/ph)
      yparam(2)=REAL(-(n_t-reps*ph)/ph)
 
 
      chknu=n_u/(1._R8-reps)/ph
      chknt=n_t/reps/ph
      chknt2=n_t*(f0+reps)/(reps*((1._R8-dil)*apwt*n_i+aiwt*zimp*n_im))
      if (lprint.eq.2) write (6,*) 'chknu,chknt,chknt2:',                &  
     &    chknu,chknt,chknt2
 
! non linear saturation rule
 
      gamma_r= 0.2_R8*3._R8/2._R8*abs(w_d)*taui   !only scaling important
      if(iglf.eq.1) gamma_r= 0.2_R8*3._R8/2._R8*abs(w_d0)*taui
!
      gamma_net=gamma-abs(alpha_star*gamma_star                          &  
     &         +alpha_e*gamma_e+alpha_mode*gamma_mode)-kdamp
 
      gamma_net=max(gamma_net,xparam(8)*gamma)
      if( gamma_net.gt.0._R8)then
       ph_m=gamma_net**(1._R8-adamp)*gamma_r**adamp/(k_m*ky)
! set flag ngrow_k_gf: found at least one unstable mode for this k
       if(ipert_gf.eq.0)ngrow_k_gf(iky0)=1
      endif
 
      if( gamma_net.le.0._R8) ph_m=0._R8
!
      if(xparam(24).gt.0._R8) then
        if(gamma.gt.0._R8)then
          ph_m=abs(gamma)**(1._R8-adamp)*gamma_r**adamp/(k_m*ky)/        &  
     &    sqrt(1._R8+(abs(alpha_star*gamma_star+                         &  
     &    alpha_e*gamma_e+alpha_mode*gamma_mode)/                        &  
     &    (abs(gamma)+.00001_R8))**xparam(24))
          if(ipert_gf.eq.0)ngrow_k_gf(iky0)=1
        else
           ph_m=0._R8
        endif
      endif
 
! 7.17.96
      if(xparam(22).gt.0) then
        if(gamma.gt.0._R8) ph_m=gamma**(1._R8-adamp-xparam(22))          &  
     &   *gamma_r**adamp/(k_m*ky)
        if(gamma.le.0._R8) ph_m=0._R8
      endif
 
         phi_norm=0
         phi_normroot(iroot)=0._R8
 
       if( ph_m.gt.0._R8) then
 
       phi_norm=ph_m*ph_m/ABS((conjg(ph)*ph))
 
       phi_normroot(iroot)=phi_norm
 
! note only real part survives in diffusivities
!    ...units are c_s*rho_s**2/a
! magnetic futter component is too small to worry about
 
      d_hat    = phi_norm*REAL(conjg(n_i)*(-xi*ky*ph))/rlni             &  
     &+d_hat
 
      d_im_hat    = phi_norm*REAL(conjg(n_im)*(-xi*ky*ph))/(rlnimp+     &  
     &epsilon)                                                           &  
     &+d_im_hat
 
      chii_hat = phi_norm*3._R8/2._R8*                                   &  
     &REAL(conjg((1._R8/3._R8)*p_par+(2._R8/3._R8)*p_per)*(-xi*ky*ph))/  &  
     & rlti                                                              &  
     &+chii_hat
      chii_hat=chii_hat + aiwt/apwt*xparam(21)*phi_norm*3._R8/2._R8*     &  
     &REAL(conjg((1._R8/3._R8)*p_im_par+(2._R8/3._R8)*p_im_per)*        &  
     &   (-xi*ky*ph))/rlti
 
      chie_hat = phi_norm*3._R8/2._R8*                                   &  
     &REAL(conjg(p_t+n_u+t_u)*(-xi*ky*ph))/rlte                         &  
     &+chie_hat
 
! electron to ion energy exchange in units n0*t0*c_s/a*(rho_a/a)**2
! note here we interpret QLT d/dt=-xi*freq dropping gamma part to
! avoid getting nonzero result for n_e->ph adiabatic limit
! ie <(d n_e/dt)*conjg(ph)>_time ave -> 0 adiabatic limit
! note  (-1) means exch_hat is electron to ion rate or
! ion heating rate, ie positive exch_hat cools electrons and heats ions
 
      exch_hat = phi_norm*(-1._R8)*                                      &  
     &REAL(conjg(-xi*freq*n_e)*ph)                                      &  
     &+exch_hat
 
      eta_par_hat=(1._R8-xparam(14))*phi_norm*                           &  
     &REAL(conjg(u_par)                                                 &  
     &*(-xi*ky*ph))/(gamma_p+epsilon)*(-gamma_p*alpha_n)                 &  
     &+xparam(14)*phi_norm*REAL(conjg(                                  &  
     &-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1*ph/(-xi*freq+gamma))          &  
     &*(-xi*ky*ph))/(gamma_p+epsilon)*(-gamma_p*alpha_n)                 &  
     &+eta_par_hat
 
      eta_per_hat = phi_norm*                                            &  
     &REAL(conjg(-ky*(ky*shat*rms_theta)*ph)*                           &  
     &   (ph+taui*((1._R8/3._R8)*p_par+(2._R8/3._R8)*p_per)))*           &  
     &   (-gamma_p*alpha_n)                                              &  
     &/(gamma_p+epsilon)                                                 &  
     &+eta_per_hat
 
       endif
       endif
      enddo
 777  continue
!
      if(ilh.eq.1) then
 
       do i=1,nmode
        yparam_k_gf(i,iky)=yparam(i)
       enddo
       xkyf_k_gf(iky)=kyf
 
       do iroot=1,4
        gamma_k_gf(iroot,iky)=gammaroot(iroot)
        freq_k_gf(iroot,iky)=freqroot(iroot)
        phi_norm_k_gf(iroot,iky)=phi_normroot(iroot)
       enddo
 
       diff_k_gf(iky)=d_hat
       diff_im_k_gf(iky)=d_im_hat
       chii_k_gf(iky)=chii_hat
       chie_k_gf(iky)=chie_hat
       exch_k_gf(iky)=exch_hat
 
 
!test       exch_gf=-freq_gf(1)/xkyf_gf*diff_gf*rlni
 
       eta_par_k_gf(iky)=eta_par_hat
       eta_per_k_gf(iky)=eta_per_hat
 
! b_pol/b_phi=rmin/(rmaj*q)
       eta_phi_k_gf(iky)=eta_par_k_gf(iky)+                              &  
     &  rmin/(rmaj*q)*eta_per_k_gf(iky)
 
      endif
!
! computed high-k ETG electron transport at each ky
! Note: not added to ITG chi-e here ... done after ky loop
      chie_e_k_gf(iky)=0._R8
      if(ilh.eq.2) then
        chie_e_k_gf(iky)=xparam(10)*chii_hat*                            &  
     &                   taui_gf**(3._R8/2._R8)/                         &  
     &                   (1836._R8*amassgas_gf)**.5_R8
      endif
! end ky loop
      enddo
 
! end big loop on ilh ... no longer used
!     enddo
!
! check to see if any unstable modes were found
!
      if(ipert_gf.eq.0)then
        do j=1,nmode
         if(ngrow_k_gf(j).ne.0)ngrow_k_gf(0)=1
        enddo
        if(ngrow_k_gf(0).eq.0)go to 888
      endif
!
!
!***********************************************************************
! initializations for summations over ky
!
      anorm_k=0._R8
      diff_gf=0._R8
      diff_im_gf=0._R8
      chii_gf=0._R8
      chie_gf=0._R8
      exch_gf=0._R8
      eta_par_gf=0._R8
      eta_per_gf=0._R8
      eta_phi_gf=0._R8
      chie_e_gf=0._R8
!
! Sum ITG and ETG transport
! over logarithmic ky grid (d ky=ky*d yk)
!
      do iky=1,ikymax_gf
       del_k=xkyf_k_gf(iky)
       anorm_k=anorm_k+del_k
       diff_gf=diff_gf+diff_k_gf(iky)*del_k
       diff_im_gf=diff_im_gf+diff_im_k_gf(iky)*del_k
       chii_gf=chii_gf+chii_k_gf(iky)*del_k
       chie_gf=chie_gf+chie_k_gf(iky)*del_k
       exch_gf=exch_gf+exch_k_gf(iky)*del_k
       eta_par_gf=eta_par_gf+eta_par_k_gf(iky)*del_k
       eta_per_gf=eta_per_gf+eta_per_k_gf(iky)*del_k
       eta_phi_gf=eta_phi_gf+eta_phi_k_gf(iky)*del_k
       chie_e_gf=chie_e_gf+chie_e_k_gf(iky)*del_k
      enddo
!
! Add ITG and ETG electron transport
!
      chie_gf=chie_gf + 1._R8*chie_e_gf
!
      diff_gf=diff_gf/anorm_k
      diff_im_gf=diff_im_gf/anorm_k
      chii_gf=chii_gf/anorm_k
      chie_gf=chie_gf/anorm_k
      exch_gf=exch_gf/anorm_k
      eta_par_gf=eta_par_gf/anorm_k
      eta_per_gf=eta_per_gf/anorm_k
      eta_phi_gf=eta_phi_gf/anorm_k
      chie_e_gf=chie_e_gf/anorm_k
!
!
! pick off maximum gamma
!
      do iroot=1,4
       gamma_k_max=-1.E6_R8
       do iky=1,ikymax_gf
        if(gamma_k_gf(iroot,iky).gt.gamma_k_max) then
         gamma_k_max=gamma_k_gf(iroot,iky)
         gamma_gf(iroot)=gamma_k_gf(iroot,iky)
         freq_gf(iroot)=freq_k_gf(iroot,iky)
         xky_gf(iroot)=xkyf_k_gf(iky)
        endif
       enddo
      enddo
!
! pick off 2nd maximum gamma
!
       gamma_k_max=-1.E6_R8
       do iky=1,ikymax_gf
        if( (gamma_k_gf(1,iky).gt.gamma_k_max) .and.                     &  
     &      (gamma_k_gf(1,iky).lt.gamma_gf(1)) ) then
         gamma_k_max=gamma_k_gf(1,iky)
         gamma_gf(2)=gamma_k_gf(1,iky)
         freq_gf(2)=freq_k_gf(1,iky)
         xky_gf(2)=xkyf_k_gf(iky)
        endif
       enddo
!
!       write(6,*) gamma_gf(1), gamma_gf(2), xky_gf(1), xky_gf(2)
!
! print to file log
!      write(*,66)chii_gf,(gamma_gf(j),j=1,4)
 66    format(f14.9,4f14.9)
 67    format(2i2,f14.9)
 
 
      if(xparam(22).gt.0._R8) then
       phi_renorm=1._R8
       gamma_gross_net=gamma_gf(1)-abs(alpha_star*gamma_star             &  
     &         +alpha_e*gamma_e+alpha_mode*gamma_mode)-kdamp
       if(gamma_gross_net.gt.0._R8)                                      &  
     &  phi_renorm=gamma_gross_net**(xparam(22)*2._R8)
       if(gamma_gross_net.le.0._R8) phi_renorm=0._R8
 
      diff_gf=diff_gf*phi_renorm
      diff_im_gf=diff_im_gf*phi_renorm
      chii_gf=chii_gf*phi_renorm
      chie_gf=chie_gf*phi_renorm
      exch_gf=exch_gf*phi_renorm
      eta_par_gf=eta_par_gf*phi_renorm
      eta_per_gf=eta_per_gf*phi_renorm
      eta_phi_gf=eta_phi_gf*phi_renorm
      chie_e_gf=chie_e_gf*1.0_R8
 
 
      endif
 
! put in cnorm_gf 12/22/96
 
 
      diff_gf=cnorm_gf*diff_gf
      diff_im_gf=cnorm_gf*diff_im_gf
      chii_gf=cnorm_gf*chii_gf
      chie_gf=cnorm_gf*chie_gf
      exch_gf=cnorm_gf*exch_gf
      eta_par_gf=cnorm_gf*eta_par_gf
      eta_per_gf=cnorm_gf*eta_per_gf
      eta_phi_gf=cnorm_gf*eta_phi_gf
      chie_e_gf=cnorm_gf*chie_e_gf
 
 
      if(lprint.gt.0) then
          write(1,*) 'gamma_gf=',  gamma_gf
          write(1,*) 'freq_gf=',  freq_gf
          write(1,*) 'ph_m=',  ph_m
          write(1,*) 'diff_gf=', diff_gf
          write(1,*) 'diff_im_gf=', diff_im_gf
          write(1,*) 'chii_gf=', chii_gf
          write(1,*) 'chie_gf=', chie_gf
          write(1,*) 'exch_gf=', exch_gf
      endif
 
 
      if (lprint.eq.98) then
        write(6,*) 'rlti,rlte,rlne,rlni,rlnimp: ',                       &  
     &    rlti,rlte,rlne,rlni,rlnimp
       write(6,*) 'chii,chie,diff,diff_im: ',                            &  
     &    chii_gf,chie_gf,diff_gf,diff_im_gf
       write(6,*) 'gamma_gf,freq_gf,ph_m: ',                             &  
     &    gamma_gf,freq_gf,ph_m
       write(6,*) 'jmax: ',                                              &  
     &    jmax
        write(2,*) 'rlti,rlte,rlne,rlni,rlnimp: ',                       &  
     &    rlti,rlte,rlne,rlni,rlnimp
       write(2,*) 'chii,chie,diff,diff_im: ',                            &  
     &    chii_gf,chie_gf,diff_gf,diff_im_gf
       write(2,*) 'gamma_gf,freq_gf,ph_m: ',                             &  
     &    gamma_gf,freq_gf,ph_m
       write(2,*) 'jmax: ',                                              &  
     &    jmax
 
        write (2,*) ' ar(j1,j2)  j2 ->'
        do j1=1,neq
          write (2,*) (ar(j1,j2),j2=1,neq)
        enddo
!
        write (2,*)
        write (2,*) ' ai(j1,j2)  j2->'
        do j1=1,neq
          write (2,*) (ai(j1,j2),j2=1,neq)
        enddo
!
        write (2,*) ' vr(j1,j2)  j2 ->'
        do j1=1,neq
          write (2,*) (vr(j1,j2),j2=1,neq)
        enddo
!
        write (2,*)
        write (2,*) ' vi(j1,j2)  j2->'
        do j1=1,neq
          write (2,*) (vi(j1,j2),j2=1,neq)
        enddo
 
      endif
 
 999  continue
      if (lprint.gt.0) close(1)
 
      return
!
! return for case with  no unstable modes
!
 888  continue
      diff_gf=0._R8
      diff_im_gf=0._R8
      chii_gf=0._R8
      chie_gf=0._R8
      exch_gf=0._R8
      eta_par_gf=0._R8
      eta_per_gf=0._R8
      eta_phi_gf=0._R8
      chie_e_gf=0._R8
      do j1=1,4
        gamma_gf(j1)=0._R8
        freq_gf(j1)=0._R8
        xky_gf(j1)=0._R8
      enddo
      return
 
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
