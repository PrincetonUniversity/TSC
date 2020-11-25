c#include "f77_dcomplx.h"
 
c
cglf2d.f 12-mar-03 Kinsey
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
 
           subroutine glf2d(iglf)
 
c
c***********************************************************************
c questions  should be addressed to
c  Ron Waltz 619-455-4584  or email: waltz@gav.gat.com
c***********************************************************************
 
 
c 2D GLF equations with massless isothermal passing electrons from
c  Waltz et al, Phys. of Plasmas 6(1995)2408
 
c In Eq. 24  p_u_par is replaced by n_u/(1-reps) +t_u and
c the isothermal conditions t_u= (betae/2)*w_s/k_par*(rlt)*a_par is
c used. Thus n_u = (1-reps)*ph (adiabatic passing electrons) for betae=0
c
c In Eq. 23 (p_u_par+p_u_per) is replaced by n_u+t_u
c using the isothermal condition the mHD beta limit is too low by
c 1/(1+ reps(-0.25+0.75rlte)/(taui*(rln+rlti)+(rln+rlte))
c It is possible to patch this up by replacing with x*n_u+y*t_u
c then solving for x and y to obtain
c universal MHD betae_crit=2*k_par**2/(w_d*(taui*(rln+rlti)+w_d*(rln+rlti)))
c beta_crit=(1+taui)*betae_crit=(1/2)(1/q**2)*(L_p/rmajor)
c 1/2 is replaced by s_hat in models with shear
 
c EVERYTHING else including normalizing units follows the paper.
 
c  unit of microlength is rho_s, unit of macrolength is a
c   a is typically rho value at separatrix
c  unit of time is a/c_s; unit of diffusion is (c_s/a)*rho_s**2
c  c_s=sqrt(Te/M_i), omega=eB/cM_i, rho_s=c_s/omega
 
c example balance equations to clarify meaning  of diffusivities
c
c       chie_hat is effective energy diffusivity
c
c  (3/2) d ne Te/ dt =
c -1/V(rho)_prime d /d rho V(rho)_prime
c             |grad rho|**2 ne chie_hat (c_s rho_s**2/a)(-d Te/ d rho)
c      -exch_hat (ne Te) c_s/a (rho_s/a)**2 +heating density
c
c and similarly for ion equation
c note that no convective part is added, ie "convection" is included
c inside chie_hat
c note no impurity energy flow is computed although it could be easily done
 
c        d_hat is effective plasma diffusivity for ions
c
c        d ni / dt =
c -1/V(rho)_prime d /d rho V(rho)_prime
c               |grad rho|**2 ne d_hat (c_s rho_s**2/a) (-d ni/ d rho)
c        + plasma source density
 
c        d_im_hat is the effective plasma diffusivity for impurities
 
c        eta_phi is the toroidal viscosity or momentum diffusivity
c
c       M ni d v_phi/ dt =
c error found 5/12/98  should be d (M ni v_phi)/ dt =
c -1/V(rho)_prime d /d rho V(rho)_prime
c      |grad rho|**2 ( M ni eta_phi_hat (c_s rho_s**2/a) (-d v_phi/ d rho)
c                    +    M v_phi ne d_hat (c_s rho_s**2/a) (-d ne/ d rho))
c        + toroidal momentum source density
c
c note that a convective part was added
c
c  eta_par_hat and eta_per_hat are diagnostic. See CAUTION on eta_phi_hat
c  at large gamma_p=  (-d v_phi/ d rho) /(c_s/a)
c
c  chie_e_gf is the eta_e mode electron transport which is te <-> ti
c  and mi <-> me isomorphic to eta_i (ITG) ion transport
c  with adiabatic electrons.
c  these mode obtain at high-k where the ions are adiabatic from
c  the gyro cut-off.
c  their wave numbers are sqrt(mi/me) larger than ITG modes and
c  since their frequencies are sqrt(mi/me) larger, they are not
c  rotationally shaer satbilized.
c  when xparam_gf(10).eq.0 xparam_gf(10)*chie_e_gf is added to
c  chie_gf and chie_e_gf is a diagnostic.
 
c input

c  eigen_gf = 0 use cgg eigenvalue solver (default)
c           = 1 use generalized tomsqz eigenvalue solver
c           = 2 use zgeev eigenvalue solver
c  nroot number of equations
c  iflagin(1:20) control flags
c   iflagin(1) 0 use ky=ky0; 1 use landau damping point
c   iflagin(2) 0. local w_d and k_par "2d"; 1 fit to trial function "3d"
c   iflagin(3) 0,1,and 2 fix up park low high beta and beta lim elong factor
c   iflagin(4) 0 trapped electron Waltz EoS 1 weiland EoS
c   iflagin(5) rms_theta 0:fixed; 1 inverse to q/2 ; 2 inverse to root q/2
c                        3: inverse to xparam(13)*(q/2-1)+1.
c              5 for retuned rms-theta
c  xparam(1:20) control parameters
c   xparam(1:2): idelta=xi*xparam(1)+xparam(2) nonadiabatic electron response
c   xparam(3) multiplier park_gf(high betae)/ park_gf(low betae) -1
c   xparam(6)+1. is enhancement of xnueff
c   xparam(7) coef of resistivity
c   xparam(8) cut off on rotational stabilization
c   xparam(9)+1. is shape (triangularity) enhancement to beta_crit
c   xparam(10) is high k electron mode enhancement
c   xparam(11:12) lamda parameters
c   xparam(13) rms_theta q-dependence
c   xparam(14)  adjustment to gamma_p avoiding negative viscosity
c   xparam(15)   (1+xparam(15)*reps trapped electron fraction
c   xparam(16) rms_theta shat dependence
c   xparam(17) ""
c   xparam(18) rms_theta betae dependence
c   xparam(19:20)  extra
c   xparam(21) 1 add impurity energy diffusivity to ion energy diffusivity
c   xparam(22) >0 keeps gamma_e from changeing spectrum
c   xparam(23) 1. kills kx**2 in k_m**2
c   xparam(24) exb damping model
c  ky0=k_theta*rho_s; k_theta= nq/r; normally 0.3
c  rms_theta width of phi**2 mode function for best fit near pi/3
c  rlti=a/L_Ti   a/L_f= sqrt(kappa) a d ln f / d rho
c  rlte=a/L_Te
c  rlne= a/L_ne
c  rlni= a/L_ni
c  rlnimp= a/L_nim
c  dil=1.-ni_0/ne_0  dilution
c  apwt = ni_0/ne_0
c  aiwt = nim_0/ne_0
c  taui=Ti/Te
c  rmin=r/a
c  rmaj=Rmaj/a
c  xnu=nu_ei/(c_s/a)
c  betae=neTe/(B**2/(8pi))  0 is electrostatic
c  shat= dlnr/drho used only for parallel dynamics part
c  alpha local shear parameter or MHD pressure grad (s-alpha diagram)
c  elong= local elongation or kappa
c  xwell amount of magnetic well xwell*min(alpha,1)
c  park=1  (0) is a control parameter to turn on (off) parallel motion
c       0.405 best at zero beta and 2.5x larger at high beta..see iflagin(3)
c  ghat=1  (0) is a control parameter to turn on (off) curvature drift
c  gchat=1 (0) is a control parameter to turn on (off) div EXB motion
c  adamp= radial mode damping exponent  1/4 < adamp < 3/4
c       0.25 from direct fit of simulations varying radial mode damping
c   but 0.75 is better fit to rlti dependence
c  alpha_star O(1-3)  gyyrobohm breaking coef for diamg. rot. shear
c  gamma_star ion diamagnetic rot shear rate in units of c_s/a
c  alpha_e O(1-3)   doppler rot shear coef
c  gamma_e    doppler rot shear rate in units of c_s/a
c  alpha_p 1.5  fit for parallel velocity shear effect at rmaj=3 and q=2
c  gamma_p    parallel velocity shear rate (-d v_phi/ drho) in units of c_s/a
c  kdamp model damping normally 0.
 
c output
 
c  yparam(20) output diagnostics
c kyf  value of ky used
c gamma   leading mode growth rate in c_s/a
c freq    leading mode freq rate in c_s/a
c ph_m    (e phi /T_e)/(rho_s/a)  saturation value
c d_hat    plasma diffusivity for ions
c d_im_hat    plasma diffusivity for impurities
c chii_hat ion energy diffusivity
c chie_hat electron energy diffusivity
c exch_hat anomalous e to i energy exchange
c eta_par_hat parallel component of toroidal momentum diffusivity
c eta_per_hat perpendicular    ""
c eta_phi_hat toroidal momentun diffusivity
 
c internal definitions
c nroot = number of equations,
c   nroot=12 full impurity dynamics
c   nroot=9 exb convective impurity dynamics
c   nroot=8 full pure plasma, nrout=6 (betae=0), nrout=5 (betae=0 and park=0)
c v(i)  12 solution vector
c   v(1)=n_i,v(2)=p_par,v(3)=p_per,v(4)=n_t,v(5)=p_t
c   v(6)=u_par, v(7)=n_u, v(8)=a_par
c   v(9)=n_im, v(10)=p_im_par,v(11)=p_im_per,v(12)=u_im_par
c -i*omega v(i)= sum_j amat(i,j) v(j) where omega=freq+xi*gamma
c quasineitrality is
c  (-idelta+(1-dil)*(1/taui)*(1-g1))*ph=f0*ph=(1-dil)*n_i-n_t-n_u
c  or (1-idelta-reps+(1-dil)*(1/taui)*(1-g1))*ph=f1*ph=(1-dil)*n_i-n_t
c  for betae=0
c k_m  inverse mixing length
 
c numerical controls and flags
c
c
c...Dimensions
c
c neq maximum number of equations
c ieq actual number used
 
c***********************************************************************
c***********************************************************************
c
      use glf
      implicit none
c
c     include 'glf.i'
c
c Glf is common block, which must contain all the _gf inputs and outputs
c
      character*1 jobvr, jobvl
      integer neq, iflagin(30), ilhmax, ilh, ikymaxtot,
     >  lprint, ieq, j1, j2, j, i, jmax,
     >  iar, ifail, jroot(4), itheta, iky, iky0, iroot, iglf
      double precision epsilon
      parameter ( neq = 12, epsilon = 1.D-34 )
c
      double precision pi, xparam(30),yparam(2*nmode),
     >  nroot,ky0,rms_theta,rlti,rlte,rlne,rlni,dil,taui,
     >  rmin,rmaj,q,rlnimp,amassimp,zimp,mimp,
     >  aikymax, aiky, apwt, aiwt,
     >  alpha_mode, gamma_mode, alpha_p, gamma_p,
     >  amassgas, chi_par_1, chi_per_1, x43,
     >  anorm, ave_g, ave_g0,
     >  ave_cos, ave_theta2, ave_kxdky2,
     >  dtheta, theta, phi2, ave_k_par, chk,
     >  alpha_n, yk, byk, fnn, fnp, fnf, fpn, fpp, fpf,
     >  xnu_col, amass_e, chkf, xadiabat,
     >  eta_par_hat, eta_per_hat, anorm_k, del_k,
     >  gamma_k_max, gamma_gross_net,
     >  xnu,betae,shat,alpha,elong,xwell,
     >  park,ghat,gchat,kdamp,
     >  adamp,alpha_star,gamma_star,alpha_e,gamma_e,
     >  kyf,gamma,freq,ph_m,d_hat,d_im_hat,
     >  chii_hat,chie_hat,exch_hat
      double COMPLEX xi, idelta,
     >  v(1:12), amat(1:12,1:12),
     >  n_i,p_par,p_per,n_t,p_t,u_par,n_u,a_par,ph,t_u,n_e,
     >  n_im,p_im_par,p_im_per
c     complex u_im_par
      double precision b0,g0,g1,g2,g3,g12,g23,
     >  b0i,g0i,g1i,g2i,g3i,g12i,g23i
      double COMPLEX f0,f1
      double precision k_par,ky,kx,k_per,k_m,
     >  w_s, w_d, w_d0, w_cd,
     >  reps,xnueff,betae0,k_par0
      double COMPLEX xmu,lamda_d,
     >  xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per
      double precision gam_par,gam_per,x_par,x_per,xt_mhd,yt_mhd,
     >  th,tc,fh,fc,
     >  phi_norm,gamma_r
      double COMPLEX chknu,chknt,chknt2
      double precision phi_renorm,gamma_net
c
c...Declarations for eigenvaluesolver
c
      double precision zgamax
c
c... solver varaibles
c
      parameter ( iar=neq )
c     integer iai, ivr, ivi, intger(neq) ! if NAG solver f02ake used
c     parameter ( iai=neq, ivr=neq, ivi=neq )
      double precision ar(iar,neq), ai(iar,neq), rr(neq), ri(neq)
     &  , vr(iar,neq), vi(iar,neq)
      double precision br(iar,neq), bi(iar,neq), beta_tom(neq), ztemp1
 
      integer matz
      double precision fv1(neq),fv2(neq),fv3(neq)
c
c amat(i,j) = complex matrix A
c zevec(j) = complex eigenvector
c
      integer lwork
      parameter ( lwork=198 )
      double complex mata(iar,neq),cvr(iar,neq),cvl(iar,neq),w(neq)
      double complex work(lwork)
      double precision rwork(2*neq)
      double COMPLEX zevec(neq,neq), zomega(neq)
      double precision gammaroot(4),freqroot(4),phi_normroot(4)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c
c   return zero if no unstable modes
c
      if(ipert_gf.eq.1.and.ngrow_k_gf(0).eq.0)go to 888   
c
c...initialize variables
c
c inputs.........................................................

      do i=1,30
       iflagin(i)=iflagin_gf(i)
       xparam(i)=xparam_gf(i)
      enddo
 
      ilhmax=1
      ikymaxtot=ikymax_gf
c     if (xparam_gf(10).gt.0.) ilhmax=2
c
c If ETG modes included, then double ky spectrum
c Inside ky loop, ilh=1 low k ion modes and ilh=2 high k electron modes
c For ETG mode, use complete te <-> ti and mi <-> me isomorphism
c with adiabatic electron ITG mode then chii_hat will be electron 
c transport in gyrobohm electron units with T_i.
c chie_e_gf converted back to c_s*rho_s**2/a units and added 
c to chie_gf after ky loop
c
      if (xparam_gf(10).gt.0.) then
        ilhmax=2
        ikymaxtot=2*ikymax_gf
      endif
c
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
c
      idelta=0.D0
c     if(ilh.eq.1) idelta=xi*xparam(1)+xparam(2)
 
c.................................................................
c
      if (lprint.gt.0) open(1)
      ieq  = nroot
c
      if (lprint.eq.99) then
      write(1,*) 'ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q: ',
     >    ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q
      write(1,*)'xnu,beta,shat,alpha,elong,xwell: ',
     >    xnu,betae,shat,alpha,elong,xwell
      write(1,*)'park, ghat, gchat: ',
     >    park, ghat, gchat
      write(1,*)'adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp: ',
     >    adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp
 
      endif
      if (lprint.eq.98) then
      write(2,*) 'ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q: ',
     >    ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q
        write(2,*)'xnu,betae,shat,alpha,elong,xwell: ',
     >    xnu,betae,shat,alpha,elong,xwell
        write(2,*)'park, ghat, gchat: ',
     >    park, ghat, gchat
        write(2,*)'adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp: ',
     >    adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp
 
      endif
c
      xi=(0.D0,1.D0)
      pi=atan2 ( 0.0D0, -1.0D0 )
 
c GLF model coefficients
 
      chi_par_1=2.D0*sqrt(2.D0/pi)
      chi_per_1=sqrt(2.D0/pi)
 
      gam_par=3.D0
      gam_per=1.D0
      x_par=2.D0
      x_per=3.D0/2.D0
 
      xmu=(0.80D0+.57D0*xi)
 
      xnu_par_par=(1.D0+xi)
      xnu_par_per=0.D0
      xnu_per_par=0.D0
      xnu_per_per=(1.D0+xi)
 
      lamda_d=(-0.7D0-0.80D0*xi)
      if(xparam(11).ne.0..or.xparam(12).ne.0.)
     >  lamda_d=xparam(11)-xi*xparam(12)
      x43=1.D0
      if (iflagin(4).eq.1) then
       lamda_d=5.D0/3.D0
       x43=4.D0/3.D0
      endif
c
c 3d trial wave function analysis
c
      if(iflagin(2).ge.1) then
       if(rms_theta.eq.0.) rms_theta=pi/3.D0
       if(iflagin(5).eq.1) rms_theta=rms_theta_gf*(2.D0/q_gf)
       if(iflagin(5).eq.2) rms_theta=rms_theta_gf*(2.D0/q_gf)**0.5D0
       if(iflagin(5).eq.3) rms_theta=rms_theta_gf/
     >      (xparam(13)*(q_gf/2.D0-1.D0)+1.D0)
     >  /sqrt(1.D0+xparam(16)*(shat_gf**2-1.D0)+xparam(17)*
     >    (shat_gf-1.D0)**2)
     >  /(1.D0+xparam(18)*sqrt(betae/.006D0))
       if(iflagin(5).eq.4) rms_theta=rms_theta_gf/
     >      (xparam(13)*(q_gf/2.D0-1.D0)+1.D0)
     >  /sqrt(1.D0+xparam(16)*(shat_gf**2-1.0D0)+xparam(17)*
     >    (shat_gf-1.D0)**2+xparam(19)*(alpha-0.5D0)**2.D0)
     >  /(1.D0+xparam(18)*sqrt(betae/.006D0))
       if(iflagin(5).eq.5) rms_theta=rms_theta_gf/
     >      (xparam(13)*((q_gf/2.D0)-1.D0)+1.D0)
     >  /sqrt(1.D0+xparam(16)*((shat_gf-
     >  xparam(19)*alpha)**2-0.5D0)+xparam(17)*
     >  (shat_gf-xparam(19)*alpha-0.5D0)**2)/
     >  taui_gf**0.25D0
c along the filed line physics with wave function
c phi=exp(-theta**2/(4.*rms_theta**2))=W_even
c ave_F= [int(0 to inf) F phi**2 d_theta]/[int(0 to inf)  phi**2 d_theta]
c ave_theta**2=(rms_theta)**2
 
c phi, densities and pressures are even functions of theta
c the odd functions like u_par can be represented by
c W_odd= W*i*theta/rms_theta*W_even
c then for W=-1, the k_par operator = i*k_par=1/(rmaj*q) d/dtheta
c becomes 1/(rmaj*q)/(2.*rms_theta)*park  (park=1) in every equation
c ie ave_k_par=1/(2.*rms_theta)
c park is tuned to best fit.
c
c parallel velocity shear gamma_p breaks parity so wave functions
c become mixed W_even->(1-xi*gamma_p*alpha_n*theta/rms_theta)*W_even
c
c gamma_p*alpha_n mustbe dimensionless and independent of norm a
c hence alpha_n=(rmaj/3)*alpha_p since gamma_p is in units
c of c_s/a  and rmaj=rmajor/a. Since in the slab limit where
c where we can have the parallel shear drive, rmaj enters only with
c parameter rmaj*q, we further assume
c alpha_n=(rmaj/3.)*(q/2)*alpha_p as the appropriate scaling
c q=2 and rmaj=3 are the norm points for alpha_p=1.5
c For the extreme toroidal limit q-> infinity where rmaj and
c q are not product associated, we will lose the instability.
c
c to first order in gamma_p*alpha_n
c this leads to a weighting  factor [gamma_p*alpha_n] in the
c xi*ky*ph1*gamma_p linear drive and in the
c eta_phi_hat=conjg(u_par)*(-xi*ky*ph)/gamma_p toroidal vocosity.
c
c the correct dependence gamma-gamma0 going like gamma_p**2 is found
c but QLT eta_phi_hat goes like gamma_p**2 also
c CAUTION: this is worked out only to first order in gamma_p*alpha_n
c small.  It seems likely that there is a higher order saturation factor
c something like 1/(1+(gamma_p*alpha_n)**2) in  eta_phi_hat
c
c doppler (EXB) rotational shear also breaks parity
c thus there should a term egamma*alpha_n_e added to gamma_p*alpha_n
c but it is unclear how to weight alpha_n_e compared to alpha_n
c
c see R.R.Dominguez and G.M Staebler Phys. Fluids B5 (1993) 3876
c for a discussion of QLT theory of anomalous momentum transport in
c slab geometry
 
c compute weight factors
c fix later so these are computed only once per j grid point
 
c      ave_theta2=rms_theta**2
 
      anorm=0.D0
      ave_g=0.D0
      ave_g0=0.D0
      ave_cos=0.D0
      ave_theta2=0.D0
      ave_kxdky2=0.D0
 
      dtheta=4.D0*rms_theta/100.D0
      theta=0.D0
c
      do itheta=1,100
       theta=theta+dtheta
       phi2=exp(-theta**2/(2.D0*rms_theta**2))
       anorm=anorm+phi2*dtheta
       ave_theta2=ave_theta2+
     >  theta**2*phi2*dtheta
       ave_g=ave_g +
     > (-xwell*min(1.D0,alpha)+cos(theta)+
     >    (shat*theta-alpha*sin(theta))*sin(theta))*phi2*dtheta
       ave_g0=ave_g0 + phi2*dtheta
       ave_kxdky2=ave_kxdky2+
     >  (abs(shat*theta-alpha*sin(theta)))**2*phi2*dtheta
       ave_cos=ave_cos +
     >  cos(theta)*phi2*dtheta
      enddo
c
      ave_theta2=ave_theta2/anorm
      ave_g=ave_g/anorm
      ave_g0=ave_g0/anorm
      ave_kxdky2=ave_kxdky2/anorm
      ave_cos=ave_cos/anorm
c 
      ave_k_par=1/(2.D0*rms_theta)
 
      chk=abs(ave_theta2-rms_theta**2)/rms_theta**2
      if (chk.gt..02) write (6,*) 'chk:', chk
 
      alpha_n=(rmaj/3.D0)*(q/2.D0)*alpha_p
 
      if(lprint.eq.2) then
       write(6,*) 'rms_theta,chk :', rms_theta, chk
       write(6,*) 'ave_theta2,ave_g,ave_k_par,ave_cos:',
     >   ave_theta2,ave_g,ave_k_par,ave_cos
      endif
      endif
      if(iflagin(2).eq.0) then
       shat=1.D0
       ave_theta2=1.D0
       ave_g=1.D0
       ave_g0=1.D0
       ave_kxdky2=1.D0
       ave_k_par=1.D0
       ave_cos=1.D0
      endif
c
c start ky loop 
c first half ITG, second half high-k ETG ... each with ikymax_gf modes
c ilh=1 low k ion modes  ilh=2 high k electron modes

      do iky0=1,ikymaxtot
c
cgms      iky=iky0
      iky = ikymax_gf+1-iky0
      ilh=1
c
c offset iky if in high-k range and set ilh=2
c
      if (iky0.gt.ikymax_gf) then
cgms         iky=iky0-ikymax_gf
         iky = ikymaxtot+1-iky0
         ilh=2
      endif
c
      if (ilh.eq.2) then
       nroot=6
       ieq=nroot
       xnu=0.D0
       betae=1.D-6
       rlte=rlti_gf
       rlti=rlte_gf
       rlne=rlni_gf
       rlni=rlne_gf
       rlnimp=epsilon
       dil=1.D0-1.D0/(1.D0-dil_gf)
       apwt=1.D0
       aiwt=0.D0
       taui=1.D0/taui_gf
       rmin=epsilon
       xparam(7)=0.D0
       xparam(6)=-1.D0
       alpha_star=0.D0
       alpha_e=0.D0
       alpha_p=0.D0
       alpha_n=0.D0
       alpha_mode=0.D0
c check this for current driven mode
      endif
c
      idelta=0.D0
      if (ilh.eq.1) idelta=xi*xparam(1)+xparam(2)
c
c logarithmic ky grid
c
      if(ikymax_gf.gt.1) then
       aikymax=ikymax_gf
       aiky=iky
       yk=aiky/aikymax
 
       byk=log(xkymax_gf/xkymin_gf)/(1.D0-1.D0/aikymax)
       ky=xkymax_gf*exp(byk*(yk-1.D0))
      endif
      if(ikymax_gf.eq.1) then
c     ky=sqrt(2.*taui)/rlti/(rmaj*q)
c     from w_star_ti=v_i_th*k_par
c  possible physics basis of q (ie current) scaling ..to be determined
       if(iflagin(1).eq.0) ky=ky0
       if(iflagin(1).eq.1) ky=ky0*sqrt(taui)*(3.D0/rlti)*
     >    (3.D0*2.D0/rmaj/q)
       if(iflagin(1).eq.2) ky=ky0*(2.D0/q)
       if(ky0.eq.0.) ky=0.3D0
      endif
 
 
      kyf=ky
c
 
      kx=ky*sqrt(ave_kxdky2)
      k_per=sqrt(ky**2+kx**2)
      k_m=sqrt(ky**2+(1.D0-xparam(23))*kx**2) ! inverse mixing length model
 
       do iroot=1,4
        gammaroot(iroot)=0.D0
        freqroot(iroot)=0.D0
        phi_normroot(iroot)=0.D0
       enddo
        d_hat=0.D0
        d_im_hat=0.D0
        chie_hat=0.D0
        chii_hat=0.D0
        exch_hat=0.D0
        eta_par_hat=0.D0
        eta_per_hat=0.D0
        jroot(1)=0
        jroot(2)=0
        jroot(3)=0
c
c skip this k for perturbation if last call was stable
c
      if(ipert_gf.eq.1.and.ngrow_k_gf(iky0).eq.0)go to 777 
c skip this k if the previous k was stable and 4 k's have been done      
      if(iky.lt.ikymax_gf-4.and.ngrow_k_gf(iky0-1).eq.0)go to 777
 
c primary ions
 
      b0=taui*k_per**2
 
c     Pade aproximates...may use gamma functions later
      g0=1.D0
      g1=1.D0/(1+b0)
      g2=1.D0/(1+b0)*g1
      g3=1.D0/(1+b0)*g2
 
      g12=(g1+g2)/2.D0
      g23=(g2+g3)/2.D0
 
c impurity ions
 
      b0i=taui*k_per**2*amassimp/amassgas/zimp**2
 
c     Pade aproximates...may use gamma functions later
      g0i=1.D0
      g1i=1.D0/(1+b0i)
      g2i=1.D0/(1+b0i)*g1i
      g3i=1.D0/(1+b0i)*g2i
 
      g12i=(g1i+g2i)/2.D0
      g23i=(g2i+g3i)/2.D0
 
      mimp=amassimp/amassgas
 
 
      w_s=ky
      w_d=(ghat*2.D0/rmaj)*ky*ave_g
      w_d0=(ghat*2.D0/rmaj)*ky*ave_g0
      w_cd=(gchat*2.D0/rmaj)*ky*ave_g
 
      k_par=park/(rmaj*q)*ave_k_par*sqrt((1.D0+elong**2)/2.D0)
 
c     sqrt((1.+elong**2)/2.) to get higher beta_crit prop to k_par**2
c     roughly same as betae-> betae/((1.+elong**2)/2.)
c     physically like shortening the connection length to good curv.
 
      if (iflagin(3).eq.2) then
       betae=betae_gf/(1.D0+xparam(3))**2/(1.D0+xparam(9))
      endif
 
 
      if (iflagin(3).eq.1) then
       k_par=park/(rmaj*q)*ave_k_par
       betae=betae_gf/(1.D0+xparam(3))**2/(1.D0+xparam(9))
     >      /((1.D0+elong**2)/2.D0)
      endif
 
c     we put the park enhancement directy into betae
c     ie park=.5 best at low beta and 2.5x.5=1.25 at high beta
 
c     option iglagin(3)=1 puts beta_crit elongation enhancement
c     directly into betae
 
c     option iflagin(3)=2 puts beta_crit elongation factor into
c     the connection length
 
c     an extra shape  factor 2 (triangularity) enhancement
c     is optained by (1.+xparam(9))=2.
c     if w_d negative flip sign of dissipative parts
      if(w_d.lt.0.) then
c error 12/21       lamda_d=-conjg(lamda_d)
       lamda_d=conjg(lamda_d)
 
       xmu=-conjg(xmu)
 
       xnu_par_par=-conjg(xnu_par_par)
       xnu_par_per=-conjg(xnu_par_per)
       xnu_per_par=-conjg(xnu_per_par)
       xnu_per_per=-conjg(xnu_par_per)
      endif
 
      reps=(1.D0+xparam(15))*
     >   sqrt((rmin/rmaj)*(1.D0+ave_cos)/(1.D0+(rmin/rmaj)*ave_cos))
 
      if(nroot.le.3) reps=0.D0
 
c fix trapped eletron MHD limit
c 3/4*reps*(1+rlte) + xt_mhd*(1-reps)+yt_mhd*rlte=1+rlte
c solve for xt_mhd and yt_mhd
c 3/4*reps+yt_mhd=1; 3/4*reps+xt_mhd*(1-reps)=1
 
      yt_mhd=(1-x43*(3.D0/4.D0)*reps)
      xt_mhd=(1.D0-x43*(3.D0/4.D0)*reps)/(1.D0-reps)
 
c collision detrapping retrapping model
 
      xnueff=(1.D0+xparam(6))*xnu/(reps**2+1.D-6)
 
c very difficult get xnueff correct hince add enhancement factor
c and fit to Beer or GKS
 
      th=4.08D0
      tc=0.918D0
      fh=0.184D0
      fc=0.816D0
 
      fnn=xnueff*((th/tc**(3.D0/2.D0))-(tc/th**(3.D0/2.D0)))/(th-tc)
      fnp=xnueff*(3.D0/2.D0)*
     >   ((1.D0/th)**(3.D0/2.D0)-(1.D0/tc)**(3.D0/2.D0))/(th-tc)
      fnf=xnueff*((fh/th**(3.D0/2.D0))+(fc/tc**(3.D0/2.D0)))
 
      fpn=xnueff*(2.D0/3.D0)*
     >   ((th/tc**(1.D0/2.D0))-(tc/th**(1.D0/2.D0)))/(th-tc)
      fpp=xnueff*((1.D0/th)**(1.D0/2.D0)-(1.D0/tc)**(1.D0/2.D0))/(th-tc)
      fpf=xnueff*(2.D0/3.D0)*((fh/th**(1.D0/2.D0))+(fc/tc**(1.D0/2.D0)))
 
c  collisional modes added with xnu_col
c  must fix for atomic mass dependence other than deuterium
      xnu_col=xparam(7)
      amass_e=2.7D-4*(2.D0/amassgas)
 
c check adiabatic property that chkf should be 1.0  (finite xnu)
 
      chkf=(fnn*fpp-fpn*fnp)/((fnf*fpp-fpf*fnp)+epsilon)
      if (lprint.eq.2) write(6,*) 'chkf:', chkf
 
      if(neq.le.3) reps=0.D0
 
      f0=-idelta+(1.D0-dil)*apwt*(1/taui)*(1.D0-g1)
     >           +zimp**2*aiwt*(1/taui)*(1.D0-g1i)
      f1=1.D0-reps + f0
 
      xadiabat=0.D0
      if(nroot.le.6) then
        betae=0.D0
        f0=f1
        xadiabat=1.D0
      endif
      if(nroot.le.5) k_par=0.D0
 
      betae0=betae+epsilon
      k_par0=k_par+epsilon
 
      if (lprint.eq.98) then
        write(2,*) 'ky,g1,g2,g3,g12,g23,w_s,w_d,w_cd: ',
     >    ky,g1,g2,g12,g23,w_s,w_d,w_cd
        write(2,*) 'f0,taui,k_par,reps: ',
     >    f0,taui,k_par,k_per,reps
        write(2,*) 'chi_par_1,chi_per_1,gam_par,gam_per:',
     >    chi_par_1,chi_per_1,gam_par,gam_per
        write(2,*) 'x_par,x_per,xmu:',
     >    x_par,x_per,xmu
        write(2,*) 'xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per:',
     >    xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per
        write(2,*) 'lamda_d,betae,xadiabat:',
     >    lamda_d,betae,xadiabat
        write(2,*) 'yt_mhd,xt_mhd:',
     >    yt_mhd,xt_mhd
 
      endif
c 
c matrix in order
c note ph=(n_i-n_t-n_u)/f0 results in (i,1)-(i,4)-(i,7) parts
c
c n_i equ #1
 
      amat(1,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,2)=
     >  +xi*w_d*taui*0.5D0
 
      amat(1,3)=
     >  +xi*w_d*taui*0.5D0
 
      amat(1,4)= -(-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,5)= 0.D0
 
      amat(1,6)=
     >  -xi*k_par
 
      amat(1,7)= -(-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,8)= 0.D0
 
      amat(1,9)= aiwt*zimp*
     >  (-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,10)=0.D0
 
      amat(1,11)=0.D0
 
      amat(1,12)=0.D0
 
c p_par equ #2
 
      amat(2,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
     >  +k_par*chi_par_1
     >  -(xi*w_d*taui*3.D0/2.D0-w_d*taui*xnu_par_par)
     >  -(xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_par_per)
 
      amat(2,2)=
     >  -k_par*chi_par_1
     >  +xi*w_d*taui*x_par +
     >  (xi*w_d*taui*3.D0/2.D0-w_d*taui*xnu_par_par)
 
      amat(2,3)=
     >  (xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_par_per)
 
      amat(2,4)= -(-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,5)=0.D0
 
      amat(2,6)=
     >  -xi*gam_par*k_par
 
      amat(2,7)= -(-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,8)=0.D0
 
      amat(2,9)= aiwt*zimp*
     >  (-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,10)=0.D0
 
      amat(2,11)=0.D0
 
      amat(2,12)=0.D0
 
c p_per equ #3
 
      amat(3,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*((rlni-rlti)*g2+2.D0*rlti*g3)+xi*x_per*w_cd*g23)/f0
     >  +k_par*chi_per_1
     >  -(xi*w_d*taui-w_d*taui*xnu_per_per)
     >  -(xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_per_par)
 
 
      amat(3,2)=
     >  +(xi*w_d*taui*1/2-w_d*taui*xnu_per_par)
 
      amat(3,3)=
     >  -k_par*chi_per_1
     >  +xi*w_d*taui*x_per   +(xi*w_d*taui-w_d*taui*xnu_per_per)
 
      amat(3,4)=
     > -(-xi*w_s*((rlni-rlti)*g2+2.D0*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,5)=0.D0
 
      amat(3,6)=
     >  -xi*gam_per*k_par
 
      amat(3,7)=
     >  -(-xi*w_s*((rlni-rlti)*g2+2.D0*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,8)=0.D0
 
      amat(3,9)= aiwt*zimp*
     >   (-xi*w_s*((rlni-rlti)*g2+2.D0*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,10)=0.D0
 
      amat(3,11)=0.D0
 
      amat(3,12)=0.D0
 
c n_t equ #4
 
      amat(4,1)=(1.D0-dil)*apwt*
     >  (-xi*w_s*rlne*reps*g0+xi*x43*3.D0/4*w_cd*reps*g0)/f0
     >  -(1.D0-dil)*apwt*(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(4,2)=0.D0
 
      amat(4,3)=0.D0
 
      amat(4,4)=-(-xi*w_s*rlne*reps*g0+xi*x43*3.D0/4*w_cd*reps*g0)/f0
     >  -((1.D0-reps)*fnn)
     >  -(-(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat))
 
      amat(4,5)=
     >  -xi*w_d*x43*3.D0/4.D0
     >  -((1.D0-reps)*fnp)
 
      amat(4,6)=0.D0
 
      amat(4,7)=-(-xi*w_s*rlne*reps*g0+xi*x43*3.D0/4*w_cd*reps*g0)/f0
     >   -(-reps*fnf)
 
      amat(4,8)=0.D0
 
      amat(4,9)=aiwt*zimp*
     >   (-xi*w_s*rlne*reps*g0+xi*x43*3.D0/4*w_cd*reps*g0)/f0
     >   -aiwt*zimp*(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(4,10)=0.D0
 
      amat(4,11)=0.D0
 
      amat(4,12)=0.D0
 
c p_t equ #5
 
      amat(5,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.D0/4*w_cd*reps*g0)/f0
     >  -(1.D0-dil)*apwt*(-reps*fpf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(5,2)=0.D0
 
      amat(5,3)=0.D0
 
      amat(5,4)=
     >  -(-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.D0/4*w_cd*reps*g0)/f0
     >            +xi*w_d*lamda_d
     >  -((1.D0-reps)*fpn)
     >  -(-(-reps*fpf*(1.D0-reps)*g0/f0*xadiabat))
 
      amat(5,5)=
     >  -xi*w_d*x43*5.D0/4.D0-xi*w_d*lamda_d
     >  -((1.D0-reps)*fpp)
 
      amat(5,6)=0.D0
 
      amat(5,7)=
     >  -(-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.D0/4*w_cd*reps*g0)/f0
     >  -(-reps*fpf)
 
      amat(5,8)=0.D0
 
      amat(5,9)= aiwt*zimp*
     >  (-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.D0/4*w_cd*reps*g0)/f0
     >  -aiwt*zimp*(-reps*fpf*(1.D0-reps)*g0/f0*xadiabat)

      amat(5,10)=0.D0

      amat(5,11)=0.D0

      amat(5,12)=0.D0

c u_par equ #6
 
      amat(6,1)=(1.D0-dil)*apwt*
     >   (-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0
     >   -(betae/2.D0)*(-xi*k_par*(2.D0/betae0)*g0)/f0)
 
      amat(6,2)=-xi*k_par*taui
 
      amat(6,3)=0.D0
 
      amat(6,4)=
     >  -(-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0)
     >  -(-(betae/2.D0)*(-xi*k_par*(2.D0/betae0)*g0)/f0)
 
      amat(6,5)=0.D0
 
      amat(6,6)=
     > +xi*w_d*(gam_par+gam_per)/2.D0 -w_d*xmu
 
      amat(6,7)=
     >  -(-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0)
     >  -(-(betae/2.D0)*(-xi*k_par*(2.D0/betae0)*g0)/f0)
     >  -(betae/2.D0)*xi*k_par*(2.D0/betae0)/(1.D0-reps)
 
      amat(6,8)=
     >  -(betae/2.D0)*(-xi*w_s*(rlni*g1+rlti*g2))
     >  -(betae/2.D0)*(-xi*w_s*rlne)
     >  +amass_e*xnu*xnu_col*k_per**2
     >  -(betae/2.D0)*(-2.D0/betae0*amass_e*xnu*xnu_col*k_per**2)
c note there is no double counting in last two terms
 
      amat(6,9)=aiwt*zimp*
     >  (-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0
     >  -(betae/2.D0)*(-xi*k_par*(2.D0/betae0)*g0)/f0)
 
      amat(6,10)=0.D0
 
      amat(6,11)=0.D0
 
      amat(6,12)=0.D0
 
c n_u equ #7
 
      amat(7,1)=(1.D0-dil)*apwt*
     >  (-xi*w_s*rlne*(1.D0-reps)*g0+xi*w_cd*
     >  (1.D0-x43*(3.D0/4.D0)*reps)*g0)/f0
     >  +(1.D0-dil)*apwt*(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(7,2)=0.D0
 
      amat(7,3)=0.D0
 
      amat(7,4)=
     >  -(-xi*w_s*rlne*(1.D0-reps)*g0+xi*w_cd*
     >  (1.D0-x43*(3.D0/4.D0)*reps)*g0)/f0
     >  +((1.D0-reps)*fnn)
     >  +(-(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat))
 
      amat(7,5)=0.D0
     >  +((1.D0-reps)*fnp)
 
      amat(7,6)=-xi*k_par
 
      amat(7,7)=
     >  -(-xi*w_s*rlne*(1.D0-reps)*g0+xi*w_cd*
     >  (1.D0-x43*(3.D0/4.D0)*reps)*g0)/f0
     >  -xi*w_d*xt_mhd
     >  +(-reps*fnf)
 
      amat(7,8)=
     >  -xi*k_par*(-k_per**2)
     >  -xi*w_d*yt_mhd*(w_s*(betae/2.D0)/k_par0*rlte)
 
      amat(7,9)=aiwt*zimp*
     >  (-xi*w_s*rlne*(1.D0-reps)*g0+xi*w_cd*
     >  (1.D0-x43*(3.D0/4.D0)*reps)*g0)/f0
     >  +aiwt*zimp*(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(7,10)=0.D0
 
      amat(7,11)=0.D0
 
      amat(7,12)=0.D0
 
c a_par equ #8
 
      amat(8,1)=(1.D0-dil)*apwt*(-xi*k_par*(2.D0/betae0)*g0/f0)
 
      amat(8,2)=0.D0
 
      amat(8,3)=0.D0
 
      amat(8,4)=-(-xi*k_par*(2.D0/betae0)*g0/f0)
 
      amat(8,5)=0.D0
 
      amat(8,6)=0.D0
 
      amat(8,7)=-(-xi*k_par*(2.D0/betae0)*g0/f0)
     >  +xi*k_par*(2.D0/betae0)/(1.D0-reps)
 
      amat(8,8)=-xi*w_s*rlne
     >  -(2.D0/betae0)*amass_e*xnu*xnu_col*(k_per**2)
 
      amat(8,9)=aiwt*zimp*(-xi*k_par*(2.D0/betae0)*g0/f0)
 
      amat(8,10)=0.D0
 
      amat(8,11)=0.D0
 
      amat(8,12)=0.D0
 
c n_im equ #9
 
      amat(9,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,2)=0.D0
 
      amat(9,3)=0.D0
 
      amat(9,4)= -(-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,5)= 0.D0
 
      amat(9,6)= 0.D0
 
      amat(9,7)= -(-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,8)= 0.D0
 
      amat(9,9) = aiwt*zimp*
     >  (-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,10)=
     >  +xi*w_d*taui*0.5D0/zimp
 
      amat(9,11)=
     >  +xi*w_d*taui*0.5D0/zimp
 
      amat(9,12)=
     >  -xi*k_par
 
c pim_par equ #10
 
      amat(10,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
 
      amat(10,2)=0.D0
 
      amat(10,3)=0.D0
 
      amat(10,4)= -(-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
 
      amat(10,5)=0.D0
 
      amat(10,6)=0.D0
 
      amat(10,7)= -(-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
 
      amat(10,8)=0.D0
 
      amat(10,9)= aiwt*zimp*
     >   (-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
     >   +k_par*chi_par_1/sqrt(mimp)
     >   -(xi*w_d*taui*3.D0/2.D0-w_d*taui*xnu_par_par)/zimp
     >   -(xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_par_per)/zimp
 
      amat(10,10)=
     >  -k_par*chi_par_1/sqrt(mimp)
     >  +xi*w_d*taui*x_par/zimp +(xi*w_d*taui*3.D0/2.D0
     >  -w_d*taui*xnu_par_par)/zimp
 
      amat(10,11)=
     >  (xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_par_per)/zimp
 
      amat(10,12)=
     >  -xi*gam_par*k_par
 
c pim_per equ #11
 
      amat(11,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*((rlnimp-rlti)*g2i+2.D0*rlti*g3i)+
     >  xi*x_per*w_cd*g23i)/f0
 
      amat(11,2)= 0.D0
 
      amat(11,3)= 0.D0
 
      amat(11,4)=
     >  -(-xi*w_s*((rlnimp-rlti)*g2i+2.D0*rlti*g3i)+
     >  xi*x_per*w_cd*g23i)/f0
 
      amat(11,5)=0.D0
 
      amat(11,6)=0.D0
 
      amat(11,7)=
     >  -(-xi*w_s*((rlnimp-rlti)*g2i+2.D0*rlti*g3i)+
     >  xi*x_per*w_cd*g23i)/f0
 
      amat(11,8)=0.D0
 
      amat(11,9)= aiwt*zimp*
     >  (-xi*w_s*((rlnimp-rlti)*g2i+2.D0*rlti*g3i)+
     >  xi*x_per*w_cd*g23i)/f0
     >  +k_par*chi_per_1/sqrt(mimp)
     >  -(xi*w_d*taui-w_d*taui*xnu_per_per)/zimp
     >  -(xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_per_par)/zimp
 
      amat(11,10)=
     >  +(xi*w_d*taui*1/2-w_d*taui*xnu_per_par)/zimp
 
      amat(11,11)=
     >  -k_par*chi_per_1/sqrt(mimp)
     >  +xi*w_d*taui*x_per/zimp 
     >  +(xi*w_d*taui-w_d*taui*xnu_per_per)/zimp
 
      amat(11,12)=
     >  -xi*gam_per*k_par
 
c uim_par equ #12
cgms 5/21/99 added gamma_p to amat(12,1),amat(12,4),amat(12,7),amat(12,9)
c    added xnu_col term to amat(12,8)
c    fixed mimp factor in amat(12,12)
 
      amat(12,1)=(1.D0/mimp)*(1.D0-dil)*apwt*
     >  ((-xi*k_par*g1i/f0)*zimp
     >  -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0
     >  -(betae/2.D0)*zimp*(-xi*k_par*(2.D0/betae0)*g0i)/f0)

      amat(12,2)=0.D0
 
      amat(12,3)=0.D0
 
      amat(12,4)=
     > -(1.D0/mimp)*(-xi*k_par*g1i/f0)*zimp
     > -(1.D0/mimp)*(-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0)
     > -(1.D0/mimp)*(-(betae/2.D0)*(-xi*k_par
     > *(2.D0/betae0)*g0i)/f0)*zimp

      amat(12,5)=0.D0
 
      amat(12,6)=0.D0
 
       amat(12,7)=
     > -(1.D0/mimp)*(-xi*k_par*g1i/f0)*zimp
     > -(1.D0/mimp)*(-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0)
     > -(1.D0/mimp)*(-(betae/2.D0)*
     >  (-xi*k_par*(2.D0/betae0)*g0i)/f0)*zimp
     > -(1.D0/mimp)*(betae/2.D0)*xi*
     >  k_par*(2.D0/betae0)/(1.D0-reps)*zimp
 
      amat(12,8)=
     >  -(1.D0/mimp)*(betae/2.D0)*(-xi*w_s*(rlnimp*g1i+rlti*g2i))
     >  -(1.D0/mimp)*(betae/2.D0)*(-xi*w_s*rlne)*zimp
     >  +(1.D0/mimp)*zimp*amass_e*xnu*xnu_col*k_per**2
     >  -(1.D0/mimp)*(betae/2.D0)*
     >   (-2.D0/betae0*amass_e*xnu*xnu_col*k_per**2)*zimp

      amat(12,9)=(1.D0/mimp)*aiwt*zimp*
     >   ((-xi*k_par*g1i/f0)*zimp
     >  -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0
     >  -(betae/2.D0)*zimp*(-xi*k_par*(2.D0/betae0)*g0i)/f0)

      amat(12,10)=-(1.D0/mimp)*xi*k_par*taui
 
      amat(12,11)=0.D0
 
      amat(12,12)= (1.D0/mimp)*(xi*w_d*(gam_par+gam_per)
     >  /2.D0/zimp -w_d*xmu/zimp)
c
c put in rot shear stabilization and possible source of gyrobohm breaking
c and model damping kdamp
c
c***********************************************************************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c solve 12x12 complex
c -xi*omega*v(i)=sum_j amat(i,j)*v(j)  omega=freq+xi*gamma
c upto nroot
c order with max gamma and find eigenvector v(i) with ant fixed norm.
c
c...Fill matricies for eigenvalue equation
c
      do j1=1,neq
        rr(j1) = 0.0D0
        ri(j1) = 0.0D0
        do j2=1,neq
          ar(j1,j2) = REAL(  amat(j1,j2) )
          ai(j1,j2) = aimag( amat(j1,j2) )
c...test tmp
c         ai(j1,j2) = 0.0
c         ar(j1,j2) = 0.0
c         if (j1.eq.j2) ar(j1,j2)=j1
c
          vr(j1,j2) = 0.0D0
          vi(j1,j2) = 0.0D0
        enddo
      enddo
c
c...diagnostic output
c
      if ( lprint .gt. 6 ) then
        write (1,*)
        write (1,*) ' ar(j1,j2)  j2 ->'
        do j1=1,neq
          write (1,192) (ar(j1,j2),j2=1,neq)
        enddo
c
        write (1,*)
        write (1,*) ' ai(j1,j2)  j2->'
        do j1=1,neq
          write (1,192) (ai(j1,j2),j2=1,neq)
        enddo
 192    format (1p8e10.2)
 193    format (1p8e12.4)
      endif
c
c..find the eigenvalues and eigenvectors 
c
c.. eigen_gf = 0 use cgg solver (default)
c..          = 1 use tomsqz solver
c..          = 2 use zgeev solver
c.. not longer used:
c        call f02ake( ar,iar,ai,iai,ieq,rr,ri,vr,ivr,vi,ivi,
c     >               intger, ifail )
c
        ifail = 0
c
        if (eigen_gf .eq. 2 ) then
c
        jobvl = 'N'
        jobvr = 'V'
        do j1=1,neq
         do j2=1,ieq
           mata(j1,j2) = cmplx(ar(j1,j2),ai(j1,j2))
         enddo
        enddo
c
c        call zgeev(jobvl,jobvr,ieq,mata,neq,w,cvl,neq,cvr,
c     &             neq,work,lwork,rwork,ifail)
        do j1=1,neq
         rr(j1) = real(w(j1))
         ri(j1) = aimag(w(j1))
         do j2=1,ieq
           vr(j1,j2) = real(cvr(j1,j2))
           vi(j1,j2) = aimag(cvr(j1,j2))
         enddo
        enddo
c
        elseif (eigen_gf .eq. 1 ) then
c
        do j2=1,neq
           do j1=1,neq
              bi(j1,j2)=0.0D0
              if(j1.eq.j2) then
                 br(j1,j2)=1.0D0
              else
                 br(j1,j2)=0.0D0
              endif
           enddo
        enddo
c
        call r8tomsqz(neq,ieq,ar,ai,br,bi, rr,ri,beta_tom, vr,vi, ifail)
c
        do j1=1,ieq
           ztemp1 = beta_tom(j1)
           if ( abs(beta_tom(j1)) .lt. epsilon ) ztemp1 = epsilon
           rr(j1)=rr(j1) / ztemp1
           ri(j1)=ri(j1) / ztemp1
        enddo
c
        else
c
        matz=1
c       write(*,*) 'neq = ',neq
c       write(*,*) 'ieq = ',ieq
c       write(*,*) 'matz = ',matz
c       write (*,*) ' ar(j1,j2)  j2 ->'
c       do j1=1,neq
c         write (*,193) (ar(j1,j2),j2=1,neq)
c       enddo
c       write (*,*) ' ai(j1,j2)  j2 ->'
c       do j1=1,neq
c         write (*,193) (ai(j1,j2),j2=1,neq)
c       enddo
c
        call cgg_glf(neq,ieq,ar,ai,rr,ri,matz,vr,vi,fv1,fv2,fv3,ifail)
c
c       write (*,*) ' wr(j1) and wi(j1)'
c       do j1=1,neq
c         write (*,193) rr(j1), ri(j1)
c       enddo
c       write (*,*) ' zr(j1,j2)  j2 ->'
c       do j1=1,neq
c         write (*,193) (vr(j1,j2),j2=1,neq)
c       enddo
c       write (*,*) ' zi(j1,j2)  j2 ->'
c       do j1=1,neq
c         write (*,193) (vi(j1,j2),j2=1,neq)
c       enddo

        endif
c
        if ( lprint .gt. 1 ) then
          write (1,*) ifail,' = ifail routine '
        endif
c
c..print eigenvalues
c
        if ( lprint .gt. 6 ) then
          write (1,121)
          do j=1,ieq
            write (1,122) rr(j), ri(j)
          enddo
 121      format (/' Solution of the eigenvalue equations'
     &     /t4,'real   ',t18,'imag   ')
 122      format (1p2e14.5)
        endif
c
c...Store the complex eigenvectors and eigenvalues
c...Note the routines here solve A.v = lambda v
c...but that the equation solved is A.v = -i omega v
c...The i-th column of v is the i-th eigenvector
c
        do j1=1,ieq
          zomega(j1) = xi*(rr(j1)+xi*ri(j1))
          do j2=1,ieq
            zevec(j2,j1) = vr(j2,j1) + xi*vi(j2,j1)
          enddo
        enddo
c
        if ( lprint .gt. 6 ) then
          write (6,123)
          do j=1,ieq
            write (6,122) REAL(zomega(j)), aimag(zomega(j))
          enddo
 123      format (/' Multiplied by i: '
     &     /t4,'zomegar',t18,'zomegai')
        endif
        do iroot=1,4
c
c
c..save growth rates and frequencies in real variables
c
        zgamax = 0.0D0
ctemp
        zgamax = -1.D10
        jmax=0
        gamma=0.D0
        do j=1,ieq
         if(j.ne.jroot(1).and.j.ne.jroot(2).and.j.ne.jroot(3)) then
          if (aimag(zomega(j)).gt. zgamax) then
            zgamax = aimag(zomega(j))
            jmax=j
          endif
         endif
        enddo
c
        freq = REAL( zomega(jmax) )
        gamma = aimag( zomega(jmax) )
c
c skip stable modes
c        if(gamma.lt.0.D0)go to 775

         jroot(iroot)=jmax
 
ctemp        if(zgamax.lt.zepsqrt) gamma=0.
 
        if (jmax.ne.0) then
         gammaroot(iroot)=gamma
         freqroot(iroot)=freq
 
 
        do j=1,12
         v(j)=0.D0
        enddo
        do j=1,ieq
          v(j) = zevec(j,jmax)
        enddo
c
c***********************************************************************
 
      n_i=0.D0
      p_par=0.D0
      p_per=0.D0
      n_t=0.D0
      p_t=0.D0
      u_par=0.D0
      n_u=0.D0
      a_par=0.D0
      n_im=0.D0
      p_im_par=0.D0
      p_im_per=0.D0
c     u_im_par=0.
 
      t_u=0.D0
 
 
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
c     if(ieq.eq.12) u_im_par=v(12)
 
      if (ieq.ge.8) ph=((1.D0-dil)*apwt*n_i+aiwt*zimp*n_im-n_t-n_u)/f0
      if (ieq.lt.8) ph= ((1.D0-dil)*apwt*n_i-n_t)/f0
      if (ieq.le.3) ph= ((1.D0-dil)*apwt*n_i)/f0
      t_u=(betae/2.D0)*w_s/k_par0*rlte*a_par
 
      n_e=(1.D0-dil)*apwt*(n_i-(g0-g1)/taui*ph)
     >     +aiwt*zimp*(n_im-zimp*(g0i-g1i)/taui*ph)
 
c impurity trace convective limit
      if (aiwt.lt.-epsilon) then
       n_im=0.D0
       do j=1,8
        n_im=n_im+amat(9,j)*v(j)/(-xi*freq+gamma)
       enddo
      endif
 
c idelta=xi*yparam(1)+yparam(2)   for trapped electrons
 
      yparam(1)=aimag(-(n_t-reps*ph)/ph)
      yparam(2)=REAL(-(n_t-reps*ph)/ph)
 
 
      chknu=n_u/(1.D0-reps)/ph
      chknt=n_t/reps/ph
      chknt2=n_t*(f0+reps)/(reps*((1.D0-dil)*apwt*n_i+aiwt*zimp*n_im))
      if (lprint.eq.2) write (6,*) 'chknu,chknt,chknt2:',
     >    chknu,chknt,chknt2
 
c non linear saturation rule
 
      gamma_r= 0.2D0*3.D0/2.D0*abs(w_d)*taui   !only scaling important
      if(iglf.eq.1) gamma_r= 0.2D0*3.D0/2.D0*abs(w_d0)*taui
c
      gamma_net=gamma-abs(alpha_star*gamma_star
     >         +alpha_e*gamma_e+alpha_mode*gamma_mode)-kdamp

      gamma_net=max(gamma_net,xparam(8)*gamma)
      if( gamma_net.gt.0.D0)then
       ph_m=gamma_net**(1.D0-adamp)*gamma_r**adamp/(k_m*ky)
c set flag ngrow_k_gf: found at least one unstable mode for this k
       if(ipert_gf.eq.0)ngrow_k_gf(iky0)=1
      endif
 
      if( gamma_net.le.0.) ph_m=0.D0
c
      if(xparam(24).gt.0.D0) then
        if(gamma.gt.0.D0)then
          ph_m=dabs(gamma)**(1.D0-adamp)*gamma_r**adamp/(k_m*ky)/
     >    dsqrt(1.D0+(dabs(alpha_star*gamma_star+
     >    alpha_e*gamma_e+alpha_mode*gamma_mode)/
     >    (dabs(gamma)+.00001D0))**xparam(24))
          if(ipert_gf.eq.0)ngrow_k_gf(iky0)=1
        else
           ph_m=0.D0 
        endif
      endif
 
c 7.17.96
      if(xparam(22).gt.0) then
        if(gamma.gt.0.) ph_m=gamma**(1.D0-adamp-xparam(22))
     >   *gamma_r**adamp/(k_m*ky)
        if(gamma.le.0.) ph_m=0.D0
      endif
 
         phi_norm=0
         phi_normroot(iroot)=0.D0
 
       if( ph_m.gt.0.) then
 
       phi_norm=ph_m*ph_m/ABS((conjg(ph)*ph))
 
       phi_normroot(iroot)=phi_norm
 
c note only real part survives in diffusivities
c    ...units are c_s*rho_s**2/a
c magnetic futter component is too small to worry about
 
      d_hat    = phi_norm*REAL(conjg(n_i)*(-xi*ky*ph))/rlni
     >+d_hat
 
      d_im_hat    = phi_norm*REAL(conjg(n_im)*(-xi*ky*ph))/(rlnimp+
     >epsilon)
     >+d_im_hat
 
      chii_hat = phi_norm*3.D0/2.D0*
     >REAL(conjg((1.D0/3.D0)*p_par+(2.D0/3.D0)*p_per)*(-xi*ky*ph))/rlti
     >+chii_hat
      chii_hat=chii_hat + aiwt/apwt*xparam(21)*phi_norm*3.D0/2.D0*
     >REAL(conjg((1.D0/3.D0)*p_im_par+(2.D0/3.D0)*p_im_per)*
     >   (-xi*ky*ph))/rlti
 
      chie_hat = phi_norm*3.D0/2.D0*
     >REAL(conjg(p_t+n_u+t_u)*(-xi*ky*ph))/rlte
     >+chie_hat
 
c electron to ion energy exchange in units n0*t0*c_s/a*(rho_a/a)**2
c note here we interpret QLT d/dt=-xi*freq dropping gamma part to
c avoid getting nonzero result for n_e->ph adiabatic limit
c ie <(d n_e/dt)*conjg(ph)>_time ave -> 0 adiabatic limit
c note  (-1) means exch_hat is electron to ion rate or
c ion heating rate, ie positive exch_hat cools electrons and heats ions
 
      exch_hat = phi_norm*(-1.D0)*
     >REAL(conjg(-xi*freq*n_e)*ph)
     >+exch_hat
 
      eta_par_hat=(1.D0-xparam(14))*phi_norm*
     >REAL(conjg(u_par)
     >*(-xi*ky*ph))/(gamma_p+epsilon)*(-gamma_p*alpha_n)
     >+xparam(14)*phi_norm*REAL(conjg(
     >-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1*ph/(-xi*freq+gamma))
     >*(-xi*ky*ph))/(gamma_p+epsilon)*(-gamma_p*alpha_n)
     >+eta_par_hat
 
      eta_per_hat = phi_norm*
     >REAL(conjg(-ky*(ky*shat*rms_theta)*ph)*
     >   (ph+taui*((1.D0/3.D0)*p_par+(2.D0/3.D0)*p_per)))*
     >   (-gamma_p*alpha_n)
     >/(gamma_p+epsilon)
     >+eta_per_hat
 
       endif
       endif
      enddo
 777  continue
c
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
 
 
ctest       exch_gf=-freq_gf(1)/xkyf_gf*diff_gf*rlni
 
       eta_par_k_gf(iky)=eta_par_hat
       eta_per_k_gf(iky)=eta_per_hat
 
c b_pol/b_phi=rmin/(rmaj*q)
       eta_phi_k_gf(iky)=eta_par_k_gf(iky)+
     >  rmin/(rmaj*q)*eta_per_k_gf(iky)
 
      endif
c
c computed high-k ETG electron transport at each ky
c Note: not added to ITG chi-e here ... done after ky loop
      chie_e_k_gf(iky)=0.D0
      if(ilh.eq.2) then
        chie_e_k_gf(iky)=xparam(10)*chii_hat*
     >                   taui_gf**(3.D0/2.D0)/
     >                   (1836.D0*amassgas_gf)**.5D0
      endif
c end ky loop
      enddo
 
c end big loop on ilh ... no longer used
c     enddo
c
c check to see if any unstable modes were found
c
      if(ipert_gf.eq.0)then
        do j=1,nmode
         if(ngrow_k_gf(j).ne.0)ngrow_k_gf(0)=1
        enddo
        if(ngrow_k_gf(0).eq.0)go to 888
      endif
c
c
c***********************************************************************
c initializations for summations over ky
c
      anorm_k=0.D0
      diff_gf=0.D0
      diff_im_gf=0.D0
      chii_gf=0.D0
      chie_gf=0.D0
      exch_gf=0.D0
      eta_par_gf=0.D0
      eta_per_gf=0.D0
      eta_phi_gf=0.D0
      chie_e_gf=0.D0
c
c Sum ITG and ETG transport
c over logarithmic ky grid (d ky=ky*d yk)
c
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
c
c Add ITG and ETG electron transport
c
      chie_gf=chie_gf + 1.D0*chie_e_gf
c
      diff_gf=diff_gf/anorm_k
      diff_im_gf=diff_im_gf/anorm_k
      chii_gf=chii_gf/anorm_k
      chie_gf=chie_gf/anorm_k
      exch_gf=exch_gf/anorm_k
      eta_par_gf=eta_par_gf/anorm_k
      eta_per_gf=eta_per_gf/anorm_k
      eta_phi_gf=eta_phi_gf/anorm_k
      chie_e_gf=chie_e_gf/anorm_k
c
c
c pick off maximum gamma
c 
      do iroot=1,4
       gamma_k_max=-1.D6
       do iky=1,ikymax_gf
        if(gamma_k_gf(iroot,iky).gt.gamma_k_max) then
         gamma_k_max=gamma_k_gf(iroot,iky)
         gamma_gf(iroot)=gamma_k_gf(iroot,iky)
         freq_gf(iroot)=freq_k_gf(iroot,iky)
         xky_gf(iroot)=xkyf_k_gf(iky)
        endif
       enddo
      enddo
c
c pick off 2nd maximum gamma
c
       gamma_k_max=-1.D6
       do iky=1,ikymax_gf
        if( (gamma_k_gf(1,iky).gt.gamma_k_max) .and.
     >      (gamma_k_gf(1,iky).lt.gamma_gf(1)) ) then
         gamma_k_max=gamma_k_gf(1,iky)
         gamma_gf(2)=gamma_k_gf(1,iky)
         freq_gf(2)=freq_k_gf(1,iky)
         xky_gf(2)=xkyf_k_gf(iky)
        endif
       enddo
c
c       write(6,*) gamma_gf(1), gamma_gf(2), xky_gf(1), xky_gf(2)
c
c print to file log
c      write(*,66)chii_gf,(gamma_gf(j),j=1,4)
 66    format(f14.9,4f14.9)
 67    format(2i2,f14.9)

 
      if(xparam(22).gt.0.) then
       phi_renorm=1.D0
       gamma_gross_net=gamma_gf(1)-abs(alpha_star*gamma_star
     >         +alpha_e*gamma_e+alpha_mode*gamma_mode)-kdamp
       if(gamma_gross_net.gt.0.)
     >  phi_renorm=gamma_gross_net**(xparam(22)*2.D0)
       if(gamma_gross_net.le.0.) phi_renorm=0.D0
 
      diff_gf=diff_gf*phi_renorm
      diff_im_gf=diff_im_gf*phi_renorm
      chii_gf=chii_gf*phi_renorm
      chie_gf=chie_gf*phi_renorm
      exch_gf=exch_gf*phi_renorm
      eta_par_gf=eta_par_gf*phi_renorm
      eta_per_gf=eta_per_gf*phi_renorm
      eta_phi_gf=eta_phi_gf*phi_renorm
      chie_e_gf=chie_e_gf*1.0D0
 
 
      endif
 
c put in cnorm_gf 12/22/96
 
 
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
        write(6,*) 'rlti,rlte,rlne,rlni,rlnimp: ',
     >    rlti,rlte,rlne,rlni,rlnimp
       write(6,*) 'chii,chie,diff,diff_im: ',
     >    chii_gf,chie_gf,diff_gf,diff_im_gf
       write(6,*) 'gamma_gf,freq_gf,ph_m: ',
     >    gamma_gf,freq_gf,ph_m
       write(6,*) 'jmax: ',
     >    jmax
        write(2,*) 'rlti,rlte,rlne,rlni,rlnimp: ',
     >    rlti,rlte,rlne,rlni,rlnimp
       write(2,*) 'chii,chie,diff,diff_im: ',
     >    chii_gf,chie_gf,diff_gf,diff_im_gf
       write(2,*) 'gamma_gf,freq_gf,ph_m: ',
     >    gamma_gf,freq_gf,ph_m
       write(2,*) 'jmax: ',
     >    jmax
 
        write (2,*) ' ar(j1,j2)  j2 ->'
        do j1=1,neq
          write (2,*) (ar(j1,j2),j2=1,neq)
        enddo
c
        write (2,*)
        write (2,*) ' ai(j1,j2)  j2->'
        do j1=1,neq
          write (2,*) (ai(j1,j2),j2=1,neq)
        enddo
c
        write (2,*) ' vr(j1,j2)  j2 ->'
        do j1=1,neq
          write (2,*) (vr(j1,j2),j2=1,neq)
        enddo
c
        write (2,*)
        write (2,*) ' vi(j1,j2)  j2->'
        do j1=1,neq
          write (2,*) (vi(j1,j2),j2=1,neq)
        enddo
 
      endif
 
 999  continue
      if (lprint.gt.0) close(1)

      return
c
c return for case with  no unstable modes
c
 888  continue  
      diff_gf=0.D0
      diff_im_gf=0.D0
      chii_gf=0.D0
      chie_gf=0.D0
      exch_gf=0.D0
      eta_par_gf=0.D0
      eta_per_gf=0.D0
      eta_phi_gf=0.D0
      chie_e_gf=0.D0
      do j1=1,4
        gamma_gf(j1)=0.D0
        freq_gf(j1)=0.D0
        xky_gf(j1)=0.D0
      enddo
      return

      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccgg
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c 
      subroutine cgg_glf(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)

      integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fv3(n)

c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.

c     on input

c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.

c        n  is the order of the matrix  a=(ar,ai).

c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.

c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.

c     on output

c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.

c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.

c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.

c        fv1, fv2, and  fv3  are temporary storage arrays.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50

   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)

      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s

c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by  cbal.

c        scale contains information determining the permutations
c          and scaling factors used by  cbal.

c        m is the number of eigenvectors to be back transformed.

c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.

c     on output

c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120

      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.000/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue

  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140

         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue

  140 continue

  200 return
      end
      subroutine cbal(nm,n,ar,ai,low,igh,scale)

      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv

c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.

c     on output

c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.

c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.

c        scale contains information determining the
c           permutations and scaling factors used.

c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.

c     note that 1 is returned for igh if igh is zero formally.

c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)

c     arithmetic is real throughout.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      radix = 16.000

      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50

      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue

      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue

   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj

         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.000 .or. ai(j,i) .ne. 0.000) go to 120
  110    continue

         m = l
         iexc = 1
         go to 20
  120 continue

      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1

  140 do 170 j = k, l

         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.000 .or. ai(i,j) .ne. 0.000) go to 170
  150    continue

         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.000
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.

      do 270 i = k, l
         c = 0.000
         r = 0.000

         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(ar(j,i)) + abs(ai(j,i))
            r = r + abs(ar(i,j)) + abs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.000 .or. r .eq. 0.000) go to 270
         g = r / radix
         f = 1.000
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.000 / f
         scale(i) = scale(i) * f
         noconv = .true.

         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue

         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue

  270 continue

      if (noconv) go to 190

  280 low = k
      igh = l
      return
      end
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci

c     complex division, (cr,ci) = (ar,ai)/(br,bi)

      double precision s,ars,ais,brs,bis
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)

      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag,dlapy3gf

c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.

c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.

c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.

c     on output

c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.

c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.

c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.

c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  sqrt(a*a + b*b) .

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1

      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.000) go to 170
         norm = dlapy3gf(hr(i,i-1),hi(i,i-1))
crew inserted norm+1.d-100
         yr = hr(i,i-1) / (norm+1.d-100)
         yi = hi(i,i-1) / (norm+1.d-100)
         hr(i,i-1) = norm
         hi(i,i-1) = 0.000

         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue

         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue

  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue

      en = igh
      tr = 0.000
      ti = 0.000
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.000 .and. xi .eq. 0.000) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.000
      yi = (hi(enm1,enm1) - si) / 2.000
      call csroot(yr**2-yi**2+xr,2.000*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.000) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.000

  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue

      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1

      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.000
         norm = dlapy3gf(dlapy3gf(hr(i-1,i-1),hi(i-1,i-1)),sr)
crew inserted norm+1.d-100
         xr = hr(i-1,i-1) / (norm+1.d-100)
         wr(i-1) = xr
         xi = hi(i-1,i-1) / (norm+1.d-100)
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.000
         hi(i,i-1) = sr / (norm+1.d-100)

         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue

  500 continue

      si = hi(en,en)
      if (si .eq. 0.000) go to 540
      norm = dlapy3gf(hr(en,en),si)
crew inserted norm+1.d-100
      sr = hr(en,en) / (norm+1.d-100)
      si = si / (norm+1.d-100)
      hr(en,en) = norm
      hi(en,en) = 0.000
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)

         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.000
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue

  600 continue

      if (si .eq. 0.000) go to 240

      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue

      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)

      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag, dlapy3gf

c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.

c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.

c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.000 for these elements.

c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.

c     on output

c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.

c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.

c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.

c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.

c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  sqrt(a*a + b*b) .

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated october 1989.

c     ------------------------------------------------------------------

      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n

         do 100 i = 1, n
            zr(i,j) = 0.000
            zi(i,j) = 0.000
  100    continue
         zr(j,j) = 1.000
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.000 .and. orti(i) .eq. 0.000) go to 140
         if (hr(i,i-1) .eq. 0.000 .and. hi(i,i-1) .eq. 0.000) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1

         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue

         do 130 j = i, igh
            sr = 0.000
            si = 0.000
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
crew inserted norm+1.d-100
            sr = sr / (norm+1.d-100)
            si = si / (norm+1.d-100)

            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue

  130    continue

  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1

      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.000) go to 170
         norm = dlapy3gf(hr(i,i-1),hi(i,i-1))
crew     inserted norm+1.d-100
         yr = hr(i,i-1) / (norm+1.d-100)
         yi = hi(i,i-1) / (norm+1.d-100)
         hr(i,i-1) = norm
         hi(i,i-1) = 0.000

         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue

         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue

         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue

  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue

      en = igh
      tr = 0.000
      ti = 0.000
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.000 .and. xi .eq. 0.000) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.000
      yi = (hi(enm1,enm1) - si) / 2.000
      call csroot(yr**2-yi**2+xr,2.000*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.000) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.000

  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue

      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1

      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.000
         norm = dlapy3gf(dlapy3gf(hr(i-1,i-1),hi(i-1,i-1)),sr)
crew inserted norm+1.d-100
         xr = hr(i-1,i-1) / (norm+1.d-100)
         wr(i-1) = xr
         xi = hi(i-1,i-1) / (norm+1.d-100)
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.000
         hi(i,i-1) = sr / (norm+1.d-100)

         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue

  500 continue

      si = hi(en,en)
      if (si .eq. 0.000) go to 540
      norm = dlapy3gf(hr(en,en),si)
crew inserted norm+1.d-100
      sr = hr(en,en) / (norm+1.d-100)
      si = si / (norm+1.d-100)
      hr(en,en) = norm
      hi(en,en) = 0.000
      if (en .eq. n) go to 540
      ip1 = en + 1

      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)

         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.000
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue

         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue

  600 continue

      if (si .eq. 0.000) go to 240

      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue

      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue

      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.000

      do 720 i = 1, n

         do 720 j = i, n
            tr = abs(hr(i,j)) + abs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue

      if (n .eq. 1 .or. norm .eq. 0.000) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.000
         hi(en,en) = 0.000
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.000
            zzi = 0.000
            ip1 = i + 1

            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue

            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.000 .or. yi .ne. 0.000) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = abs(hr(i,en)) + abs(hi(i,en))
            if (tr .eq. 0.000) go to 780
            tst1 = tr
            tst2 = tst1 + 1.000/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue

  780    continue

  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840

         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue

  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)

         do 880 i = low, igh
            zzr = 0.000
            zzi = 0.000

            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue

            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue

      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)

      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale,pythag,dlapy3gf

c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.

c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.

c     on output

c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.

c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.

c     calls pythag for  sqrt(a*a + b*b) .

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200

      do 180 m = kp1, la
         h = 0.000
         ortr(m) = 0.000
         orti(m) = 0.000
         scale = 0.000
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))

         if (scale .eq. 0.000) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue

         g = sqrt(h)
         f = dlapy3gf(ortr(m),orti(m))
         if (f .eq. 0.000) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.000 + g) * ortr(m)
         orti(m) = (1.000 + g) * orti(m)
         go to 105

  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.000
            fi = 0.000
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue

            fr = fr / h
            fi = fi / h

            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue

  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.000
            fi = 0.000
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue

            fr = fr / h
            fi = fi / h

            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue

  160    continue

         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue

  200 return
      end
      subroutine csroot(xr,xi,yr,yi)
      double precision xr,xi,yr,yi

c     (yr,yi) = complex sqrt(xr,xi)
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)

      double precision s,tr,ti,pythag,dlapy3gf
      tr = xr
      ti = xi
      s = sqrt(0.5d0*(dlapy3gf(tr,ti) + abs(tr)))
      if (tr .ge. 0.000) yr = s
      if (ti .lt. 0.000) s = -s
      if (tr .le. 0.000) yi = s
      if (tr .lt. 0.000) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.000) yi = 0.5d0*(ti/yr)
      return
      end
      double precision function pythag(a,b)
      double precision a,b

c     finds sqrt(a**2+b**2) without overflow or destructive underflow

      double precision p,r,s,t,u
crew changed dmax1 to max
      p = max(abs(a),abs(b))
      if (p .eq. 0.000) go to 20
crew changed dmin1 to min
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.000 + r
c        write(*,*) 't = ',t
         if (abs(t-4.000) .lt. 1.e-5) go to 20
         s = r / t
         u = 1.000 + 2.000*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c
      DOUBLE PRECISION FUNCTION DLAPY3GF( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
*     ..
*
*  Purpose
*  =======
*
*  DLAPY3GF returns sqrt(x**2+y**2+z**2), taking care not to cause
*  unnecessary overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*  Z       (input) DOUBLE PRECISION
*          X, Y and Z specify the values x, y and z.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      Z = 0
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
         DLAPY3GF = ZERO
      ELSE
         DLAPY3GF = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY3GF
c
      end
c@callglf2d.f
c 12-feb-03 Kinsey, v1.61 retuned GLF23 model
c 11-apr-01 Kinsey, updated for v1.50
c 29-aug-00 Kinsey, added ave_ve for 3pt smoothing of ve (0 for none)
c 05-nov-99 Kinsey, precall routine for glf2d.f
c added i_dengrad switch for dilution
************************************************************************
       subroutine callglf2d(
     >                 !INPUTS
     > leigen,         ! eigenvalue solver
     >                 ! 0 for cgg (default), 1 for tomsqz, 2 for zgeev
     > nroot,          ! no. roots,8 for default, 12 for impurity dynamics
     > iglf,           ! 0 for original GLF23, 1 for retuned version
     > jshoot,         ! jshoot=0 time-dep code;jshoot=1 shooting code
     > jmm,            ! grid number;jmm=0 does full grid jm=1 to jmaxm-1
     > jmaxm,          ! profile grids 0 to jmaxm
     > itport_pt,      ! 1:5 transport flags
     > irotstab,       ! 0 to use egamma_exp; 1 use egamma_m
     > te_m,           ! 0:jmaxm te Kev           itport_pt(2)=1 transport
     > ti_m,           ! 0:jmaxm ti Kev           itport_pt(3)=1 transport
     > ne_m,           ! 0:jmaxm ne 10**19 1/m**3
     > ni_m,           ! 0:jmaxm ni 10**19 1/m**3 itport_pt(1)=1 transport
     > ns_m,           ! 0:jmaxm ns 10**19 1/m**3 
     > i_grad,         ! default 0, for D-V method use i_grad=1 to input gradients
     > idengrad,       ! default 2, for simple dilution
     > zpte_in,        ! externally provided log gradient te w.r.t rho (i_grad=1)
     > zpti_in,        ! externally provided log gradient ti w.r.t rho
     > zpne_in,        ! externally provided log gradient ne w.r.t rho
     > zpni_in,        ! externally provided log gradient ni w.r.t rho
     > angrotp_exp,    ! 0:jmaxm exp plasma toroidal angular velocity 1/sec
     >                 ! if itport_pt(4)=0 itport_pt(5)=0
     > egamma_exp,     ! 0:jmaxm exp exb shear rate in units of csda_exp
     >                 ! if itport_pt(4)=-1 itport_pt(5)=0
     > gamma_p_exp,    ! 0:jmaxm exp par. vel. shear rate in units of csda_exp
     >                 ! if itport_pt(4)=-1 itport_pt(5)=0
     > vphi_m,         ! 0:jmaxm toroidal velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=0 otherwise output
     > vpar_m,         ! 0:jmaxm parallel velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=1 otherwise output
     > vper_m,         ! 0:jmaxm perp. velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=1 otherwise output
     > zeff_exp,       ! 0:jmaxm ne in 10**19 1/m**3
     > bt_exp,         ! vaccuum axis toroidal field in tesla
     > bt_flag,        ! switch for effective toroidal field use in rhosda
     > rho,            ! 0:jmaxm 0 < rho < 1 normalized toroidal flux (rho=rho/rho(a))
     > arho_exp,       ! rho(a), toroidal flux at last closed flux surface (LCFS)
     >                 !   toroidal flux= B0*rho_phys**2/2 (m)
     >                 !   B0=bt_exp, arho_exp=rho_phys_LCFS
     > gradrho_exp,    ! 0:jmaxm dimensionless <|grad rho_phys |**2>
     > gradrhosq_exp,  ! 0:jmaxm dimensionless <|grad rho_phys |>
     >                 !NOTE:can set arho_exp=1.,if gradrho_exp=<|grad rho |>
     >                 !                 and gradrhosq_exp = <|grad rho |**2>
     > rmin_exp,       ! 0:jmaxm minor radius in meters
     > rmaj_exp,       ! 0:jmaxm major radius in meters
     > rmajor_exp,     ! axis major radius
     > zimp_exp,       ! effective Z of impurity
     > amassimp_exp,   ! effective A of impurity
     > q_exp,          ! 0:jmaxm safety factor
     > shat_exp,       ! 0:jmaxm magnetic shear, d (ln q_exp)/ d (ln rho)
     > alpha_exp,      ! 0:jmaxm MHD alpha from experiment
     > elong_exp,      ! 0:jmaxm elongation
     > amassgas_exp,   !  atomic number working hydrogen gas
     > alpha_e,        ! 1 full (0 no) no ExB shear stab
     > x_alpha,        ! 1 full (0 no) alpha stabilization  with alpha_exp
     >                 !-1 full (0 no) self consistent alpha_m stab.
     > i_delay,        !i_delay time delay for ExB shear should be non-zero only
     >                 ! once per step and is less than or equal 10
     >                 !OUTPUTS
     > diffnem,        ! ion plasma diffusivity in m**2/sec
     > chietem,        ! electron ENERGY diffuivity in m**2/sec
     > chiitim,        ! ion      ENERGY diffuivity in m**2/sec
     > etaphim,        ! toroidal velocity diffusivity in m**2/sec
     > etaparm,        ! parallel velocity diffusivity in m**2/sec
     > etaperm,        ! perpendicular velocity diffusivity in m**2/sec
     > exchm,          ! turbulent electron to ion ENERGY exchange in MW/m**3
     >                 ! 0:jmaxm values
     >  diff_m,
     >  chie_m,
     >  chii_m,
     >  etaphi_m,
     >  etapar_m,
     >  etaper_m,
     >  exch_m,
     >
     >  egamma_m,      !0:jmaxm exb shear rate in units of local csda_m
     >  egamma_d,      !0:jmaxm exb shear rate delayed by i_delay steps
     >  gamma_p_m,     !0:jmaxm par. vel. shear rate in units of local csda_m
     >  anrate_m,      !0:jmaxm leading mode rate in unints of local csda_m
     >  anrate2_m,     !0:jmaxm 2nd mode rate in units of local csda_m
     >  anfreq_m,      !0:jmaxm leading mode frequency
     >  anfreq2_m      !0:jmaxm 2nd mode frequency
     > )
c************************************************************************
 
c see glf2d documentation for use of diffusivities in transport equations
c ITER definitions of diffusivity used.
c must add neoclassical diffusion
 
c.......begin common block ....................
c      only common used for communication with glf2d....designation by xxxxx_gf
 
      use glf
      implicit none

c     include 'glf.i'
 
c.......end common block....................
 
c.......begin dimensions.................
 
      double precision epsilon, zeps, zpi
      parameter ( epsilon = 1.D-34, zeps = 1.D-6 )
 
c external arrays
 
      integer jmaxm, jshoot, jmm, i_grad, idengrad, itport_pt(1:5),
     >  i_delay, j, jin, jout, jm, irotstab, iglf,
     >  jpt, jptf, jptt, jj, ii, ik, bt_flag, leigen, nroot
      double precision alpha_e, x_alpha, xalpha,
     >  diffnem, chietem, chiitim, etaphim,
     >  etaparm, etaperm, exchm,
     >  rmajor_exp, zimp_exp, amassimp_exp, 
     >  bt_exp, arho_exp, amassgas_exp,
     >  cbetae, cxnu, relx, cmodel, drho, zeff_e,
     >  zpmte, zpmti, zpmne, zpmni, vstar_sign,
     >  egeo_local, pgeo_local, rdrho_local, rdrho_local_p1, fc,
     >  akappa1, alpha_neo, alpha_neo_hold,
     >  zpmne_q, zpmni_q, zpmnimp, gfac
      double precision te_m(0:jmaxm),ti_m(0:jmaxm),
     >  ne_m(0:jmaxm),ni_m(0:jmaxm), ns_m(0:jmaxm),
     >  vphi_m(0:jmaxm),angrotp_exp(0:jmaxm),
     >  egamma_exp(0:jmaxm),gamma_p_exp(0:jmaxm),
     >  vpar_m(0:jmaxm),vper_m(0:jmaxm),
     >  rho(0:jmaxm),rmin_exp(0:jmaxm),rmaj_exp(0:jmaxm),
     >  gradrho_exp(0:jmaxm),gradrhosq_exp(0:jmaxm),
     >  zeff_exp(0:jmaxm),q_exp(0:jmaxm),shat_exp(0:jmaxm),
     >  bteff_exp(0:jmaxm)
      double precision alpha_exp(0:jmaxm),elong_exp(0:jmaxm),
     >  diff_m(0:jmaxm),chie_m(0:jmaxm),chii_m(0:jmaxm),
     >  etaphi_m(0:jmaxm),
     >  etapar_m(0:jmaxm),etaper_m(0:jmaxm), exch_m(0:jmaxm),
     >  egamma_m(0:jmaxm),egamma_d(0:jmaxm,10),gamma_p_m(0:jmaxm),
     >  anrate_m(0:jmaxm), anrate2_m(0:jmaxm),
     >  anfreq_m(0:jmaxm), anfreq2_m(0:jmaxm)
 
c internal arrays (which can be converted to externals)
 
      double precision zpte_in, zpti_in, zpne_in, zpni_in,
     >  zpte_m(0:jmaxm),zpti_m(0:jmaxm),
     >  zpne_m(0:jmaxm),zpni_m(0:jmaxm),
     >  drhodr(0:jmaxm),drhodrrrho(0:jmaxm),geofac(0:jmaxm),
     >  rhosda_m(0:jmaxm),csda_m(0:jmaxm),cgyrobohm_m(0:jmaxm),
     >  betae_m(0:jmaxm),xnu_m(0:jmaxm),
     >  alpha_m(0:jmaxm),vstarp_m(0:jmaxm)
 
c working arrays and variables
 
      double precision ve(0:jmaxm),vpar(0:jmaxm)
c     real*8 vmode(0:jmaxm)
c     real*8 kevdsecpmw
 
c diagnostic arrays (if used)
 
c     real*8 vstar_m(0:jmaxm),vexb_m(0:jmaxm)
c     real*8 vmode_m(0:jmaxm)
c     real*8 gamma_mode_m(0:jmaxm)
c     real*8 gamma_k_j(20,0:jmaxm),freq_k_j(20,0:jmaxm)
c     real*8 chie_k_j(20,0:jmaxm),chii_k_j(20,0:jmaxm)
c     real*8 vnewstare_m(0:jmaxm),vnewstari_m(0:jmaxm)
c     real*8 ky_j(0:jmaxm)
c     real*8 gamma_j(0:jmaxm,1:4),freq_j(0:jmaxm,1:4)
c     real*8 phi_norm_j(0:jmaxm,1:4)
c     real*8 dnrate_m(0:jmaxm), dtnrate_m(0:jmaxm)
c     real*8 dnfreq_m(0:jmaxm)
 
c some internals

      double precision tem,tim,nem,nim,nsm,zeffm, aiwt_jp1, 
     >       xnimp_jp1, xnimp, vnewk3x

c..Controls for diagnostic printout implemented by Bateman, 15 March 1009

c  initial = 0 for the first step and then is incremented upward
c  iprint  > 0 for diagnostic printout

      INTEGER initial, iprint

      data initial /0/
      data iprint /0/
 
c.......end   dimensions.................
 
c    jm is local grid 0 < jin_m < j < jout_m < jmaxm
c    jm must be greater than 0 and less than jmaxm
c
c
c...constants
      zpi = atan2 ( 0.0D0, -1.0D0 )
c
cdmc      write(6,7701) zpi
cdmc 7701 format(' zpi = ',1pe17.10)
cdmc      write(6,7702) jshoot,jmm,jmaxm,(itport_pt(j),j=1,5)
cdmc 7702 format(/' jshoot,jmm,jmaxm = ',3(1x,i6)/' itport_pt = ',5(1x,i3))
cdmc      call echo('te_m',te_m,jmaxm)
cdmc      call echo('ti_m',ti_m,jmaxm)
cdmc      call echo('ne_m',ne_m,jmaxm)
cdmc      call echo('ni_m',ni_m,jmaxm)
c
cdmc      write(6,7703) i_grad,zpte_in,zpti_in,zpne_in,zpni_in
cdmc 7703 format(/' igz: ',i1,1x,4(1pe17.10))
c
cdmc      call echo('angrotp_exp',angrotp_exp,jmaxm)
cdmc      call echo('egamma_exp',egamma_exp,jmaxm)
cdmc      call echo('gamma_p_exp',gamma_p_exp,jmaxm)
cdmc      call echo('vphi_m',vphi_m,jmaxm)
cdmc      call echo('vpar_m',vpar_m,jmaxm)
cdmc      call echo('vper_m',vper_m,jmaxm)
cdmc      call echo('zeff_exp',zeff_exp,jmaxm)
c
cdmc      write(6,7705) bt_exp
cdmc 7705 format(/' bt_exp = ',1pe17.10)
c
cdmc      call echo('rho',rho,jmaxm)
c
cdmc      write(6,7706) arho_exp
cdmc 7706 format(/' arho_exp = ',1pe17.10)
c
cdmc      call echo('gradrho_exp',gradrho_exp,jmaxm)
cdmc      call echo('gradrhosq_exp',gradrhosq_exp,jmaxm)
cdmc      call echo('rmin_exp',rmin_exp,jmaxm)
cdmc      call echo('rmaj_exp',rmaj_exp,jmaxm)
c
cdmc      write(6,7707) rmajor_exp
cdmc 7707 format(/' rmajor_exp = ',1pe17.10)
c
cdmc      call echo('q_exp',q_exp,jmaxm)
cdmc      call echo('shat_exp',shat_exp,jmaxm)
cdmc      call echo('alpha_exp',alpha_exp,jmaxm)
cdmc      call echo('elong_exp',elong_exp,jmaxm)
c
cdmc      write(6,7708) amassgas_exp
cdmc 7708 format(/' atomic no.:  ',1pe17.10)
c
cdmc      write(6,7709) alpha_e,x_alpha
cdmc 7709 format(/' alpha_e, x_alpha = ',2(1x,1pe17.10))
c
cdmc      write(6,7710) i_delay
cdmc 7710 format(/' i_delay = ',i5)
c
c...initialize variables
c
      do j=0,jmaxm
        zpte_m(j) = 0.D0
        zpti_m(j) = 0.D0
        zpne_m(j) = 0.D0
        zpni_m(j) = 0.D0
 
        betae_m(j)     = 0.D0
        xnu_m(j)       = 0.D0
        cgyrobohm_m(j) = 0.D0
        rhosda_m(j)    = 0.D0
        csda_m(j)      = 0.D0
 
        geofac(j)      = 0.D0
        drhodr(j)      = 0.D0
        drhodrrrho(j)  = 0.D0
 
        gamma_p_m(j)   = 0.D0
        egamma_m(j)    = 0.D0
        ve(j)          = 0.D0
        vper_m(j)      = 0.D0
        vpar(j)        = 0.D0
c        vphi_m(j)      = 0.D0
        vstarp_m(j)    = 0.D0
        alpha_m(j)     = 0.D0
 
        anrate_m(j)    = 0.D0
        anrate2_m(j)   = 0.D0
        anfreq_m(j)    = 0.D0
        anfreq2_m(j)   = 0.D0
 
        exch_m(j)      = 0.D0
        diff_m(j)      = 0.D0
        chie_m(j)      = 0.D0
        chii_m(j)      = 0.D0
        etaphi_m(j)    = 0.D0
        etapar_m(j)    = 0.D0
        etaper_m(j)    = 0.D0
      enddo
c
c diagnostic arrays (not used)
c
c     do j=0,jmaxm
c       vstar_m(j)     = 0.
c       vexb_m(j)      = 0.
c       dnrate_m(j)    = 0.
c       dtnrate_m(j)   = 0.
c       dnfreq_m(j)    = 0.
c      do k=1,4
c       gamma_j(j,k)    = 0.
c       freq_j(j,k)     = 0.
c       phi_norm_j(j,k) = 0.
c      enddo
c     enddo
c
c     do j=0,jmaxm
c      do ik=1,20
c       gamma_k_j(ik,j) = 0.D0
c       freq_k_j(ik,j)  = 0.D0
c       chie_k_j(ik,j)  = 0.D0
c       chii_k_j(ik,j)  = 0.D0
c      enddo
c      vnewstare_m(j) = 0.
c      vnewstari_m(j) = 0.
c     enddo
c
c
************************************************************************
cmnt   profiles of quantities derived from model profiles of te,ti,ne
cmnt   revised derivedmodlocal to center between jm-jptf and jm+ptf
************************************************************************
 
c.......begin switches and settings......
 
c**********************************************************************
c     glf23 parameters
cxx      settings same as function glf23_v4_1_10
 
ck      write(6,*)  'jmm=',jmm,'going in'

      eigen_gf = leigen
c      nroot_gf=8     ! 8 for pure plasma, 12 for full impurity dynamics
      nroot_gf=nroot
      iflagin_gf(1)=0
      iflagin_gf(2)=1
      iflagin_gf(3)=1
      iflagin_gf(4)=0
      iflagin_gf(5)=3
 
      xparam_gf(1)=0.D0
      xparam_gf(2)=0
      xparam_gf(3)=.7D0
      xparam_gf(4)=0.D0
      xparam_gf(6)=0.D0
      xparam_gf(7)=1.D0
      xparam_gf(8)=0.D0
      xparam_gf(9)=1.D0
      xparam_gf(10)=0.D0
      xparam_gf(11)=0.D0
      xparam_gf(12)=0.D0
      xparam_gf(13)=0.2D0
      xparam_gf(14)=1.D0
      xparam_gf(15)=-0.1D0
      xparam_gf(16)=0.D0
      xparam_gf(17)=0.1D0
      xparam_gf(18)=.0D0
      xparam_gf(19)=0.D0
      xparam_gf(20)=0.D0
      xparam_gf(21)=0.D0
      xparam_gf(22)=0.D0
      xparam_gf(23)=1.D0
      xparam_gf(24)=0.D0
      xparam_gf(25)=0.D0
      xparam_gf(26)=0.D0
      xparam_gf(27)=0.D0
      xparam_gf(28)=0.D0
      xparam_gf(29)=0.D0
      xparam_gf(30)=0.D0
c
      xky0_gf= .2D0
      rms_theta_gf=zpi/3.D0
      park_gf  =0.7D0
      ghat_gf  =1.D0
      gchat_gf =1.D0
 
      adamp_gf=.50D0
      alpha_star_gf  =0.D0
      alpha_mode_gf=0.D0
      gamma_e_gf  =0.D0
ctemp
      gamma_e_gf  =-.000000000001D0
      xkdamp_gf     =0.D0
 
      alpha_p_gf=0.50D0
 
c   cbetae=1 is full electromagetic
      cbetae=1.D-6
c      cbetae=1.D0
c   full collisionality
      cxnu=1.D0
 
      cnorm_gf=100.D0
 
      ikymax_gf=10
      xkymin_gf=.02D0
      xkymax_gf=.5D0
 
c# non glf23 paramerter
 
       cmodel=1.D0
       xalpha=x_alpha
 
c      ialphastab=1
c      ineo=-2
 
c      iexch_m=1
c      iexp_exch=-1
 
c      i_dengrad=2
c      iexp_imp=1
c      igeo_m=3
 
c      irotstab=1
c      irot1=1
c      irot2=1
c
cxx      endf
 
c......begin important optional settings
 
c turn on self-consistant alpha-stabilization
c       ialphastab=1
 
c       turn on EXB shear stabilization
c       alpha_e_gf=1. full on ExB shear
        alpha_e_gf=alpha_e
 
c       turn on self consistant EXB stabilization
c       irotstab=1
 
c itport_pt(1)=1 plasma transport on; itport_pt(2:3)=1; electron and ion heat on
c itport_pt(4)=1 ;itport_pt(5)=0  transport vphi with neoclassical determining vtheta
c itport_pt(4)=1 ;itport_pt(5)=1  transport vphi and vtheta with fast time scale
c       neoclassical drag built into vphi and vtheta transport equations...
c       consult G.M. Staebler
 
c if only vphi_exp is available itport_pt(4)=0
 
c      grid-centering in computing EXB shear  span jm-jptf to jm+jpt
 
c      turn on high-k eta-e modes
        xparam_gf(10)=1.D0
c
c    relaxation turned off relx=0.
c    relaxation can be turned on for one call per step
       relx=0.D0
c
c settings for retuned GLF23 model
c
       if (iglf.eq.1) then      ! retuned model
         cnorm_gf=50.D0         ! ITG normalization (via GYRO runs)
         xparam_gf(10)=12.D0    ! ETG normalization (cnorm*xparam(10))
         iflagin_gf(5)=5        ! rms theta fit formula
         xparam_gf(13)=0.15     ! rms_theta q-dependence
         xparam_gf(16)=0.15     ! rms_theta shat dependence
         xparam_gf(17)=0.25     ! rms_theta shat dependence
         xparam_gf(19)=1.0      ! rms_theta alpha dependence
         adamp_gf=.70D0         ! radial mode damping exponent
         alpha_p_gf=0.35D0      ! parallel velocity shear fit
         park_gf=0.8D0          ! parallel ion motion fit
         bt_flag=1              ! use real geometry ExB shear
       endif
c
c.......end important optional settings
 
c.......end switches and settings......
 
c.......start setups.....................
c
 
c*********************************************************************************
c GEOMETRY FACTORS NEEDED FOR SETUP
c external
c      rho(jm)    :toroidal flux co-ordinate 0:50 grids 0->1
c      gradrhosq_exp(jm) : <|grad rho_phys |**2> toroidal flux= B0*rho_phys**2/2
c                 rho_phys=rho*arho_exp
c                 hence   gradrhosq_exp ->1 for a circle
c      gradrho_exp(jm)  : <|grad rho_phys |>
c internal
c      drhodr(jm)
c      drhodrrrho(jm)
c      geofac(jm)
c
c        geofac(j)=gradrho_exp(j)*(rho(j+1)-rho(j))*arho_exp
c     >   /(rmin_exp(j+1)-rmin_exp(j))/gradrhosq_exp(j)
c
c        drhodr(j)=(rho(j+1)-rho(j))*arho_exp/
c     >   (rmin_exp(j+1)-rmin_exp(j))
c
c        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/arho_exp/rho(j)
 
c  surface factor for power flow
c        sfactor(j)=2.*pi_m*arho_exp*rho(j)*h_exp(j)*2.*pi_m*rmajor_exp
c        h_exp(j-1)=hcap_d(j) hcap in ONETWO
 
c******************************************************************************

      if(jmm.gt.0) then
       jin=jmm
       jout=jmm
      endif
      if(jmm.eq.0) then
       jin=1
       jout=jmaxm-1
      endif
      do jm=jin,jout
 
c time dependent codes jshoot=0
c diffusion coefficients and gradients between jm+1 and jm
 
c outside to inside shooting codes jshoot=1 diffusion coefficient at jm
c gradient is forward jm to jm-1
c backward gradient is jm to jm+1
c shear is between forward and backward gradient.
c forward is implicit and backward is already updated
 
       tem=(te_m(jm+1-jshoot)+te_m(jm))/2.D0
       tim=(ti_m(jm+1-jshoot)+ti_m(jm))/2.D0
       nem=(ne_m(jm+1-jshoot)+ne_m(jm))/2.D0
       nim=(ni_m(jm+1-jshoot)+ni_m(jm))/2.D0
       nsm=(ns_m(jm+1-jshoot)+ns_m(jm))/2.D0
       zeffm=(zeff_exp(jm+1-jshoot)+zeff_exp(jm))/2.D0
 
       betae_m(jm) = 400.D0*nem*tem/(1.D5*bt_exp**2)
c      betai_m(jm) = 400.*nim*tim/(1.e5*bt_exp**2)
 
 
crew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
       vnewk3x=
     >   0.117D0*nem*tem**(-1.5D0)/(tim**0.5D0)*arho_exp*
     >   (amassgas_exp/2.D0)**0.5D0
crew   as used in gks multiply by 1/2 and takout any 1/2 factor in solfp
crew          vnewk3x=vnewk3x/2.
       xnu_m(jm) =vnewk3x/(2.D0*tem/tim)**0.5D0
crew  10/25/95 fixed zeff+1 factor: zeff col with ions;1 col with elecs.
       zeff_e=0.D0
       xnu_m(jm) = xnu_m(jm)*(zeff_exp(jm)+zeff_e)
 
c..Bateman, 14 March 2009, diagnostic printout

      IF ( iprint .gt. 0 ) THEN

      WRITE(*,*) ' vnewk3x=',vnewk3x
     & ,' nem=',nem
     & ,' tem=',tem
     & ,' tim=',tim

      WRITE(*,*) ' arho_exp=',arho_exp
     & ,' amassgas_exp=',amassgas_exp
     & ,' taui_gf=',taui_gf
     & ,' zeff_exp=',zeff_exp(jm)

      WRITE(*,*) ' jshoot = ', jshoot,' jm = ',jm
     & ,' te_m(jm)=', te_m(jm)
     & ,' te_m(jm+1-jshoot)=',te_m(jm+1-jshoot)

      ENDIF
 
c      vnewstare_m(jm)=zeff_exp(jm) *2.91e-6*nem*1.e13*15./
c    >  (tem*1.e3)**2*rmaj_exp(jm)*100.*q_exp(jm)
c    >     /(rmin_exp(jm)/rmaj_exp(jm)+epsilon)**1.5/4.19e7
c
c      vnewstari_m(jm)=4.78e-8*nem*1.e13*15./
c    >  (tim*1.e3)**2*rmaj_exp(jm)*100.*q_exp(jm)
c    >     /(rmin_exp(jm)/rmaj_exp(jm)+epsilon)**1.5/9.79e5
 
c
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
 
c  note: dependence of shear's  on zpxx_m(jm-1) and zpxx_m(jm+1)
c  and alpha(jm) depends on zpxx_m(jm+1)
 200  format(i2,2x,0p1f4.2,0p5f10.5)
c
      do j=jm-jptf,jm+jpt
c
c... some geometric factors
c
        geofac(j)=gradrho_exp(j)*(rho(j)-rho(j-1))*arho_exp
     >   /(rmin_exp(j)-rmin_exp(j-1)+epsilon)/gradrhosq_exp(j)
 
        drhodr(j)=(rho(j)-rho(j-1))*arho_exp/
     >   (rmin_exp(j)-rmin_exp(j-1)+epsilon)
 
        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/
     >   arho_exp/(rho(j)+epsilon)

c..Bateman, 20 March 2009, diagnostic printout

      IF ( iprint .gt. 0 ) THEN

        WRITE(*,*)
     &   ' j   z_geo_factor  gradRho  drhodr zarho_exp gradRhoSq'
     &   , ' rmin'
        WRITE(*,*)
     &   j, geofac(j), gradrho_exp(j), drhodr(j), arho_exp
     &    , gradrhosq_exp(j), rmin_exp(j)

      ENDIF
c
c... local rate unit
c
       csda_m(j)=9.79D5*(te_m(j)*1.D3)**.5D0/
     >    (arho_exp*100.D0)/amassgas_exp**0.5D0
c
c... local rho_star
c Note: use effective B-field if bt_flag > 0
c
       if (bt_flag .gt. 0) then
         bteff_exp(j)=bt_exp*rho(j)*arho_exp/
     >         rmin_exp(j)*drhodr(j)
         rhosda_m(j)=((1.02D2*(te_m(j)*1.D3)**.5D0)/bteff_exp(j)
     >         /1.D4)*amassgas_exp**.5D0/(arho_exp*100.D0)
       else
         rhosda_m(j)=((1.02D2*(te_m(j)*1.D3)**.5D0)/bt_exp/1.D4)
     >         *amassgas_exp**.5D0/(arho_exp*100.D0)
       endif
      enddo
c
c   local gyrobohm unit of diffusion
c
       cgyrobohm_m(jm)=1.D-4*
     >  9.79D5*(tem*1.D3)**.5D0/(arho_exp*100.D0)
     >  *(1.02D2*(tem*1.D3)**.5D0/bt_exp/1.D4)**2*amassgas_exp**.5D0
c 
      do j=jm-jptf, jm+jpt
        drho=rho(j-1)-rho(j)+epsilon
        zpte_m(j)=-(dlog(te_m(j-1))-dlog(te_m(j)))/drho
        zpti_m(j)=-(dlog(ti_m(j-1))-dlog(ti_m(j)))/drho
        zpne_m(j)=-(dlog(ne_m(j-1))-dlog(ne_m(j)))/drho
        zpni_m(j)=-(dlog(ni_m(j-1))-dlog(ni_m(j)))/drho
      enddo
 
        zpmte=zpte_m(jm+1-jshoot)
        zpmti=zpti_m(jm+1-jshoot)
        zpmne=zpne_m(jm+1-jshoot)
        zpmni=zpni_m(jm+1-jshoot)
c
c... check on zero norm gradients:
c
        if (DABS(zpmti).lt.zeps) zpmti=zeps
        if (DABS(zpmte).lt.zeps) zpmte=zeps
        if (DABS(zpmne).lt.zeps) zpmne=zeps
        if (DABS(zpmni).lt.zeps) zpmni=zeps
c
        if(i_grad.eq.1) then
         zpmte=zpte_in
         zpmti=zpti_in
         zpmne=zpne_in
         zpmni=zpni_in
        endif
 
c MHD alpha parameter
 
         alpha_m(jm)=drhodr(jm)*
     >    q_exp(jm)**2*rmaj_exp(jm)/arho_exp*
     >    betae_m(jm)*((tim/tem*nim/nem)*
     >   (zpmni+zpmti)
     >    +zpmne+zpmte)
 
c vstarp_m is diamagnetic part of egamma_m (doppler shear rate)
c vstar_sign is negative for negative vstar_i. Thus for co-injection or positive
c angrot toroidal rotation cancels the diamgnetic rotation
 
       vstar_sign=-1.D0
 
        j=jm
          rho(j-jptf)=rho(j-jptf)+epsilon
          rho(j+jpt)=rho(j+jpt)+epsilon
          egeo_local=1.D0
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
          rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
 
        vstarp_m(jm)=(
     > +egeo_local*vstar_sign*
     > (rho(j+jpt)*rdrho_local_p1+rho(j-jptf)*rdrho_local)/2.D0*(
     >(ti_m(j+jpt)/te_m(j+jpt))*csda_m(j+jpt)*
     > (zpti_m(j+jpt)+zpni_m(j+jpt))
     >  *pgeo_local*rhosda_m(j+jpt)/rho(j+jpt)/rdrho_local_p1
     >-(ti_m(j-jptf)/te_m(j-jptf))*csda_m(j-jptf)*
     > (zpti_m(j-jptf)+zpni_m(j-jptf))
     >  *pgeo_local*rhosda_m(j-jptf)/rho(j-jptf)/rdrho_local
     >  )/(rho(j+jpt)-rho(j-jptf)+epsilon)/csda_m(j)
     >  )
c
       do jj=1,2
 
        if(jj.eq.1) j=jm-jptf
 
        if(jj.eq.2) j=jm+jpt

c banana regime ie collisionless limit formulas
 
        fc=1-1.46D0*(rmin_exp(j)/rmaj_exp(j))**0.5D0+
     >      0.46D0*(rmin_exp(j)/rmaj_exp(j))**1.5D0
        akappa1=0.8839D0*fc/(0.3477D0+0.4058D0*fc)
        alpha_neo=-akappa1+1.D0
        alpha_neo_hold=alpha_neo
 
cxc angrotp is plasma rotation
cxc angrot is impurity rotation from experiment
cxc if angrotp_exp is not suppied must insert
cx
cx        akappa2=(1.-fc)/(1.+1.1671*fc)
cxc trace impurity limit to go from impurity rotation angrot_exp to
cxc plasma rotation angrotp_exp
cx        angrotp_exp(j)=angrot_exp(j)+
cx     >     akappa2*3./2.*csda_exp(j)*zpti_exp(j)*(ti_exp(j)/te_exp(j))
cx     >     /rho(j)*q_exp(j)*rhosda_exp(j)*pgeo_local/rdrho_local
cx        if(angrot_exp(j).eq.0.) angrotp_exp(j)=0.
cx        angrotp_exp(j)=corot*angrotp_exp(j)
 
          egeo_local=1.D0
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j)/arho_exp/(rho(j)+epsilon)
 
        ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     > -rho(j)*rdrho_local*
     >  arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)
 
        vpar(j)=rmajor_exp*angrotp_exp(j)-vstar_sign*
     > (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >((alpha_neo-1.D0)*zpti_m(j))*rho(j)*rdrho_local*
     > arho_exp/rmajor_exp/q_exp(j)
 
c        vmode(j)=anfreq_m(j)/(ky_j(j)+epsilon)*
c     >    csda_m(j)*arho_exp*rhosda_m(j)
 
       if(itport_pt(4).eq.0.and.itport_pt(5).eq.0) then
         vphi_m(j)=rmajor_exp*angrotp_exp(j)
 
         vpar_m(j)=vpar(j)
 
         vper_m(j)=ve(j)
     >  +(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
       if(abs(itport_pt(4)).eq.1.and.itport_pt(5).eq.1) then
         vpar(j)=vpar_m(j)
       endif
 
       if(abs(itport_pt(4)).eq.1.and.itport_pt(5).eq.0) then
c this option vpar is vphi vexb from neo+vphi
         ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >  (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     >  -rho(j)*rdrho_local*
     >   arho_exp/rmajor_exp/q_exp(j)*vphi_m(j)
 
         vpar(j)=vphi_m(j)-vstar_sign*
     >   (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >   ((alpha_neo-1.D0)*zpti_m(j))*rho(j)*rdrho_local*
     >   arho_exp/rmajor_exp/q_exp(j)
 
         vpar_m(j)=vpar(j)
 
         vper_m(j)=ve(j)
     >   +(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
       if(itport_pt(5).eq.1) then
c this option vexb from vper and vpar with neo dampng built into
c vpar and vper transport equations
         ve(j)=vper_m(j)
     >   -(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
c       vexb_m(j)=ve(j)
c       vmode_m(j)=vmode(j)
c       vstar_m(j)=(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
c    > (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
 
      enddo
c 
c compute shears from outside minus inside
c 
        j=jm
        rho(j-jptf)=rho(j-jptf)+epsilon
        rho(j+jpt)=rho(j+jpt)+epsilon
        egeo_local=1.D0
        pgeo_local=drhodr(j)
        rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
        rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
c
        egamma_m(jm)=relx*egamma_m(jm)+(1.D0-relx)*(
     >  egeo_local*drhodrrrho(j)*
     >  (rho(j-jptf)+rho(j+jpt))/(q_exp(j-jptf)+q_exp(j+jpt))*
     >  (ve(j+jpt)*q_exp(j+jpt)/rho(j+jpt)/rdrho_local_p1-
     >  ve(j-jptf)*q_exp(j-jptf)/rho(j-jptf)/rdrho_local)/
     >  (rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
     >  )
c
c       write(*,*) jm, rho(jm), egamma_m(jm), ' egamma'
c
c       gamma_mode_m(jm)=relx*gamma_mode_m(jm)+(1.-relx)*(
c    >  egeo_local*drhodrrrho(j)*
c    >  (rho(j-jptf)+rho(j+jpt))/2.*
c    >  (vmode(j+jpt)/rho(j+jpt)/rdrho_local_p1-
c    >  vmode(j-jptf)/rho(j-jptf)/rdrho_local)/
c    >  (rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
c    >  )
 
        gamma_p_m(jm)=relx*gamma_p_m(jm)+(1.D0-relx)*(
     >   -drhodr(j)*
     >   (vpar(j+jpt)-vpar(j-jptf))
     >   /(rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
     >  )
 
       if (jm.eq.1.and.gamma_p_m(jm).gt.10)
     >    gamma_p_m(jm)=10.D0
        alpha_neo=alpha_neo_hold
 
c.......end   setups...........
 
c***********************************************************************
cvv      subroutine model
************************************************************************
 
c   units:
c     diffusion (m**2/sec) note: computed in derived modprofiles
c     density (10**13 cm**-3 or 10**19 m**-3)
c     arho_exp and rmajor_exp (m)
c     power (MW)
c     flow (MW/kev=kA)
c
c   kev/sec per MW
c   kevdsecpmw=1.6022e-19*1.0e3*1.e-6
c
c       cgyrobohm_m(jm)=1.e-4*
c    >  9.79e5*(tem*1.e3)**.5/(arho_exp*100.)
c    >  *(1.02e2*(tem*1.e3)**.5/bt_exp/1.e4)**2*(amassgas_exp)**.5
c
c   sfactor(j)=
c     >      2.*pi_m*arho_exp*rho(j)*h_exp(j)*2.*pi_m*rmajor_exp
c   units of meter**2
c
c   9000 format (1i6,7e10.3)
 
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cmnt                         the model
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cmnt
cmnt supply model for chietem,chietim,chienem
cmnt                  chiitem,chiitim,chiinem
cmnt                  difftem,difftim,diffnem
cmnt
cmnt    chi_s and diff_s must be in meters**2/sec units
cmnt     and recall chi_s refer to total energy flow
cmnt
cmnt    if the model chi_s refer to "heat conduction" flow
cmnt    then a convection term xconv*3./2.*t_m*flow_exp is added.
cmnt    normally input xconv=0. otherwise xconv=1. or 5./3.
cmnt
cmnt    it is also possible to build convection into the model
cmnt    with "aconv".  aconv and xconv should not be double counted.
cmnt
cmnt note: can use diagonal forms with off diagonal dependence on
cmnt zpmte,zpmti,zpmne intrinsic to the diagonal elements as in sample
cmnt normall models are written in diagonal for with dependence on off
cmnt diagonal gradients implicit
cmnt
cmnt note: when flow is large anomalous e-i exchange should be added
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cxx      if (imodel.eq.8) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 2DGLF quasilinear model  GLF23
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c  test effect of canonical density gradients
c note pzmn_sol only in zpmne_q and zpmni_q which go into
c glf2d drivers rlne_gf and rlni_gf
 
c       zpmne_q=(1.+xparam_gf(20))*zpmne*pzmn_sol(jm)-xparam_gf(19)*
c     > (alog(q_exp(jm-1))-alog(q_exp(jm)))/(rho(jm-1)-rho(jm))
c       zpmni_q=(1.+xparam_gf(20))*zpmni*pzmn_sol(jm)-xparam_gf(19)*
c     > (alog(q_exp(jm-1))-alog(q_exp(jm)))/(rho(jm-1)-rho(jm))
 
       zpmne_q= zpmne
       zpmni_q= zpmni
 
       rmaj_gf=rmaj_exp(jm)/arho_exp
       rmin_gf=rmin_exp(jm)/arho_exp
       q_gf=q_exp(jm)
       betae_gf=dmax1(cbetae*betae_m(jm), 1.D-6)
       shat_gf=shat_exp(jm)*drhodrrrho(jm)
 
       if(xalpha.lt.0.) alpha_gf=-xalpha*alpha_m(jm)
       if(xalpha.gt.0.) alpha_gf=xalpha*alpha_exp(jm)
       if(alpha_gf.gt.4.D0) alpha_gf=4.D0
 
       elong_gf=elong_exp(jm)
 
       xnu_gf=cxnu*xnu_m(jm)
       taui_gf=tim/tem
       amassgas_gf=amassgas_exp
 
       apwt_gf=1.D0
c impurity dynamics not turned on by default 
c and simple dilution included (idengrad=2, dil_gf=1-nim/nem)
c to turn on impurity dynamics need to change number of roots
c supply zimp_exp, amassimp_exp, and fractional density weights
c apwt_gf and aiwt_gf
       dil_gf=0.D0
       aiwt_gf=0.D0
c      zimp_gf=6.D0
c      amassimp_gf=12.D0
       zimp_gf=zimp_exp
       amassimp_gf=amassimp_exp
       rlnimp_gf=1.D0
       zpmnimp=1.D0
       if (idengrad.eq.2) dil_gf=1.D0-nim/nem
       if (idengrad.eq.3) then
         apwt_gf=nim/nem
         aiwt_jp1=(zeffm*nem-ni_m(jm+1)
     >            -ns_m(jm+1))/(zimp_gf**2*nem)
         xnimp_jp1=aiwt_jp1*ne_m(jm+1)
         aiwt_gf=(zeffm*nem
     >           -nim-ns_m(jm))/(zimp_gf**2*ne_m(jm))
         xnimp=aiwt_gf*ne_m(jm)
         zpmnimp=-(dlog(xnimp_jp1)-dlog(xnimp))/
     >           (rho(jm+1)-rho(jm))
         rlnimp_gf=zpmnimp*elong_exp(jm)**0.5
       endif
c       write(*,66) rho(jm),nem,nim,rlnimp_gf,zpmnimp,amassimp_gf,
c     >             xnimp,zimp_gf,dil_gf
c 66    format(2x,0p1f4.2,0p9f10.5)
 
         rlte_gf=zpmte*drhodr(jm)
         rlti_gf=zpmti*drhodr(jm)
         rlne_gf=zpmne_q*drhodr(jm)
         rlni_gf=zpmni_q*drhodr(jm)
         rlnimp_gf=zpmnimp*drhodr(jm)
c        write(*,200) jm, rho(jm), rlti_gf, rlte_gf, rlne_gf, rlni_gf
 
        gamma_star_gf=vstarp_m(jm)
        gamma_e_gf=egamma_m(jm)
        gamma_p_gf=gamma_p_m(jm)
        if(itport_pt(4).eq.-1) then
          gamma_e_gf=egamma_exp(jm)
          gamma_p_gf=gamma_p_exp(jm)
        endif
c..jek 8/15/00
        if (irotstab.eq.0) then
          gamma_e_gf=egamma_exp(jm)
          gamma_p_gf=gamma_p_exp(jm)
          egamma_m(jm)=egamma_exp(jm)
        endif
c       gamma_mode_gf=gamma_mode_m(jm)
        gamma_mode_gf=0.0D0
 
        if(i_delay.ne.0) then
c   i_delay should be negative on any intermediate step
 
         if(i_delay.gt.1) then
          do ii=1,i_delay-1
           egamma_d(jm,ii)=egamma_d(jm,ii+1)
          enddo
 
          egamma_d(jm,i_delay)=egamma_m(jm)
         endif
          gamma_e_gf=egamma_d(jm,1)
        endif

c..Bateman 14 March 2009, diagnostic printout

      IF ( iprint .gt. 0 ) THEN

        IF ( initial .lt. 1 ) THEN

          initial = initial + 1

          WRITE (*,*) ' Output from sglf.f'
          WRITE (*,*) ' iflagin_gf = ', (iflagin_gf(j), j=1,30)

          WRITE (*,*) ' nroot_gf=',nroot_gf,' jeigen=',jeigen
     &   ,' lprint_gf=',lprint_gf,' ikymax_gf=',ikymax_gf
          WRITE (*,*) 'eigen_gf=',eigen_gf,' i_err=', i_err
     &   ,' first_order_gf=',first_order_gf,' ipert_gf=',ipert_gf
          WRITE (*,*) ' xparam_gf(30)'
          DO j=1,30
            WRITE(*,*) j, xparam_gf(j)
          ENDDO
          WRITE (*,*)

        ENDIF

        WRITE (*,*) ' Output from sglf.f'
        WRITE (*,*)
     & ' xky0_gf, rms_theta_gf, rlti_gf, rlte_gf, rlne_gf, rlni_gf'
        WRITE (*,'(1X,7ES11.3)')
     &  xky0_gf,rms_theta_gf,rlti_gf,rlte_gf,rlne_gf,rlni_gf
        WRITE (*,*)
     & ' rlnimp_gf, dil_gf, apwt_gf, aiwt_gf, taui_gf, rmin_gf, rmaj_gf'
        WRITE (*,*)
     &  rlnimp_gf,dil_gf,apwt_gf,aiwt_gf,taui_gf,rmin_gf,rmaj_gf
        WRITE (*,*)
     & ' q_gf, xnu_gf, betae_gf, shat_gf, alpha_gf, elong_gf, xwell_gf'
        WRITE (*,*)
     &  q_gf,xnu_gf,betae_gf,shat_gf,alpha_gf,elong_gf,xwell_gf
        WRITE (*,*)
     & ' park_gf, ghat_gf, gchat_gf, adamp_gf'
     & ,' alpha_star_gf, gamma_star_gf'
        WRITE (*,*)
     &  park_gf,ghat_gf,gchat_gf,adamp_gf,alpha_star_gf,gamma_star_gf
        WRITE (*,*)
     & ' alpha_e_gf, gamma_e_gf, alpha_mode_gf, gamma_mode_gf'
     & ,' alpha_p_gf, gamma_p_gf, xkdamp_gf'
        WRITE (*,*)
     &  alpha_e_gf,gamma_e_gf,alpha_mode_gf,gamma_mode_gf
     & , alpha_p_gf,gamma_p_gf,xkdamp_gf
        WRITE (*,*)
     & ' xkyf_gf, cnorm_gf, xkymin_gf, xkymax_gf'
        WRITE (*,*)
     &  xkyf_gf,cnorm_gf,xkymin_gf,xkymax_gf
        WRITE (*,*)
     & ' amassgas_gf, amassimp_gf, zimp_gf'
        WRITE (*,*)
     &  amassgas_gf,amassimp_gf,zimp_gf

        WRITE (*,*)

      ENDIF
 
c.......THE  CALL TO GLF23
 
        call glf2d(iglf)
 
c.......POST CALL TO GLF23
 
c...diagnostic arrays (not presently used)
c
c      ky_j(jm)=xkyf_gf
c      gamma_j(jm,1)=gamma_gf(1)
c      gamma_j(jm,2)=gamma_gf(2)
c      gamma_j(jm,3)=gamma_gf(3)
c      gamma_j(jm,4)=gamma_gf(4)
c
c      freq_j(jm,1)=freq_gf(1)
c      freq_j(jm,2)=freq_gf(2)
c      freq_j(jm,3)=freq_gf(3)
c      freq_j(jm,4)=freq_gf(4)
c
c      phi_norm_j(jm,1)=phi_norm_gf(1)
c      phi_norm_j(jm,2)=phi_norm_gf(2)
c      phi_norm_j(jm,3)=phi_norm_gf(3)
c      phi_norm_j(jm,4)=phi_norm_gf(4)
c
c      do ik=1,ikymax_gf
c       gamma_k_j(ik,jm)= gamma_k_gf(1,ik)
c       freq_k_j(ik,jm) = freq_k_gf(1,ik)
c       chie_k_j(ik,jm) = chie_k_gf(ik)
c       chii_k_j(ik,jm) = chii_k_gf(ik)
c      enddo
 
       anrate_m(jm)=gamma_gf(1)
       anrate2_m(jm)=gamma_gf(2)
c      dnrate_m(jm)=0.
c      dtnrate_m(jm)=0.
       anfreq_m(jm)=freq_gf(1)
       anfreq2_m(j)=freq_gf(2)
c      dnfreq_m(jm)=0.
c       xkymax_m(jm)=xky_gf(1)
c       xkymax2_m(jm)=xky_gf(2)

       gfac=geofac(jm)
c
c exch_m in MW/m**3
c   kev/sec per MW
ck       kevdsecpmw=1.6022e-19*1.0e3*1.e-6
 
ck      exch_m(jm)=1.e19*
ck     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2*exch_gf*cmodel
      exch_m(jm)=1.D19*1.6022D-19*1.0D3*1.D-6*
     > nem*tem*csda_m(jm)*rhosda_m(jm)**2*exch_gf*cmodel
 
      exchm=exch_m(jm)
 
c exch_m is directly related to the flow
c for a single mode branch exch_gf=-(-freq_gf(1)/xkyf_gf)*diff_gf*rln_gf.
c we can  not expect to find exch_m without knowing flow_exp as input.
c and solving self consistant flow eq. flown=flow_exp for density
c density profile.
c
c however, knowing freq_gf(1) from the gf model we can compute exch_exp
c from flow_exp using
c       flowm=kevdsecpmw*1.*nem*1.e19/arho_exp*gradrhosq_exp(jm)*
c     >       sfactor(jm)*(difftem*zpmte+difftim*zpmti+diffnem*zpmne)
c we have:
 
c       diffgb_local=flow_exp(jm)/
c     > (kevdsecpmw*1.*nem*1.e19/arho_exp*gradrhosq_exp(jm)*sfactor(jm)*
c     > zpmne_q)/cgyrobohm_m(jm)
 
c       exchgb_local=-(-freq_gf(1)/xkyf_gf)*diffgb_local*rlni_gf
 
c       exch_exp(jm)=1.e19*
c     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2*exchgb_local
 
c       exch_exp(jm)=flow_exp(jm)*tem*(-1.)*(-freq_gf(1)/xkyf_gf)*
c     > sqrt(elong_exp(jm))*arho_exp/gradrhosq_exp(jm)/sfactor(jm)
c     >/arho_exp(jm)**2
 
 
c   note electron(ion) wave freq > 0(<0) cool(heat) electrons
c (-1) denotes electron to ion
 
c to emphasize, we can not know exch_exp better than we know flow_exp
 
c chietem, chiitim, diffen in m**2/sec
c ITER definition of "chi" assumes will be proceeded with gradrhosq factor
 
      chietem=cmodel*gfac*chie_gf*cgyrobohm_m(jm)
      chiitim=cmodel*gfac*chii_gf*cgyrobohm_m(jm)
      diffnem=cmodel*gfac*diff_gf*cgyrobohm_m(jm)
 
      etaphim=cmodel*gfac*eta_phi_gf*cgyrobohm_m(jm)
      etaparm=cmodel*gfac*eta_par_gf*cgyrobohm_m(jm)
      etaperm=cmodel*gfac*eta_per_gf*cgyrobohm_m(jm)
 
cxx endif
c
      if ( itport_pt(1) .eq. 0 ) diffnem = 0.D0
      if ( itport_pt(2) .eq. 0 ) chietem = 0.D0
      if ( itport_pt(3) .eq. 0 ) chiitim = 0.D0
      if ( (itport_pt(4) .eq. 0) .and. (itport_pt(5) .eq. 0) ) then
         etaphim=0.D0
         etaparm=0.D0
         etaperm=0.D0
      endif
c
      diff_m(jm)=diffnem
      chie_m(jm)=chietem
      chii_m(jm)=chiitim
 
      etaphi_m(jm)=etaphim
      etapar_m(jm)=etaparm
      etaper_m(jm)=etaperm
c
c     write(6,*) jm,rho(jm),zpmte, zpmti, zpmne, zpmni

!..Bateman, 14 March 2009, diagnostic printout

      IF ( iprint .gt. 0 ) THEN

        WRITE (*,*)
     & ' cgyrobohm_m(jm), chii_gf, chie_gf, diff_gf, diff_im_gf'
        WRITE (*,*)
     &   cgyrobohm_m(jm), chii_gf, chie_gf, diff_gf, diff_im_gf

        WRITE (*,*)
     & ' eta_phi_gf, eta_par_gf, eta_per_gf, chie_e_gf, exch_gf'
        WRITE (*,*)
     &   eta_phi_gf, eta_par_gf, eta_per_gf, chie_e_gf, exch_gf

        WRITE (*,*)

      ENDIF
 
        enddo
 
c      if(jmm.eq.0)  write(6,*) 'jmm=', jmm,'going out'
c
cdmc        write(6,7801) diffnem,chietem,chiitim,etaphim,etaparm,etaperm,
cdmc     >     exchm
cdmc 7801   format(//' OUTPUTS...'/
cdmc     >'  diffnem = ',1pe17.10,' chietem = ',1pe17.10/
cdmc     >'  chiitim = ',1pe17.10,' etaphim = ',1pe17.10/
cdmc     >'  etaparm = ',1pe17.10,' etaperm = ',1pe17.10/
cdmc     >'  exchm = ',1pe17.10)
c
cdmc        call echo('diff_m output',diff_m,jmaxm)
cdmc        call echo('chie_m output',chie_m,jmaxm)
cdmc        call echo('chii_m output',chii_m,jmaxm)
cdmc        call echo('etaphi_m output',etaphi_m,jmaxm)
cdmc        call echo('etapar_m output',etapar_m,jmaxm)
cdmc        call echo('etaper_m output',etaper_m,jmaxm)
cdmc        call echo('exch_m output',exch_m,jmaxm)
c
cdmc        call echo('egamma_m output',egamma_m,jmaxm)
cdmc        call echo('egamma_d output',egamma_d,jmaxm)
cdmc        call echo('gamma_p_m output',gamma_p_m,jmaxm)
cdmc        call echo('anrate_m output',anrate_m,jmaxm)
cdmc        call echo('anrate2_m output',anrate2_m,jmaxm)
cdmc        call echo('anfreq_m output',anfreq_m,jmaxm)
cdmc        call echo('anfreq2_m output',anfreq2_m,jmaxm)
c
       return
       end
c
cdmc      subroutine echo(lbl,arr,jmaxm)
cdmc      character*(*) lbl
cdmc      real*8 arr(0:jmaxm)
c
cdmc      write(6,1001) lbl
cdmc 1001 format(/5x,a)
cdmc      write(6,1002) (arr(j),j=0,jmaxm)
cdmc 1002 format(4(1x,1pe17.10))
c
cdmc      return
cdmc      end
!dmc -- real*8 version generated using `fgtok'.
c#include "f77_dcomplx.h"
!
      SUBROUTINE R8TOMSQZ(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL)
C-----------------------------------------------------------------------
C TOMSQZ  written by P. Strand 27-apr-98,       elfps@elmagn.chalmers.se
C-----------------------------------------------------------------------
C CQZHES, CQZVEC, and CQZVAL: Fortran subroutines implementing the QZ
C algorithm for solving the generalized eigenvalue problem for complex
C matrices. (See B.S.C Garbow, ACM TOMS 4 (1978) pp. 404-410.).
C-----------------------------------------------------------------------
C
C ON INPUT THE GENERALIZED EIGENVALUE PROBLEM IS DEFINED THROUGH THE
C COMPLEX MATRICES
C
C       A = cmplx (AR, AI)  AND   B = cmplx (BR, BI)
C
C WHERE LEADING DIMENSION N IS AS DEFINED IN THE CALLING ROUTINE AND
C WHERE NA IS THE ROW  RANGE IN THE CURRENT PROBLEM. THE EIGENVALUE
C PROBLEM IS THEN DEFINED THROUGH
C
C       A x = w B x
C
C WHERE  THE COMPLEX EIGENVECTORS
C
C       x = cmplx (ZVR, ZVI)
C
C TOGETHER WITH THE COMPLEX EIGENVALUE
C
C        w = cmplx(alfr, alfi)/beta
C
C IS OUTPUT FROM THE ROUTINE
C
C IFAIL WILL BE NONZERO IF CONVERGENCE HAS NOT BEEN REACH WITHIN 50
C ITERATIONS
C-----------------------------------------------------------------------
C DECLARATIONS FOR INPUT VARIABLES
C-----------------------------------------------------------------------
 
      IMPLICIT NONE
      INTEGER N, NA
      double precision AR(N,NA),AI(N,NA),BR(N,NA),BI(N,NA)
 
C-----------------------------------------------------------------------
C DECALRATIONS FOR OUTPUT VARIABLES
C-----------------------------------------------------------------------
 
      double precision ALFR(N),ALFI(N),BETA(N)
      double precision ZVR(N,NA), ZVI(N,NA)
      INTEGER IFAIL
 
C-----------------------------------------------------------------------
C LOCAL VARIABLES
C-----------------------------------------------------------------------
 
      LOGICAL WANTX
      double precision EPS1
      double precision ZERO,ZONE
 
C-----------------------------------------------------------------------
C START OF ACTUAL CODING
C-----------------------------------------------------------------------
 
      WANTX = .TRUE.
      EPS1  = -0.0D0
 
      CALL R8CQZHES(N,NA,AR,AI,BR,BI,WANTX,ZVR,ZVI)
      CALL R8CQZVAL(N,NA,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,WANTX,
     &            ZVR,ZVI,IFAIL)
      CALL R8CQZVEC(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI)
      RETURN
      END
 
C     ------------------------------------------------------------------
C
      SUBROUTINE R8CQZHES(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI)
C
      IMPLICIT NONE
      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1
      double precision AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),
     &       ZR(NM,N),ZI(NM,N)
      double precision R,S,T,TI,U1,U2,XI,XR,YI,YR,RHO,U1I
      LOGICAL MATZ
      double precision ZERO,ZONE
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
C     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
C     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
C     CQZVAL  AND POSSIBLY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
C          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
C          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
C          OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z **********
 
      ZERO = 0.0D0
      ZONE = 1.0D0
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 I = 1, N
C
         DO 2 J = 1, N
            ZR(I,J) = 0.0D0
            ZI(I,J) = 0.0D0
    2    CONTINUE
C
         ZR(I,I) = 1.0D0
    3 CONTINUE
C     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
C                TEMPORARILY REAL DIAGONAL ELEMENTS **********
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0.0D0
C
         DO 20 I = L, N
            S = S + ABS(BR(I,L)) + ABS(BI(I,L))
   20    CONTINUE
C
         IF (S .EQ. ZERO) GO TO 100
         RHO = 0.0D0
C
         DO 25 I = L, N
            BR(I,L) = BR(I,L) / S
            BI(I,L) = BI(I,L) / S
            RHO = RHO + BR(I,L)**2 + BI(I,L)**2
   25    CONTINUE
C
         R = SQRT(RHO)
         XR = abs(CMPLX(BR(L,L),BI(L,L)))
         IF (XR .EQ. ZERO) GO TO 27
         RHO = RHO + XR * R
         U1 = -BR(L,L) / XR
         U1I = -BI(L,L) / XR
         YR = R / XR + 1.0D0
         BR(L,L) = YR * BR(L,L)
         BI(L,L) = YR * BI(L,L)
         GO TO 28
C
   27    BR(L,L) = R
         U1 = -1.0D0
         U1I = 0.0D0
C
   28    DO 50 J = L1, N
            T = 0.0D0
            TI = 0.0D0
C
            DO 30 I = L, N
               T = T + BR(I,L) * BR(I,J) + BI(I,L) * BI(I,J)
               TI = TI + BR(I,L) * BI(I,J) - BI(I,L) * BR(I,J)
   30       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 40 I = L, N
               BR(I,J) = BR(I,J) - T * BR(I,L) + TI * BI(I,L)
               BI(I,J) = BI(I,J) - T * BI(I,L) - TI * BR(I,L)
   40       CONTINUE
C
            XI = U1 * BI(L,J) - U1I * BR(L,J)
            BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
            BI(L,J) = XI
   50    CONTINUE
C
         DO 80 J = 1, N
            T = 0.0D0
            TI = 0.0D0
C
            DO 60 I = L, N
               T = T + BR(I,L) * AR(I,J) + BI(I,L) * AI(I,J)
               TI = TI + BR(I,L) * AI(I,J) - BI(I,L) * AR(I,J)
   60       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 70 I = L, N
               AR(I,J) = AR(I,J) - T * BR(I,L) + TI * BI(I,L)
               AI(I,J) = AI(I,J) - T * BI(I,L) - TI * BR(I,L)
   70       CONTINUE
C
            XI = U1 * AI(L,J) - U1I * AR(L,J)
            AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
            AI(L,J) = XI
   80    CONTINUE
C
         BR(L,L) = R * S
         BI(L,L) = 0.0D0
C
         DO 90 I = L1, N
            BR(I,L) = 0.0D0
            BI(I,L) = 0.0D0
   90    CONTINUE
C
  100 CONTINUE
C     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
C                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
      DO 160 K = 1, NM1
         K1 = K + 1
C     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
         IF (AI(N,K) .EQ. ZERO) GO TO 105
         R = abs(CMPLX(AR(N,K),AI(N,K)))
         U1 = AR(N,K) / R
         U1I = AI(N,K) / R
         AR(N,K) = R
         AI(N,K) = 0.0D0
C
         DO 103 J = K1, N
            XI = U1 * AI(N,J) - U1I * AR(N,J)
            AR(N,J) = U1 * AR(N,J) + U1I * AI(N,J)
            AI(N,J) = XI
  103    CONTINUE
C
         XI = U1 * BI(N,N) - U1I * BR(N,N)
         BR(N,N) = U1 * BR(N,N) + U1I * BI(N,N)
         BI(N,N) = XI
  105    IF (K .EQ. NM1) GO TO 170
         NK1 = NM1 - K
C     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     ********** ZERO A(L+1,K) **********
            S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K)
            IF (S .EQ. ZERO) GO TO 150
            U1 = AR(L,K) / S
            U1I = AI(L,K) / S
            U2 = AR(L1,K) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            AR(L,K) = R * S
            AI(L,K) = 0.0D0
            AR(L1,K) = 0.0D0
C
            DO 110 J = K1, N
               XR = AR(L,J)
               XI = AI(L,J)
               YR = AR(L1,J)
               YI = AI(L1,J)
               AR(L,J) = U1 * XR + U1I * XI + U2 * YR
               AI(L,J) = U1 * XI - U1I * XR + U2 * YI
               AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110       CONTINUE
C
            XR = BR(L,L)
            BR(L,L) = U1 * XR
            BI(L,L) = -U1I * XR
            BR(L1,L) = -U2 * XR
C
            DO 120 J = L1, N
               XR = BR(L,J)
               XI = BI(L,J)
               YR = BR(L1,J)
               YI = BI(L1,J)
               BR(L,J) = U1 * XR + U1I * XI + U2 * YR
               BI(L,J) = U1 * XI - U1I * XR + U2 * YI
               BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  120       CONTINUE
C     ********** ZERO B(L+1,L) **********
            S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) + ABS(BR(L1,L))
            IF (S .EQ. ZERO) GO TO 150
            U1 = BR(L1,L1) / S
            U1I = BI(L1,L1) / S
            U2 = BR(L1,L) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            BR(L1,L1) = R * S
            BI(L1,L1) = 0.0D0
            BR(L1,L) = 0.0D0
C
            DO 130 I = 1, L
               XR = BR(I,L1)
               XI = BI(I,L1)
               YR = BR(I,L)
               YI = BI(I,L)
               BR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               BI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               BR(I,L) = U1 * YR - U1I * YI - U2 * XR
               BI(I,L) = U1 * YI + U1I * YR - U2 * XI
  130       CONTINUE
C
            DO 140 I = 1, N
               XR = AR(I,L1)
               XI = AI(I,L1)
               YR = AR(I,L)
               YI = AI(I,L)
               AR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               AI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               AR(I,L) = U1 * YR - U1I * YI - U2 * XR
               AI(I,L) = U1 * YI + U1I * YR - U2 * XI
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = 1, N
               XR = ZR(I,L1)
               XI = ZI(I,L1)
               YR = ZR(I,L)
               YI = ZI(I,L)
               ZR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               ZI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               ZR(I,L) = U1 * YR - U1I * YI - U2 * XR
               ZI(I,L) = U1 * YI + U1I * YR - U2 * XI
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
  170 RETURN
C     ********** LAST CARD OF CQZHES **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE R8CQZVAL(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,
     X                                       MATZ,ZR,ZI,IERR)
C
      IMPLICIT NONE
      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,
     X        ENM2,IERR,LOR1,ENORN
      double precision AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),
     X       ALFR(N),ALFI(N),
     X       BETA(N),ZR(NM,N),ZI(NM,N)
      double precision R,S,A1,A2,EP,SH,U1,U2,XI,XR,YI,YR,
     X       ANI,A1I,A33,A34,A43,A44,
     X       BNI,B11,B33,B44,SHI,U1I,A33I,A34I,A43I,A44I,B33I,B44I,
     X       EPSA,EPSB,EPS1,ANORM,BNORM,B3344,B3344I
      INTEGER max
      LOGICAL MATZ
      double COMPLEX Z3
C
      double precision ZERO,ZONE,ZTWO
C
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
C     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
C     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
C     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
C     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
C     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
C     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
C          WITH REAL SUBDIAGONAL ELEMENTS,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
C
C        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
C          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
C          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
C
C        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
C          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
C          THE RATIOS ((ALFR+I*ALFI)/BETA),
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF AR(J,J-1) HAS NOT BECOME
C                     ZERO AFTER 50 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      ZTWO = 2.0D0
      ZONE = 1.0D0
      ZERO = 0.0D0
 
      IERR = 0
C     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0.0D0
      BNORM = 0.0D0
C
      DO 30 I = 1, N
         ANI = 0.0D0
         IF (I .NE. 1) ANI = ABS(AR(I,I-1))
         BNI = 0.0D0
C
         DO 20 J = I, N
            ANI = ANI + ABS(AR(I,J)) + ABS(AI(I,J))
            BNI = BNI + ABS(BR(I,J)) + ABS(BI(I,J))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
C
      IF (ANORM .EQ. ZERO) ANORM = 1.0D0
      IF (BNORM .EQ. ZERO) BNORM = 1.0D0
      EP = EPS1
      IF (EP .GT. ZERO) GO TO 50
C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = 1.0D0
   40 EP = EP / ZTWO
      IF (ZONE + EP .GT. ZONE) GO TO 40
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     ********** REDUCE A TO TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = N
C     ********** BEGIN QZ STEP **********
   60 IF (EN .EQ. 0) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
   70 DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (ABS(AR(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 AR(L,LM1) = 0.0D0
C     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
   95 B11 = abs(CMPLX(BR(L,L),BI(L,L)))
      IF (B11     .EQ. ZERO) GO TO 98
      U1 = BR(L,L) / B11
      U1I = BI(L,L) / B11
C
      DO 97 J = L, ENORN
         XI = U1 * AI(L,J) - U1I * AR(L,J)
         AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
         AI(L,J) = XI
         XI = U1 * BI(L,J) - U1I * BR(L,J)
         BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
         BI(L,J) = XI
   97 CONTINUE
C
      BI(L,L) = 0.0D0
   98 IF (L .NE. EN) GO TO 100
C     ********** 1-BY-1 BLOCK ISOLATED **********
      ALFR(EN) = AR(EN,EN)
      ALFI(EN) = AI(EN,EN)
      BETA(EN) = B11
      EN = NA
      GO TO 60
C     ********** CHECK FOR SMALL TOP OF B **********
  100 L1 = L + 1
      IF (B11 .GT. EPSB) GO TO 120
      BR(L,L) = 0.0D0
      S = ABS(AR(L,L)) + ABS(AI(L,L)) + ABS(AR(L1,L))
      U1 = AR(L,L) / S
      U1I = AI(L,L) / S
      U2 = AR(L1,L) / S
      R = SQRT(U1*U1+U1I*U1I+U2*U2)
      U1 = U1 / R
      U1I = U1I / R
      U2 = U2 / R
      AR(L,L) = R * S
      AI(L,L) = 0.0D0
C
      DO 110 J = L1, ENORN
         XR = AR(L,J)
         XI = AI(L,J)
         YR = AR(L1,J)
         YI = AI(L1,J)
         AR(L,J) = U1 * XR + U1I * XI + U2 * YR
         AI(L,J) = U1 * XI - U1I * XR + U2 * YI
         AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
         XR = BR(L,J)
         XI = BI(L,J)
         YR = BR(L1,J)
         YI = BI(L1,J)
         BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         BR(L,J) = U1 * XR + U1I * XI + U2 * YR
         BI(L,J) = U1 * XI - U1I * XR + U2 * YI
         BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110 CONTINUE
C
      LM1 = L
      L = L1
      GO TO 90
C     ********** ITERATION STRATEGY **********
  120 IF (ITS .EQ. 50) GO TO 1000
      IF (ITS .EQ. 10) GO TO 135
C     ********** DETERMINE SHIFT **********
      B33 = BR(NA,NA)
      B33I = BI(NA,NA)
      IF (abs(CMPLX(B33,B33I)) .GE. EPSB) GO TO 122
      B33 = EPSB
      B33I = 0.0D0
  122 B44 = BR(EN,EN)
      B44I = BI(EN,EN)
      IF (abs(CMPLX(B44,B44I)) .GE. EPSB) GO TO 124
      B44 = EPSB
      B44I = 0.0D0
  124 B3344 = B33 * B44 - B33I * B44I
      B3344I = B33 * B44I + B33I * B44
      A33 = AR(NA,NA) * B44 - AI(NA,NA) * B44I
      A33I = AR(NA,NA) * B44I + AI(NA,NA) * B44
      A34 = AR(NA,EN) * B33 - AI(NA,EN) * B33I
     X    - AR(NA,NA) * BR(NA,EN) + AI(NA,NA) * BI(NA,EN)
      A34I = AR(NA,EN) * B33I + AI(NA,EN) * B33
     X     - AR(NA,NA) * BI(NA,EN) - AI(NA,NA) * BR(NA,EN)
      A43 = AR(EN,NA) * B44
      A43I = AR(EN,NA) * B44I
      A44 = AR(EN,EN) * B33 - AI(EN,EN) * B33I - AR(EN,NA) * BR(NA,EN)
      A44I = AR(EN,EN) * B33I + AI(EN,EN) * B33 - AR(EN,NA) * BI(NA,EN)
      SH = A44
      SHI = A44I
      XR = A34 * A43 - A34I * A43I
      XI = A34 * A43I + A34I * A43
      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 140
      YR = (A33 - SH) / 2.0D0
      YI = (A33I - SHI) / 2.0D0
      Z3 = sqrt(CMPLX(YR**2-YI**2+XR,2.0D0*YR*YI+XI))
      U1 = REAL(Z3)
      U1I = AIMAG(Z3)
      IF (YR * U1 + YI * U1I .GE. ZERO) GO TO 125
      U1 = -U1
      U1I = -U1I
  125 Z3 = (CMPLX(SH,SHI) - CMPLX(XR,XI) / CMPLX(YR+U1,YI+U1I))
     X   / CMPLX(B3344,B3344I)
      SH = REAL(Z3)
      SHI = AIMAG(Z3)
      GO TO 140
C     ********** AD HOC SHIFT **********
  135 SH = AR(EN,NA) + AR(NA,ENM2)
      SHI = 0.0D0
C     ********** DETERMINE ZEROTH COLUMN OF A **********
  140 A1 = AR(L,L) / B11 - SH
      A1I = AI(L,L) / B11 - SHI
      A2 = AR(L1,L) / B11
      ITS = ITS + 1
      IF (.NOT. MATZ) LOR1 = L
C     ********** MAIN LOOP **********
      DO 260 K = L, NA
         K1 = K + 1
         K2 = K + 2
         KM1 = max(K-1,L)
C     ********** ZERO A(K+1,K-1) **********
         IF (K .EQ. L) GO TO 170
         A1 = AR(K,KM1)
         A1I = AI(K,KM1)
         A2 = AR(K1,KM1)
  170    S = ABS(A1) + ABS(A1I) + ABS(A2)
         U1 = A1 / S
         U1I = A1I / S
         U2 = A2 / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
C
         DO 180 J = KM1, ENORN
            XR = AR(K,J)
            XI = AI(K,J)
            YR = AR(K1,J)
            YI = AI(K1,J)
            AR(K,J) = U1 * XR + U1I * XI + U2 * YR
            AI(K,J) = U1 * XI - U1I * XR + U2 * YI
            AR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            AI(K1,J) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(K,J)
            XI = BI(K,J)
            YR = BR(K1,J)
            YI = BI(K1,J)
            BR(K,J) = U1 * XR + U1I * XI + U2 * YR
            BI(K,J) = U1 * XI - U1I * XR + U2 * YI
            BR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            BI(K1,J) = U1 * YI + U1I * YR - U2 * XI
  180    CONTINUE
C
         IF (K .EQ. L) GO TO 240
         AI(K,KM1) = 0.0D0
         AR(K1,KM1) = 0.0D0
         AI(K1,KM1) = 0.0D0
C     ********** ZERO B(K+1,K) **********
  240    S = ABS(BR(K1,K1)) + ABS(BI(K1,K1)) + ABS(BR(K1,K))
         U1 = BR(K1,K1) / S
         U1I = BI(K1,K1) / S
         U2 = BR(K1,K) / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
         IF (K .EQ. NA) GO TO 245
         XR = AR(K2,K1)
         AR(K2,K1) = U1 * XR
         AI(K2,K1) = -U1I * XR
         AR(K2,K) = -U2 * XR
C
  245    DO 250 I = LOR1, K1
            XR = AR(I,K1)
            XI = AI(I,K1)
            YR = AR(I,K)
            YI = AI(I,K)
            AR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            AI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            AR(I,K) = U1 * YR - U1I * YI - U2 * XR
            AI(I,K) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(I,K1)
            XI = BI(I,K1)
            YR = BR(I,K)
            YI = BI(I,K)
            BR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            BI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            BR(I,K) = U1 * YR - U1I * YI - U2 * XR
            BI(I,K) = U1 * YI + U1I * YR - U2 * XI
  250    CONTINUE
C
         BI(K1,K1) = 0.0D0
         BR(K1,K) = 0.0D0
         BI(K1,K) = 0.0D0
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            XR = ZR(I,K1)
            XI = ZI(I,K1)
            YR = ZR(I,K)
            YI = ZI(I,K)
            ZR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            ZI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            ZR(I,K) = U1 * YR - U1I * YI - U2 * XR
            ZI(I,K) = U1 * YI + U1I * YR - U2 * XI
  255    CONTINUE
C
  260 CONTINUE
C     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
      IF (AI(EN,NA) .EQ. ZERO) GO TO 70
      R = abs(CMPLX(AR(EN,NA),AI(EN,NA)))
      U1 = AR(EN,NA) / R
      U1I = AI(EN,NA) / R
      AR(EN,NA) = R
      AI(EN,NA) = 0.0D0
C
      DO 270 J = EN, ENORN
         XI = U1 * AI(EN,J) - U1I * AR(EN,J)
         AR(EN,J) = U1 * AR(EN,J) + U1I * AI(EN,J)
         AI(EN,J) = XI
         XI = U1 * BI(EN,J) - U1I * BR(EN,J)
         BR(EN,J) = U1 * BR(EN,J) + U1I * BI(EN,J)
         BI(EN,J) = XI
  270 CONTINUE
C
      GO TO 70
C     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
C                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
 1000 IERR = EN
C     ********** SAVE EPSB FOR USE BY CQZVEC **********
 1001 IF (N .GT. 1) BR(N,1) = EPSB
      RETURN
C     ********** LAST CARD OF CQZVAL **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE R8CQZVEC(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI)
C
      IMPLICIT NONE
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN
      double precision AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),
     X       ALFR(N),ALFI(N),
     X       BETA(N),ZR(NM,N),ZI(NM,N)
      double precision R,T,RI,TI,XI,ALMI,ALMR,BETM,EPSB
      double COMPLEX Z3
      double precision ZERO
C
C
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
C     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
C     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
C     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
C          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
C          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
C
C        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
C
C        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        A IS UNALTERED,
C
C        B HAS BEEN DESTROYED,
C
C        ALFR, ALFI, AND BETA ARE UNALTERED,
C
C        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
C          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      ZERO = 0.0D0
      IF (N .LE. 1) GO TO 1001
      EPSB = BR(N,1)
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         NA = EN - 1
         ALMR = ALFR(EN)
         ALMI = ALFI(EN)
         BETM = BETA(EN)
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            R = 0.0D0
            RI = 0.0D0
            M = I + 1
C
            DO 610 J = M, EN
               T = BETM * AR(I,J) - ALMR * BR(I,J) + ALMI * BI(I,J)
               TI = BETM * AI(I,J) - ALMR * BI(I,J) - ALMI * BR(I,J)
               IF (J .EQ. EN) GO TO 605
               XI = T * BI(J,EN) + TI * BR(J,EN)
               T = T * BR(J,EN) - TI * BI(J,EN)
               TI = XI
  605          R = R + T
               RI = RI + TI
  610       CONTINUE
C
            T = ALMR * BETA(I) - BETM * ALFR(I)
            TI = ALMI * BETA(I) - BETM * ALFI(I)
            IF (T .EQ. ZERO .AND. TI .EQ. ZERO) T = EPSB
            Z3 = CMPLX(R,RI) / CMPLX(T,TI)
            BR(I,EN) = REAL(Z3)
            BI(I,EN) = AIMAG(Z3)
  700    CONTINUE
C
  800 CONTINUE
C     ********** END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 2 DO -- **********
      DO 880 JJ = 2, N
         J = N + 2 - JJ
         M = J - 1
C
         DO 880 I = 1, N
C
            DO 860 K = 1, M
               ZR(I,J) = ZR(I,J) + ZR(I,K) * BR(K,J) - ZI(I,K) * BI(K,J)
               ZI(I,J) = ZI(I,J) + ZR(I,K) * BI(K,J) + ZI(I,K) * BR(K,J)
  860       CONTINUE
C
  880 CONTINUE
C     ********** NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1 **********
      DO 950 J = 1, N
         T = 0.0D0
C
         DO 930 I = 1, N
            R = abs(CMPLX(ZR(I,J),ZI(I,J)))
            IF (R .GT. T) T = R
  930    CONTINUE
C
         DO 940 I = 1, N
            ZR(I,J) = ZR(I,J) / T
            ZI(I,J) = ZI(I,J) / T
  940    CONTINUE
C
  950 CONTINUE
C
 1001 RETURN
C     ********** LAST CARD OF CQZVEC **********
      END
