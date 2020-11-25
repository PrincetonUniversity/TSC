subroutine w20diff( zz, zvr, zvi, G_ave, kps, diffs, conv )
  use w20data
  implicit none

  integer, parameter :: ichn = 6

  complex*16, intent(in), dimension(neq) :: zz
  real*8, intent(in), dimension(neq,neq) :: zvr, zvi
  real*8, intent(inout), dimension(6) :: diffs, conv
  real*8, intent(in) :: G_ave, kps

  !Internal variables
  integer :: i,j                        !Looping indeces

  !Printing flag
  integer :: iw = 0

  !Pointer to axis point (temporary, must undo all this)
  integer :: ik = 1

  !Unstable mode
  complex*16 :: wu                    !Unstable mode complex frequency
  real*8 :: wu_sq                     !Magnitude of wu squared
  real*8 :: re_wu, im_wu              !Real and imaginary part of unstable mode
  real*8 :: im_zz                     !Imaginary part of unstable mode without ExB shear stabilization
  real*8 :: im_wu_inv                 !Inverse of the growth rate
  real*8 :: wu_res                    !Residual part of omega unstable?
  integer :: iunst                    !Number of unstable modes

  !Some physics quantities
  real*8 :: d1

  !Scaling factor for geometric correction (1/gdro2 in Jetto code)
  real*8 :: gk = 1.0

  !Fractions
  real*8 :: thrd2
  real*8 :: thrd = 0.3333333333

  !Random internal variables of unknown origin, destination, or meaning
  real*8 :: hx

  !Local quantities copied from module and bounded
  real*8 :: sn

  !eta_i threshold
  real*8 :: eith

  !Various quantities related to diffusivities
  real*8, dimension(ichn) :: xi, xe, xd, xq, xdq, xtv, hpt, dm, dmt, &
       chi, che, d, chq, dhq
  real*8 :: hp, shp, shpe, dmi, hpe

  !Local diffusivities and accumulators for diffusivities
  real*8 :: chic, schef, def, dqeff, chqeff, schq, ceft, deft
  real*8 :: schi, gci, sche, gce, sd, gd, sdm, gcq, sdq, gnq, gei
  real*8 :: xih, xqh, xeh, xdh, xhh

  real*8 :: pe, pq, pn, pnq,tiq, pi, piq

  !Possibly, non-diagonal elements and pinches from EM effects
  real*8 :: dms1, dms2, dms21, dms22, dmi21, dmi22, smp, smpt, bv, dmip, adisp, bdisp, &
       dt, dmit, dms, dmst, dmd, dmdn, dmdt, dmdn2, dmdt2, dmi1, dmi2, dadk, dbdk

  !Variables for output of poloidal, toroidal, perpendicular, thermoelectric
  !components of the pinch velocities
  real*8::v_p_tor,v_p_pol,v_p_per,v_p_thr

  !Collision related parameters
  complex*16 :: ga, gb, gm2
  real*8 :: re_ga, im_ga, re_gb, im_gb
  real*8 :: hr, xde, xdi, yda

  !Related to momentum transport
  real*8 :: dmh, dmef, dmef1, dmef2, dmef21, dmef22, dmdiag, dmdiagt, dmeftest
  real*8 :: dmeff, dmtef, dmteff
  real*8 :: sp, smef1, smef21, smef22, sre
  real*8 :: tsour
  real*8 :: kxi
  real*8 :: ainf, imf, tc, dp1, gp1, gp2
  real*8 :: elfs

  !Lengthscales related to momentum transport not yet purified
  real*8 :: lvf, lvft

  !Variables related to the parallel wavenumber in momentum transoprt computation
  real*8 :: Kkap, Hm, kap1, acn, RDexp, ks

  !Magnitude square, real and imaginary parts of various perturbations
  real*8 :: nehat_sq, re_nehat, im_nehat, nzhat_sq, re_nzhat, im_nzhat
  real*8 :: re_eln, im_eln
  real*8 :: dni, dne, dnz, di, de, dz, dnin
  real*8 :: dznt, dzn1, dzn2
  real*8 :: dtimp
  real*8 :: dqn, dqt, dtiq
  real*8 :: phs

  !More unknown variables
  !Random pieces of equations which need variables
  real*8 :: tt, et, ht, h, k, t, ts, t1, t2, a3
  real*8 :: c1, c2, e1, e2, h1, h2, at, ct, a2, a1, a, b, c, e, dh, f
  real*8 :: stf, eqh

  !Divergences of stuff???
  real*8 :: divgat, divat, devgat, devat, diva3, devga1, devga2, devga3, deva1, deva2, deva3
  real*8 :: divga1, divga2, divga3, diva1, diva2, deva, devb
  real*8 :: divga, divgb, diva, divb, devga, devgb

  !Sources of stuff????
  real*8 :: svat, sva1, sva2, sva3, sva, svb

  !Effective collisionality normalized by electron perturbation magnitude square
  real*8 :: vefn

  !Complex variables from Weiland's code
  COMPLEX*16 :: Fm,hf,EMP,EMPT,VPL,Vg,WRES
  COMPLEX*16 :: FIH,NEF,NRAT,ELN,TEF,FRAT,AV,NIF,MRT,HQ,CDexp
  COMPLEX*16 :: RNI,TII,RTI,TICF,CN,AVRAT,MRT1,MRT2,NII
  COMPLEX*16 :: nehat,nzhat,BT2
  COMPLEX*16 :: IU,HC,WZ,WZP
  COMPLEX*16 :: DMSP,KPF,kpf1,kpf2,KPX,ELMS,KPDX,TCONT,NCONT,EMCONT

  wde = 2.0 * csound * kyrho / rmaj0

  !Lengthscales of velocity as defined by J. Weiland
  !Note: these are replaced by gradients (gvt, gvp) in most places

  lvf   = rmaj / ( amin * gvp )
  lvft  = rmaj / ( amin * gvt )

  !Initialize some variables
  diffs = 0.0
!  vconv = 0.0
  iunst = 0
  wu    = ( 0.0, 0.0 ) 

  THRD2=THRD*THRD
  IU=(0.D0,1.D0)
  SHP=0.D0
  CHIC=0.
  SHPE=0.
  SCHEF=0.
  DEF=0.
  EQH=(gtz/gnz)-STR+FTR*(2.0/gnz)
  d1 = 6.616*dsqrt(ahyd)/(rmaj*btor**2)

  !Two missing variables
  RDexp = 1.0

  !Initialize diffusion matrix elements 
  chi(1:6)=0.0 ; che(1:6)=0.0 ; d(1:6)  =0.0 ; chq(1:6)=0.0 
  dm(1:6) =0.0 ; dhq(1:6)=0.0 ; dmt(1:6)=0.0 ; hpt(1:6)=0.0

  !Initialize effective diffusivity
  SCHI=0.D0 ; SCHE=0.D0 ; SD  =0.D0 ; SDM =0.D0 ; SDQ =0.D0 ; XHH =0.D0

  !Main loop over unstable modes
  !(further down is second loop for momentum transport)
  !Thermal and particle diffusivities are computed here

  main_modes: do j=1,neq

     re_wu=DREAL(zz(j))
     im_wu=DIMAG(zz(j))

     !Cycle when growth rate is negative (stable mode)
     if (im_wu.le.0.0) cycle main_modes

     !Notice different substraction of ExB shear for
     !ion modes and electron modes
     if (re_wu.lt.0.0) then
        im_wu = im_wu - dabs(wexb)
     else if ( dreal(zz(j)).gt.0.0) then
        wu_res = im_wu**2-0.25*wexb**2
        if (wu_res .lt. 0.0) cycle
        im_wu = sqrt(wu_res)
     end if

     wu=re_wu+(0.0,1.0)*im_wu

     if(im_wu .le. 1.0E-3) then
        cycle main_modes
     end if

     !Define complex frequency of mode after ExB shear effect
     wu=CMPLX(re_wu,im_wu)
     WRES=wu+FTRT
     gm2=1.0+0.5*gte/(wu-1.0+IU*vef)
     BT2=BT-2.5*GM2
     HC=wu-FTR+TVR*BT1
     nehat=wu*wu-2.D0*FTR*wu+FTR+IU*VEF*HC
     wu_sq=re_wu*re_wu+im_wu*im_wu
     im_wu_inv=1.0/im_wu
     !
     !   ******   COLLISION PARAMETERS  *******
     !
     GA=wu-FTR+TVR*BT2
     GB=(wu-FTR)/(wu-1.D0+IU*VEF)
     re_ga=DREAL(GA)
     im_ga=DIMAG(GA)
     re_gb=DREAL(GB)
     im_gb=DIMAG(GB)

     HR=re_wu-FTR+TVR*BT1
     XDE=wu_sq-FTR*re_wu
     XDI=wu_sq+FTR*tau_inv*re_wu
     YDA=re_wu*(1.D0-(2.0/gne))+(gte/gne)-STR
     !   ***************************************
     !

     re_nehat=DREAL(nehat)
     im_nehat=DIMAG(nehat)
     nehat_sq=(re_nehat)**2+(im_nehat)**2
     !
     !      IF(VEF.EQ.0.) GO TO 25
     DIVGA=FTRT*((2.0/gne)*(re_ga*im_nehat-im_ga*re_nehat) &
          -YDA*im_wu+im_wu*HR*(1.-(2.0/gne)))

     DIVGB=FTRT*(re_gb*im_nehat-im_gb*re_nehat-im_wu)

     DIVA=XDI*((2.0/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1.-(2.0/gne))-YDA*HR)

     DIVB=XDI*(re_gb*re_nehat+im_gb*im_nehat-HR)

     DEVGA=BT1*((YDA-VEF*im_ga*(2.0/gne))*im_nehat    &
          -(im_wu*(1.-(2.0/gne))+VEF*re_ga*(2.0/gne)) &
          *re_nehat) &
          +FTR*(im_wu*YDA+(2.0/gne)*(im_ga*re_nehat-re_ga*im_nehat)-im_wu*HR*(1.-(2.0/gne)))

     DEVGB=BT1*((1.-VEF*im_gb)*im_nehat-VEF*re_gb*re_nehat) &
          +FTR*(im_wu+im_gb*re_nehat-re_gb*im_nehat)

     DEVA=XDE*((2.0/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1.-(2.0/gne))-YDA*HR) &
          -BT1*(re_wu-FTR)*((YDA-VEF*im_ga*(2.0/gne))*re_nehat+(im_wu*(1.-(2.0/gne)) &
          +VEF*re_ga*(2.0/gne))*im_nehat)

     DEVB=XDE*(re_gb*re_nehat+im_gb*im_nehat-HR)-BT1*(re_wu-FTR)*((1.- &
          VEF*im_gb)*re_nehat+VEF*re_gb*im_nehat)

     SVA=(2.0/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1.-(2.0/gne))-YDA*HR

     SVB=re_gb*re_nehat+im_gb*im_nehat-HR

     VEFN=VEF/nehat_sq

     DIVGA=DIVGA*VEFN
     DIVGB=DIVGB*VEFN
     DEVGA=DEVGA*VEFN
     DEVGB=DEVGB*VEFN
     DIVA=DIVA*VEFN
     DIVB=DIVB*VEFN
     DEVA=DEVA*VEFN
     DEVB=DEVB*VEFN
     SVA=SVA*VEFN
     SVB=SVB*VEFN
     !   25 CONTINUE
     !
     STF=STR-FTR*(2.0/gne)
     A1=wu_sq*((2.0/gne)-1.D0)+2.D0*re_wu*STF
     A=wu_sq*(A1+FTR*(STR*(2.0/gne)-11.D0*THRD-tau_inv*(1.D0- &
          FTR*(2.0/gne))))+FTR*FTR*tau_inv*(2.D0*re_wu*(1.D0-(2.0/gne))-STF)
     A=A/nehat_sq
     B=(wu_sq*(2.D0*(re_wu-FTR)+FTR*tau_inv)- &
          FTR*FTR*tau_inv)/nehat_sq
     C=wu_sq*(A1+TVR*FTR*((2.0/gne)-4.D0))+ &
          FTR*FTR*(2.D0*re_wu*((2.0/gne)-1.D0)+STF)
     C=C/nehat_sq
     DH=(wu_sq*(2.D0*re_wu-5.D0)+FTR*FTR)/nehat_sq
     E=wu_sq*(1.D0-(2.0/gne))-2.D0*re_wu*STF+FTR*(11.D0*THRD &
          -STR*(2.0/gne))
     E=E/nehat_sq
     F=2.D0*(-re_wu+FTR)/nehat_sq
     !
     re_nzhat=re_wu**2-im_wu**2+2.D0*FTR*ztauz*re_wu+FTR*ztauz*ztauz
     im_nzhat=2.D0*im_wu*(re_wu+FTR*ztauz)
     nzhat_sq=re_nzhat**2+im_nzhat**2
     !
     !   ****   SPLITTING IN  EN  AND EE  ************
     !
     !
     !  Ion thermal conductivity
     !
     DIVGA1=FTRT*((STR-re_wu)*im_wu+im_wu*HR)*VEFN
     DIVGA2=FTRT*(re_ga*im_nehat-im_ga*re_nehat+re_wu*im_wu-im_wu*HR)*VEFN
     DIVGA3=-FTRT*im_wu*VEFN

     DIVA1=XDI*(-im_wu*im_wu+(STR-re_wu)*HR)*VEFN
     DIVA2=XDI*(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)*VEFN
     DIVA3=-XDI*HR*VEFN

     !
     !   Electron thermal conductivity
     !
     DEVGA1=(BT1*((re_wu-STR)*im_nehat-im_wu*re_nehat)+FTR*(im_wu*(re_wu-STR)-im_wu*HR))&
          *VEFN
     DEVGA2=(BT1*((-re_wu-VEF*im_ga)*im_nehat+(im_wu-VEF*re_ga)*re_nehat)&
          +FTR*(-im_wu*re_wu+im_ga*re_nehat-re_ga*im_nehat+im_wu*HR))*VEFN
     DEVGA3=(BT1*im_nehat+FTR*im_wu)*VEFN

     DEVA1=(XDE*(-im_wu*im_wu-(re_wu-STR)*HR)-BT1*(re_wu-FTR)*((re_wu-STR)*re_nehat&
          +im_wu*im_nehat))*VEFN
     DEVA2=(XDE*(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)-BT1*(re_wu-FTR)&
          *((-re_wu-VEF*im_ga)*re_nehat-(im_wu-VEF*re_ga)*im_nehat))*VEFN
     DEVA3=(-XDE*HR-BT1*(re_wu-FTR)*re_nehat)*VEFN

     SVA1=(-im_wu*im_wu-(re_wu-STR)*HR)*VEFN
     SVA2=(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)*VEFN
     SVA3=-HR*VEFN

     !------------------------------------------------------------------------ 
     ! Contributions to the ion conductivity 
     !------------------------------------------------------------------------ 

     A2=(wu_sq*wu_sq+FTR*((STR+FTR*tau_inv-2.D0*re_wu)*wu_sq+FTR*tau_inv*(FTR-&
          2.D0*re_wu)))/nehat_sq
     A3=(wu_sq*(-wu_sq+2.D0*STR*re_wu-FTR*(11.D0*THRD+tau_inv))+FTR*FTR&
          *tau_inv*(2.D0*re_wu-STR))/nehat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the electron conductivity 
     !------------------------------------------------------------------------ 

     C1=(nehat_sq-wu_sq*wu_sq+14.D0*THRD*wu_sq*re_wu-40.D0*THRD2*wu_sq&
          -50.D0*THRD2*re_wu+175.D0/27.D0)/nehat_sq
     C2=(wu_sq*wu_sq-10.D0*THRD*wu_sq*re_wu+10.D0*THRD2*wu_sq&
          +50.D0*THRD2*re_wu-125.D0/27.D0)/nehat_sq

     !------------------------------------------------------------------------ 
     !   Contributions to the electron diffusivity 
     !------------------------------------------------------------------------ 

     E1=(wu_sq+11.D0*THRD*FTR-2.D0*re_wu*STR)/nehat_sq
     E2=-(-2.D0*re_wu*FTR+(wu_sq+STR*FTR))/nehat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the main ion Conductivity 
     !------------------------------------------------------------------------ 

     H1=(wu_sq*(-wu_sq-2.D0*ztauz*re_wu*STR&
          +FTR*ztauz*ztauz*(-11.D0*THRD)+FTR*tau_inv*ztauz)&
          +FTR*FTR*tau_inv*ztauz*ztauz*(2.D0*re_wu+STR*ztauz))/nzhat_sq

     H2=(wu_sq*(wu_sq+2.D0*ZTAUZ*re_wu*FTR+FTR*ZTAUZ*ZTAUZ*STR-FTR*tau_inv*ZTAUZ*FTR)&
          -FTR*FTR*tau_inv*ZTAUZ*ZTAUZ*(2.D0*re_wu+FTR*ZTAUZ))/nzhat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the impurity conductivity 
     !------------------------------------------------------------------------ 

     T1=(wu_sq*(-wu_sq-2.D0*ZTAUZ*re_wu*STR-FTR*ZTAUZ*ZTAUZ*8.D0*THRD)&
          +FTR*FTR*ZTAUZ**3*(2.D0*re_wu+ZTAUZ*STR))/nzhat_sq
     T2=(wu_sq*(wu_sq+2.D0*ZTAUZ*re_wu*FTR+FTR*ZTAUZ*ZTAUZ*TVR)&
          -FTR*FTR*ZTAUZ**3*(2.d0*re_wu+ZTAUZ*FTR))/nzhat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the impurity diffusivity 
     !------------------------------------------------------------------------ 

     DZN1=(-re_nzhat+2.D0*(re_wu+ZTAUZ*STR)*(re_wu+FTR*ZTAUZ))/nzhat_sq
     DZN2=(re_nzhat-2.D0*(re_wu+ZTAUZ*FTR)*(re_wu+FTR*ZTAUZ))/nzhat_sq
     !
     DIVGAT=DIVGA1+(2.0/gne)*DIVGA2+(gte/gne)*DIVGA3
     DIVAT=DIVA1+(2.0/gne)*DIVA2+(gte/gne)*DIVA3
     DEVGAT=DEVGA1+(2.0/gne)*DEVGA2+(gte/gne)*DEVGA3
     DEVAT=DEVA1+(2.0/gne)*DEVA2+(gte/gne)*DEVA3
     SVAT=SVA1+(2.0/gne)*SVA2+(gte/gne)*SVA3
     AT=(2.0/gne)*A2+A3
     CT=C1-1.+(2.0/gne)*C2
     ET=E1+(2.0/gne)*E2
     HT=H1+(2.0/gnz)*H2
     TT=T1+(2.0/gnz)*T2
     DZNT=DZN1+(2.0/gne)*DZN2
     !
     ! **** IMPURITIES *****
     !
     nzhat_sq=(re_wu*(re_wu+2.*FTR*ZTAUZ)-im_wu*im_wu+FTR*ZTAUZ*ZTAUZ)**2 &
          +4.*im_wu*im_wu*(re_wu+FTR*ZTAUZ)**2
     DTIMP=(wu_sq*(wu_sq*((2.0/gnz)-1.)+2.*ZTAUZ*re_wu*EQH+FTR*ZTAUZ &
          *ZTAUZ*(2.*(gtz/gnz)-11./3.+STR*(2.0/gnz))+FTR*tau_inv*ZTAUZ*(1. &
          +(gtz/gnz)-FTR*(2.0/gnz)))+FTR*tau_inv*(2.*FTR*ZTAUZ*ZTAUZ*re_wu*(1.-(2.0/gnz)) &
          -FTR*ZTAUZ**3*EQH))/nzhat_sq
     !  *************
     !
     DNI=(re_wu+FTR*tau_inv)**2+im_wu*im_wu
     DNE=(re_wu-FTR)**2+im_wu*im_wu
     DNZ=(re_wu+FTR*ZTAUZ)**2+im_wu*im_wu
     DI=(gti/gni)*kyrho*DNI
     DE=(gte/gne)*kyrho*DNE
     DZ=(gtz/gnz)*kyrho*DNZ
     !
     ! **** IMPURITIES *****
     !
     nzhat_sq=re_nzhat**2+im_nzhat**2
     !
     H=(wu_sq*(wu_sq*((2.0/gnz)-1.)-2.*ZTAUZ*re_wu*(STR-FTR*(2.0/gnz)) &
          +FTR*ZTAUZ*ZTAUZ*(-11.*THRD+STR*(2.0/gnz))+FTR*tau_inv*ZTAUZ*(1.-FTR &
          *(2.0/gnz)))+FTR*FTR*tau_inv*ZTAUZ*ZTAUZ*(2.*re_wu*(1.-(2.0/gnz))+(STR-FTR &
          *(2.0/gnz))*ZTAUZ))/nzhat_sq
     K=ZTAUZ*(FTR*FTR*tau_inv*ZTAUZ*ZTAUZ-wu_sq*(2.*re_wu &
          +FTR*(2.*ZTAUZ+tau_inv)))/nzhat_sq
     !
     !  *************
     !
     T=(wu_sq*(wu_sq*((2.0/gnz)-1.)-2.*ZTAUZ*re_wu*(STR-FTR*(2.0/gnz)) &
          +FTR*ZTAUZ*ZTAUZ*(-8.*THRD+TVR*(2.0/gnz)))+FTR*FTR*(ZTAUZ)**3 &
          *(2.*re_wu*(1.-(2.0/gnz))+ZTAUZ*(STR-FTR*(2.0/gnz))))/nzhat_sq
     TS=ZTAUZ*(FTR*FTR*(ZTAUZ)**3-wu_sq*(2.*re_wu+5.*ZTAUZ))/nzhat_sq
     !
     T=(GNZ/GNE)*T
     TS=(GNZ/GNE)*TS
     DQN=((2.0/gnz)-1.)*re_nzhat+2.*((1.-(2.0/gnz))*re_wu+ZTAUZ*(STR-FTR*(2.0/gnz)) &
          )*(re_wu+FTR*ZTAUZ)
     DQT=2.*ZTAUZ*(re_wu+FTR*ZTAUZ)
     DQN=DQN/nzhat_sq
     DQT=DQT/nzhat_sq
     !
     DTIQ=H-(gtz/gnz)*K
     !      DTQQ=T-EQ*TS
     !
     PHS=(re_wu-FTR)*DREAL(BT2)+im_wu*DIMAG(BT2)
     XDH=D1*(TE**1.5D0)*im_wu**3/kyrho
     XIH=XDH/DNI
     XEH=XDH/DNE
     XQH=XDH/DNZ

     !----------------------
     ! Ion transport channel
     !----------------------
     ! Diffusive part
     !
     XI(1)=XIH
     XI(2)=tvr*fte*xih*ne_nh*(B-DIVGA3-DIVGB-im_wu_inv*(DIVA3+DIVB))
     XI(3)=-tvr*xih*(gni/gne+fte*ne_nh*(A3+DIVGA1+DIVA1*im_wu_inv))
     !      XI(4)=-XIH*TVR*nz_ne*Z*ne_nh*K*TAUZ*tau_inv
     XI(4)=0.
     XI(5)=tvr*xih*zimp*h1/nh
     XI(6)=0.d0
     !
     ! Convective part
     !
     PI=fte*TVR*XIH*ne_nh*(A2+DIVGA2+im_wu_inv*DIVA2)*amin/rmaj0
     PIQ=-TVR*XIH*zimp*H2*nz_ne*ne_nh*amin/rmaj0
     HP=XIH*ne_nh*tau_inv*TVR*FTR*(1.D0-fte)*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GCI = XI(1) &
          + XI(2)*gte/gti &
          + XI(3)*gne/gti &
          + XI(5)*nz*gne/gti &
          - (HP+PI+PIQ)*2.0*gk*rmaj/amin/gti

     !-----------------------------------
     ! Electron thermal transport channel
     !-----------------------------------
     ! Diffusive part
     !
     XE(1)=0.D0
     XE(2)=fte*XEH*(1.D0+TVR*(DH-DEVGB-DEVGA3-im_wu_inv*(DEVB+DEVA3)))
     XE(3)=-TVR*fte*XEH*(C1+DEVGA1+im_wu_inv*DEVA1)
     XE(4)=0.D0
     XE(5)=0.D0
     XE(6)=0.D0
     !
     ! Convective part
     !
     PE=fte*TVR*XEH*(C2+DEVGA2+im_wu_inv*(DEVA2+VEF*PHS))*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GCE = XE(2) &
          + XE(3)*gne/gte &
          - PE*gk*2.0*rmaj/amin/gte

     !------------------------------------
     ! Electron particle transport channel
     !------------------------------------
     ! Diffusive part
     !
     XD(1)=0.D0
     XD(2)=-fte*XDH*(ne/te)*(F+im_wu_inv*(SVB+SVA3))
     XD(3)=fte*XDH*(E1-im_wu_inv*SVA1)
     XD(4)=0.D0
     XD(5)=0.D0
     XD(6)=0.D0
     !
     ! Convective part
     !
     PN=-fte*XDH*(E2-im_wu_inv*SVA2)*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GD = XD(2)*te/ne*gte/gne &
          + XD(3) &
          -  PN*gk*2.0*rmaj/amin/gne

     !-----------------------------------
     ! Impurity thermal transport channel
     !-----------------------------------
     ! Diffusive part
     !
     XQ(1)=0.D0
     XQ(2)=0.D0
     XQ(3)=0.D0
     XQ(4)=1.D0
     XQ(5)=0.D0
     XQ(6)=0.D0
     !
     ! Convective part
     !
     PQ=TVR*XQH*T2*2/rmaj
     !
     ! Effective diffusivity
     !
     !      GCQ=XQ(4)+XQ(5)*NQ/(TZ*(gtz/gnz))-PQ*GK*2.*LTZ

     GCQ = -XDH*DQT*NZ/TZ !XDQ(4)
     !------------------------------------
     ! Impurity particle transport channel
     !------------------------------------
     ! Diffusive part
     !
     XDQ(1)=0.D0
     XDQ(2)=0.D0
     XDQ(3)=0.D0
     XDQ(4)=-XDH*DQT*NZ/TZ
     XDQ(5)=XDH*DZN1
     XDQ(6)=0.D0
     !
     ! Convective part
     !
     PNQ=-XDH*DZN2*2/rmaj
     !
     ! Effective diffusivity
     !
     GNQ = XDQ(4)*TZ*(gtz/gnz)/NZ &
          + XDQ(5) &
          - PNQ*GK*rmaj/gnz

     XTV(1)=0.D0
     XTV(2)=0.D0
     XTV(3)=0.D0
     XTV(4)=0.D0
     XTV(5)=0.D0
     XTV(6)=1.D0
     !
     TIQ=1.D0/(TAUZ*(GNI/GNZ))

     SHP=SHP+HP
     !
     HPT(1)=HPT(1)+HP+PI+PIQ
     HPT(2)=HPT(2)+PE
     HPT(3)=HPT(3)+PN
     !      HPT(4)=HPT(4)+PQ
     HPT(5)=HPT(5)+PNQ
     !
     SCHI=SCHI+GCI
     SCHE=SCHE+GCE
     SD=SD+GD
     SDM=SDM+GCQ
     SDQ=SDQ+GNQ
     !
     !
     !      CHQEFF=D1*TE**1.5*WI**3*((gtz/gnz)-TVR-TVR*DTQQ)/DQ
     chqeff=0.0
     schq  =0.0
     DQEFF=XDH*(DQN-(gtz/gnz)*DQT)

     XHH=XHH+XIH

     Do i=1,ichn
        CHI(I)=CHI(I)+XI(I)
        CHE(I)=CHE(I)+XE(I)
        D(I)=D(I)+XD(I)
        !      CHQ(I)=CHQ(I)+XQ(I)
        DHQ(I)=DHQ(I)+XDQ(I)
     end do

  End Do main_modes
  ! 
  !  IF electromagnetic effects or collisions on free electrons are included
  !  the transport coefficients are corrected for this.
  !------------------------------------------------------------
  DT=D1*TE**1.5
  SHPE=0.D0
  SCHEF=0.D0
  DEF=0.D0
  DMI=0.D0
  TSOUR=0.D0
  DMIT=0.D0
  DMS=0.D0
  DMST=0.D0
  DMSP=(0.D0,0.D0)
  DMD=0.D0
  DMDT=0.D0
  DMDT2=0.D0
  DMI1=0.D0
  DMI2=0.D0
  DMS1=0.D0
  DMS2=0.D0
  DMI21=0.D0
  DMI22=0.D0
  SMP=0.D0
  SMPT=0.D0
  BV=0.D0
  DMIP=0.D0    !!  Used for summing up curvature pinch terms
  v_p_pol=0.0
  v_p_tor=0.0
  v_p_per=0.0
  v_p_thr=0.0
  adisp=((2.0/gne)-1.D0-flh*(FTRT*(2.0/gne)-tau_inv*(1.D0+(gti/gni))))/(1.D0+flh)
  bdisp=tau_inv*(2.0/gne)/(1.D0+fte)
  dadk=-2.D0*kyrho*rhos*(FTRT*(2.0/gne)-tau_inv*(1.D0+(gti/gni))+((2.0/gne)-1.D0)/(1.D0+flh))
  dadk=dadk/(1.D0+flh)
  dbdk=-2.D0*(2.0/gne)*tau_inv*kyrho*rhos/(1.D0+flh)**2
  EITH=TVR+10.D0*(2.0/gne)*tau_inv/9.D0
  KXI=2.D0*q*kyrho/(eps*rmaj0)
  !------------------------------------------- loop

  main_momentum: DO J=1,NEQ

     re_wu=DREAL(zz(j))
     im_wu=DIMAG(zz(j))

     if ( im_wu.lt.0.0 ) cycle

     if ( re_wu.lt.0.0 ) then !Ion mode
        im_wu = im_wu - dabs(wexb)
     else if ( re_wu .gt.0.0 ) then
        wu_res =im_wu**2-0.25*wexb**2
        if (wu_res.lt.0.0) cycle
        im_wu = dsqrt(wu_res)
     end if

     wu=re_wu+(0.0,1.0)*im_wu

     if(im_wu .le. 1.0E-2) then
        cycle
     end if

     wu=CMPLX(re_wu,im_wu)

     im_zz = max( 0.01, dimag(zz(j)) )

     ! --- contr. to chii from em free electr. ----
     !
     DNI=(re_wu+FTR*tau_inv)**2+im_wu*im_wu
     DNIN=(re_wu+FTR*tau_inv)**2+im_zz*im_zz
     wu_sq=re_wu*re_wu+im_wu*im_wu
     XDH=DT*im_wu**3/kyrho
     XIH=XDH/DNI
     FIH=DCMPLX(ZVR(1,J),ZVI(1,J))
     ! F. Halpern: place a lower bound on the absolute value of FIH 
     IF ( abs(fih) .lt. 0.001 ) then
        fih = (0.001,0.0)
     end if

     !-- Here NEF for disp 9 is defined --
     AV=DCMPLX(ZVR(8,J),ZVI(8,J))
     NEF=FIH-(ZZ(J)-0.5*gne)*AV/KPC
     !--------------------------------------------

     IF(ABS(FIH).LT.0.0001) FIH=(0.0001,0.)
     NRAT=NEF/FIH
     ELN=NRAT-1.
     re_eln=DREAL(ELN)
     im_eln=DIMAG(ELN)
     AINF=TVR*(FTR*tau_inv*re_eln+im_eln*(im_wu*im_wu+re_wu*(re_wu+FTR*tau_inv))/im_zz)
     HPE=XIH*ne_nh*(1.-fte)*eps0*AINF

     SHPE=SHPE+HPE
     ! ****  Free electron heat flux *********************
     !
     !-----------------------------------------------------
     IF(NEQ.EQ.11) GO TO 10099
     IF(NEQ.EQ.9) GO TO 10098
     !--- Here TEF for disp10 is defined ---
     TEF=DCMPLX(ZVR(6,J),ZVI(6,J))
     GO TO 10099
10098 CONTINUE
     !-- Here TEF for disp9 is defined ---
     !      TEF=(gte/gne)*0.5*gne*AV/KPC
     TEF = 0.5*gte*AV/KPC
     !-------------------------------------------------------------
10099 CONTINUE
     FRAT=TEF/FIH
     IMF=-IMAG(FRAT)*DNIN/DNI
     CEFT=(1.-fte)*IMF*DT*(2.0/gte)*im_wu**3/(kyrho*im_zz)
     SCHEF=SCHEF+CEFT
     !**********************************************************
     ! ----Free electron particle flux -----------
     GEI=-IMAG(NRAT)/im_wu
     DEFT=(1.D0-fte)*GEI*(2.0/gne)*XDH
     DEF=DEF+DEFT
     ! ----Momentum transport ------
     TII=DCMPLX(ZVR(2,J),ZVI(2,J))
     NIF=DCMPLX(ZVR(3,J),ZVI(3,J))
     AVRAT=AV/FIH
     MRT=(TII+NIF)/FIH
     MRT1=0.5*(gti*gne/gni-gne*tvr)/(wu+ftrt*G_ave)
     TICF=TII/FIH
     NII=NIF/FIH
     MRT2=NIF*(1.+TVR*wu/(wu+FTRT*G_ave))/FIH
     DMS21=im_wu*im_wu*tau_inv*DREAL(MRT1)
     DMS22=im_wu*im_wu*tau_inv*DREAL(MRT2)
     DMS1=im_wu*im_wu
     DMS2=im_wu*im_wu*tau_inv*DREAL(MRT)
     DMS=DMS1+DMS2
     DMI1=DMI1+DMS1
     DMI2=DMI2+DMS2
     DMI21=DMI21+DMS21    !! part of temp. pert 
     DMI22=DMI22+DMS22    !! density pert including part from temp. pert.
     DMI=DMI+DMS
     DMD=DMD+im_wu**3/wu_sq
     !Comparison between old and new diagonal element
     DMDN=(re_wu+2.*tau_inv*G_ave)**2+im_wu*im_wu
     DMDT=DMDT+im_wu**3/DMDN  !! Diagonal element for toroidal momentum transport
     TC=(re_wu+2.*tau_inv*G_ave)/DMDN  !! Correl. time for TEP and Termoel. mom pinch (Hahm)
     DP1=2.*DT*im_wu*im_wu/kyrho

! By F.Halpern: 18-Feb-2009: Introduce gp1 and gp2 as stated in Jan Weiland's MOMCODE from 25-Jan-2009

! By A Kritz: 31-Mar-2009: Moved line defining ELMS before it is used in computing gp2
     ELMS = (wu+0.5*tau_inv*(gne+gti*gne/gni))*AVRAT/KPC ! From Weiland's momcode

     gp1 = - im_wu * dp1 * vtor * G_ave / dmdn
     gp2 = dreal( IU * ( MRT - tau*em2*ELMS )/( wu + 2.0*tau_inv*G_ave ) ) * dp1 * vtor * G_ave

!     GP1=0.5*TC*DP1*vtor  !!  TEP  momentum pinch
!     GP2=DREAL(TICF/(wu+2.*tau_inv*G_ave))*DP1*vtor  !!  Termoel. momentum pinch


     ELFS = (1.D0+0.5*tau_inv*(gne+gti*gne/gni)/wu)*AVRAT/KPC ! From Weiland's momcode

! By F. Halpern: Add factor of average curvature as per 25-Jan-2009 MOMCODE
     DMSP=im_wu*im_wu*(IU/(wu+2.*tau_inv*G_ave))*(1+tau_inv*MRT-em2*ELMS) !! due to v_par incl TC
!     DMSP=im_wu*im_wu*(IU/(wu+2.*tau_inv))*(1+tau_inv*MRT-em2*ELMS) !! due to v_par incl TC

     !----- Computation of parallel wavenumber k_\parallel

     !K in Weiland's 2006 EPS paper
     !F. Halpern: replaced toroidal velocity gradient with parallel gradient
     !      Kkap = rmaj0/rmaj*q*kyrho*vtor*gvt
     Kkap = rmaj0/rmaj*q*kyrho*vpar*gvpar

     !\kappa_1 in Weiland's 2006 EPS paper
     !F. Halpern: replaced toroidal velocity gradient with parallel gradient

     kap1 = vpol*gvp/(gne*sign(max(abs(shear),0.1),shear)*kyrho)
     !      kap1 = vpol*gvpar/(gne*sign(max(abs(shear),0.1),shear)*kyrho*dsqrt(kappa))
     kap1 = sign(min(abs(kap1),1.0),kap1)

     !ks in Weiland's Jan-2009 MOMCODE
     ks   = -2.0*vtor*kyrho*q/(tau_inv*tau_inv)

     Fm=wu*(1.+FTRT)+tau_inv*0.5*(gti*gne/gni-TVR*gne)
     Hm=FTRT*(1.+tau_inv)
     hf=4.0*flh*q**2*wu*0.5*gne

     !Flux average of k_\parallel normalized, plus separate expressions 
     !poloidal and toroidal driving terms
     KPF=-(0.5*(wu+ftr)*(Kkap+ks) + IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0) !!This is dimk_par
     kpf1=-(0.5*(wu+ftr)*(Kkap+ks))/((Fm+Hm)*Q*Rmaj0)
     kpf2=-(IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0)

!     KPF=-(0.5*(wu+ftr)*Kkap+IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0) !!This is dimk_par
!     kpf1=-(0.5*(wu+ftr)*Kkap)/((Fm+Hm)*Q*Rmaj0)

     !Normalization is removed here
     KPF=KPF*CSound/Wde
     kpf1=kpf1*CSound/Wde
     kpf2=kpf2*CSound/Wde

     !Compute the fluxes here
     !Split the different pinches due to velocity shear
     ! v_p_per ---> perpendicular part
     ! v_p_pol ---> poloidal part of parallel velocity pinch
     ! v_p_tor ---> toroidal part of parallel velocity pinch

     v_p_per = v_p_per -2.D0*d1*te**1.5*(eps/q)*dms
     v_p_pol = v_p_pol +2.D0*d1*te**1.5*dreal(DMSP*kpf2)/kyrho
     v_p_tor = v_p_tor +2.D0*d1*te**1.5*dreal(DMSP*kpf1)/kyrho

     !Toroidal momentum pinch term
     DMST=-(eps/q)*dms +dreal(DMSP*kpf)/kyrho
     DMIT=DMIT+DMST

     !---- Summation of curvature pinch term fluxes here
     DMIP=DMIP+GP1+GP2     !!  Summation of curvature pinch terms

     v_p_thr = dmip

     !---- Quantities related to poloidal momentum transport

     !---- Obtained from Jan Weiland's momcode dated from Feb 4 2008
     !---- Not debugged or tested

     Vg=-((wu+ftrt)*dadk+0.5*(gti*gne/gni-EITH*gne)*dbdk) & 
          /(2.0*(wu+ftrt)+0.5*gne*adisp)
     RDexp=DREAL( 0.25*(gne*gti-eith*gne**2)*Vg/ (wu+ftrt)**2 )

     TSOUR=TSOUR+4.*RDexp*im_wu*im_wu*diffBohm*KXI*KPS/rmaj0**2
     BV=BV+im_wu*im_wu*DREAL(IU*NII*CONJG(NII-TICF))

     CN=wu+FTRT
     ACN=ABS(CN)

     ! End of looping through unstable modes
     ! -------
  End Do main_momentum

  !---- Final computation of momentum fluxes
  !---- First, some unknown computations

  CHIC=-SHPE*GK*2.*(rmaj/gti)
  HPT(1)=HPT(1)+EM*SHPE
  CHE(2)=CHE(2)+EM*SCHEF
  D(3)=D(3)+EM1*DEF
  !
  DMH=-2.D0*D1*TE**1.5*eps0*(Csound/CSOUND0)/vpol
  DMEF=DMH*DMI*LVF
  DMEF1=DMH*DMI1*LVF
  DMEF2=DMH*DMI2*LVF
  DMEF21=DMH*DMI21*LVF  !! splitting of DMEF2
  DMEF22=DMH*DMI22*LVF   !! splitting of DMEF2

  !----- Toroidal momentum transport

  !Diagonal element of momentum diffusivity
  DMDIAGT=D1*TE**1.5*DMDT/kyrho
  DMT(6)=DMDIAGT

  !Momentum pinches expressed as an effective diffusivity
  DMTEF = rmaj / rmaj0 * &
       ( DMIP  &
       + 2.0*d1*te**1.5*DMIT ) &
       /sign( max(abs(gvt*vtor),0.01) , gvt*vtor )

  !Total momentum pinch expressed as a convection velocity in m/s
  !Diagnostic only, not output in interface
  HPT(6) = ( dmip + 2.0*d1*te**1.5*dmit ) / sign( max( abs(vtor) ,0.01), vtor )

  !Parts of momentum pinches expressed as a convective velocity
  !These are output through the convective velocity array

  !Sum of thermoelectric and turbulent equipartition pinches
  v_p_thr = DMIP / sign( max(abs(rmaj0*vtor), 0.01), rmaj0*vtor )
!  v_p_thr = DMIP / sign( max(abs(vtor), 0.01), vtor )

  !Sum of parallel velocity and Reynolds stress pinches
  v_p_tor = 2.0*d1*te**1.5*DMIT / sign( max(abs(rmaj0*vtor), 0.01), rmaj0*vtor )

  !Total flux expressed as an effective diffusivity
  !WARNING: it is better to use diagonal element + convection
  !         Momentum pinches can produce large negative diffusion
  DMTEFF = DMT(6) + DMTEF

  !----- Poloidal momentum diffusivity
  ! WARNING: These computations have not been tested
  !          (There is no poloidal momentum eq. in PTRANSP,
  !           where the code was implemented and tested)

  !Diagonal element of poloidal momentum diffusivity
  !This should be fine
  DMDIAG=D1*TE**1.5*DMD/kyrho
  DM(4)=DMDIAG

  !Poloidal momentum velocity pinches
  !This looks like the Reynolds source terms
  HPT(4)=0.5D0*DMH*(DMI1+DMI21+DMI22)

  !Poloidal momentum pinches expressed as a diffusivity
  DMEFTEST=-2.*HPT(4)*LVF

  !Total flux expressed as a diffusivity 
  DMEFF=DM(4) -HPT(4)*2.*LVF

  !Unintelligible stuff. Quite likely, the Reynolds stress term
  !is being expressed as a source term here
  SP=1.
  IF(vpol.LT.0.) SP=-1.
  SMEF1=SP*D1*TE**1.5*DMI1/rmaj
  SMEF21=SP*D1*TE**1.5*DMI21/rmaj   !! splitting of DMEF2
  SMEF22=SP*D1*TE**1.5*DMI22/rmaj   !! splitting of DMEF2
!  Vconv(IK)=2.*(SMEF1+SMEF21+SMEF22)*(Csound/CSOUND0)
!  SRE=Vconv(IK)/(rmin*LVF)

  !End of momentum transport computations
  !-----------------------------------------------------------
  !Write effective diffusion and convection to returned arrays
  !

  diffs(1) = SCHI   !Ion heat diffusivity
  diffs(2) = SD     !Electron particle diffusivity?
  diffs(3) = SCHE   !Electron heat diffusivity
  diffs(4) = SDQ    !Impurity particle diffusivity
  diffs(5) = DMT(6) !Toroidal momentum diffusivity
  diffs(6) = DM(4)  !Poloidal momentum diffusivity
                    !Note that with the curvature effect omitted
                    !DM(4) is also the toroidal momentum diffusivity
  conv(1) = 0.0 
  conv(2) = 0.0 
  conv(3) = 0.0 
  conv(4) = v_p_tor ! Sum of parellel velocity + Reynolds stress pinch
  conv(5) = v_p_thr ! Thermoelectric pinch
  conv(6) = hpt(4) ! Sum of poloidal momentum related convection

  return
end subroutine w20diff
