Module modmmm7_1

!                 >>>>>   N O T E   <<<<<
! {{{ and }}} are VIM folding markers. DON NOT REMOVE THEM.

!----- DESCRIPTION -------------------------------------------------{{{
!
! This Fortran 90 module contains a subroutine called mmm and its
! supporting codes. Together these codes provide the capability for
! computing plasma transport coefficients using the Multi-Mode transport
! model (MMM). The current version of MMM includes contributions from
! four transport models: Weiland (W20), Drift-resistive-inertial
! Ballooning Mode (DRIBM or DBM), Horton ETG (ETG) and Paleoclassical
! (PLC).
!
! This module contains two public subroutines:
!
!    * Subroutine mmm
!    * Subroutine set_mmm_switches
!
! Detailed description can be found in each subroutine.
!
! Revision History
! ----------------
! Oct 11, 2010    Lixiang Luo
!                 V7.1 Original Release
!
!-------------------------------------------------------------------}}}

!----- MODULE SPECIFICATIONS ---------------------------------------{{{
Implicit None

Private ! All definitions are assumed private to avoid naming
        ! conflicts with external codes.

!.. Define a better Real*8 type
Integer, Parameter :: R8 = KIND( 0D0 )

!.. Physical constants

Real(R8), Parameter :: &
   zpi    = 3.1415926535898D0 ,&! Pi
   zcc    = 299792458D0       ,&! Speed of light                 [m/sec]
   zcmu0  = 4D-7 * zpi        ,&! Vacuum magnetic permeability   [henrys/m]
   zceps0 = 1D0/(zcc**2*zcmu0),&! Vacuum electrical permittivity [farads/m]
   zckb   = 1.602176487D-16   ,&! Energy conversion factor       [Joule/keV]
   zcme   = 9.10938215D-31    ,&! Electron mass                  [kg]
   zcmp   = 1.672621638D-27   ,&! Proton mass                    [kg]
   zce    = 1.602176487D-19     ! Electron charge                [Coulomb]

!.. Machine Epsilon (smallest number so that 1.0+zepslon>1.0)
Real(R8), Parameter :: &
   zepslon = 2D0**(-53) ,&! As defined by IEEE 754-2008
   zepsqrt = 2D0**(-26)   ! Square root of Epsilon, approximated

!----- PUBLIC ENTRIES ------------------------------------------
Integer, Parameter :: &
   MMM_NCH   = 6, &! Maximum number of transport channels
   MMM_NMODE = 4, &! Number of internal models
   MAXNOPT   = 10  ! Maximum number of internal switches

!.. Identifiers for internal models
Integer, Parameter :: &
   KW20 = 1 ,&! Weiland 20
   KDBM = 2 ,&! DRIBM
   KETG = 3 ,&! Horton ETG
   KPLC = 4   ! Paleoclassical

Public :: mmm7_1
Public :: set_mmm7_1_switches

!!! END OF MODULE SPECIFICATION }}}

Contains

Subroutine mmm7_1 (                                            &
   rmin,    rmaj,    elong,                                    &
   ne,      ni,      nz,      nf,                              &
   xzeff,   te,      ti,      q,       btor,                   &
   zimp,    aimp,    ahyd,    aimass,  wexbs,                  &
   gne,     gni,     gnh,     gnz,     gte,     gti,     gq,   &
   gvrin,   vtorin,  gvpin,   vpolin,  gvparin, vparin,        &
   eta,     csnd0,   dkdeps,                                   &
   thiig,   thdig,   theig,   thzig,   thtig,   thttig,        &
   xkiW20,  xdiW20,  xkeW20,                                   &
   xkiDBM,  xkhDBM,  xkeDBM,                                   &
   xkeETG,  xkePLC,                                            &
   velthi,  vflux,   yvelpa,                                   &
   npoints, lprint,  nprout,  nerr,                            &
   cmodel,  cswitch, lswitch )

!----- DESCRIPTION --------------------------------------------------{{{
!
!  Adjustable switches
!
!    They are passed by lswitch and cswitch. Their minimum size of the
!    first dimension required is 2 in the current version of MMM. However,
!    this is subjected to change in the future versions.
!
!  List
!  -------------------
!  W20   L1    0    Use older EPS2006 model
!        C1   1.0   Momentum pinch scaling factor
!
!  DBM   L1    0    Enable limitation of gradients
!        L2    0    Enable ExB effects
!        C1   0.0   Lower bound of magnetic shear
!
!  ETG   L1    0    Use Jenko threshold
!        C1   1.0   CEES/CEEM scaling factor
!
!--------------------------------------------------------------------}}}

!----- SUBROUTINE SPECIFICATIONS ------------------------------------{{{
Implicit None

! All the following 1-D arrays are assumed to be defined on flux
! surfaces called zone boundaries where the transport fluxes are
! to be computed.  The number of flux surfaces is given by npoints
! (see below).  For example, if you want to compute the transport
! on only one flux surface, set npoints = 1.

Real(R8), Intent(In), Dimension(1:) :: &
   rmin, rmaj, elong, ne, ni, nz, nf, xzeff, &
   te, ti, q, btor, zimp, aimp, ahyd, aimass, wexbs, &
   gne, gni, gnh, gnz, gte, gti, gq
!
! rmin(jz)   = minor radius (half-width) of zone boundary [m]
! rmaj(jz)   = major radius to geometric center of zone boundary [m]
! elong(jz)  = local elongation of zone boundary
!
! ne(jz)     = electron density [m^-3]
! ni(jz)     = sum over thermal hydrogenic ion densities [m^-3]
! nz(jz)     = sum over impurity ion densities [m^-3]
! nf(jz)     = electron density from fast (non-thermal) ions [m^-3]
!
! xzeff(jz)  = Z_eff
! te(jz)     = T_e (electron temperature) [keV]
! ti(jz)     = T_i (temperature of thermal ions) [keV]
! q(jz)      = magnetic q-value
! btor(jz)   = ( R B_tor ) / rmaj(jz)  [tesla]
!
! zimp(jz)   = average density weighted charge of impurities
!            = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
!              sum_imp = sum over impurity ions with charge state Z_imp
!
! aimp(jz)   = average density weighted atomic mass of impurities
!            = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where
!              sum_imp = sum over impurity ions, each with mass M_imp
!
! ahyd(jz)   = average density weighted atomic mass of hydrogen ions
!            = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
!              sum_hyd = sum over hydrogenic ions, each with mass M_hyd
!
! aimass(jz) = mean atomic mass of thermal ions [AMU]
!            = ( sum_i n_i M_i ) / ( sum_i n_i ) where
!              sum_i = sum over all ions, each with mass M_i
!
! wexbs(jz)  = ExB shearing rate in [rad/s].  See  K.H. Burrell,
!              "Effects of {ExB} velocity shear and magnetic shear
!              on turbulence and transport in magnetic confinement
!              devices", Phys. of Plasmas, 4, 1499 (1997).
!
!   All of the following normalized gradients are at zone boundaries.
!   r = half-width, R = major radius to center of flux surface
!
! gne(jz) = -R ( d n_e / d r ) / n_e
! gni(jz) = -R ( d n_i / d r ) / n_i
! gnh(jz) = -R ( d n_h / d r ) / n_h
! gnz(jz) = -R ( d Z n_Z / d r ) / ( Z n_Z )
! gte(jz) = -R ( d T_e / d r ) / T_e
! gti(jz) = -R ( d T_i / d r ) / T_i
! gq(jz)  =  R ( d q   / d r ) / q    related to magnetic shear
!
! where:
!   n_i  = thermal ion density (sum over hydrogenic and impurity)
!   n_h  = thermal hydrogenic density (sum over hydrogenic species)
!   n_Z  = thermal impurity density,  Z = average impurity charge
!                     summed over all impurities

!.. Variables related to momentum transport in 2006 Weiland model
Real(R8), Intent(In), Optional, Dimension(1:) :: &
   gvrin, vtorin, gvpin, vpolin, gvparin, vparin
!
! gvrin(jz)  = Normalized toroidal velocity gradient (r/v_tor)*(dv_tor/dr)
! vtorin(jz) = Toroidal velocity [m/s]
! gvpin(jz)  = Normalized poloidal velocity gradient (r/v_pol)*(dv_pol/dr)
! vpolin(jz) = Poloidal velocity [m/s]
! gvparin(jz)= Normalized parallel velocity gradient (r/v_par)*(dv_par/dr)
! vparin(jz) = Parallel velocity [m/s]

Real(R8), Intent(In), Optional, Dimension(1:) :: &
   eta
!
!  eta(jz) = Spitzer/Neoclassical resistivity [Ohm-m]
!
! This is an optional argument:
! > If not associated:
!   Spitzer resistivity will be calculated internally and used for both
!   DRIBM and Paleoclassical.
! > If associated:
!   eta will be used as the resistivity for Paleoclassical. Spitzer
!   resistivity will be calculated internally and used for DRIBM.

Real(R8), Intent(In), Optional :: &
   csnd0 ! Sound speed at the magnetic axis
! This is an optional argument:
! > If not associated:
!   The sound speed at the magnetic axis will be calculated internally
!   using other profiles at the magnetic axis.
! > If associated:
!   csnd0 will be used as the sound speed at the magnetic axis

Real(R8), Intent(In), Optional, Dimension(1:) :: &
   dkdeps ! Elongation gradient
!
!  dkdeps = d elong / d r*
!    where
!  elong is elongation and r*(jz)=rmin(jz)/rmaj(jz)
!
! > If not associated:
!   The elongation gradient will be calculated internally using other profiles.
! > If associated:
!   The elongation gradient will be assigned by dkdeps.

! Some output arguments are optional. The calling program is not required to
! use them. When they are associated, the actual arguments have to be
! allocated with enough space to store all return values (npoints is the
! minimal dimension) in advance.

Real(R8), Intent(Out), Dimension(1:) :: & ![m^2/s]
   thiig, thdig, theig, thzig, thtig, thttig
!
!  thiig(jz) = ion thermal diffusivity from the Weiland model
!  thdig(jz) = hydrogenic ion diffusivity from the Weiland model
!  theig(jz) = electron thermal diffusivity from the Weiland model
!  thzig(jz) = impurity ion diffusivity from the Weiland model
!  thtig(jz) = Toroidal momentum diffusivity
!  thttig(jz)= Poloidal momentum diffusivity

! The following component output arrays give the separate contribution from
! each internal model. Note that the momentum diffusivities are only provided
! by the Weiland model. Generally, these arrays are used for diagnostic
! output only.
!
Real(R8), Intent(Out), Optional, Dimension(1:) :: & ![m^2/s]
   xkiW20, xdiW20, xkeW20, xkiDBM, xkhDBM, xkeDBM, xkeETG, xkePLC
!
! xkiW20(jz) = W20 ion thermal diffusivity
! xdiW20(jz) = W20 particle diffusivity
! xkeW20(jz) = W20 electron thermal diffusivity
! xkiDBM(jz) = DBM ion thermal diffusivity
! xkhDBM(jz) = DBM hydrogenic ion diffusivity
! xkeDBM(jz) = DBM electron thermal diffusivity
! xkeETG(jz) = ETG electron thermal diffusivity
! xkePLC(jz) = PLC electron thermal diffusivity
!
! The component output arrays are optional.

Real(R8), Intent(Out), Dimension(1:,1:) :: &
   velthi, &! Convective velocities [m/s]
   vflux    ! Flux matrix

Real(R8), Intent(Out), Optional, Dimension(1:) :: &
   yvelpa ! Electron thermal heat pinch [m/s]
! This argument is optional.

Integer, Intent(In) :: &
   npoints ! Number of values in all of the 1-D arrays

Integer, Intent(In) :: &
   lprint, &! Verbose level
   nprout   ! Output unit number for long printout

Integer, Intent(Out) :: &
   nerr ! Error return value

Real(R8), Intent(In), Optional, Dimension(1:) :: &
   cmodel ! Weights of internal models
!
! This is an optional argument. Only the first four elements are used.
! > If associated:
!   cmodel(1)~cmodel(4) are assigned as the weights for Weiland20, DRIBM,
!   ETG and Paleoclassical, respectively.
! > If not associated:
!   Equivalent to cmodel=(/1.0,1.0,1.0,1.0/)

Real(R8), Intent(In), Optional, Dimension(1:,1:) :: &
   cswitch ! Real adjustable parameters
!
! This is an optional argument. The second dimension is corresponding to
! the index of an internal model, and the first dimension to the index
! of a real adjustable parameter for that particular model. For example,
! the first real adjustable parameter for the Weiland model would be
! stored in cswitch(1,KW20) or cswitch(1,1).
!
! > If associated:
!   The internal parameters are assigned to the given values of the
!   actual argument.
! > If not associated:
!   all actual values default to the internally set values.
!
! An up-to-date list of parameters can be found at the header of
! this subroutine, along with the default values.

Integer, Intent(In), Optional, Dimension(1:,1:) :: &
   lswitch ! Integral adjustable parameters
!
! This is an optional argument. The second dimension is corresponding to
! the index of an internal model, and the first dimension to the index
! of a integral adjustable parameter for that particular model. For example,
! the second integral adjustable parameter for the DRIBM model would be
! stored in cswitch(2,KDBM) or cswitch(2,2).
!
! > If associated:
!   The internal switches are assigned to the given values of the
!   actual argument.
! > If not associated:
!   all actual values default to the internally set values.
!
! An up-to-date list of parameters can be found at the header of
! this subroutine, along with the default values.

!----- LOCAL VARIABLES  ---------------------------------------------------

Real(R8) :: &
   cscal(MMM_NMODE)! Weights for individual models

Real(R8) :: &
   csw(MAXNOPT,MMM_NMODE) ! Real adjustable parameters

Integer :: &
   lsw(MAXNOPT,MMM_NMODE) ! Integral adjustable parameters

Integer :: &
   jz ! Loop counter

Real(R8) :: &
   zep,    &! Inverse aspect ratio
   zgyrfe, &! Electron gyrofrequency
   zgyrfi, &! Ion gyrofrequency
   zbeta,  &! Thermal/magnetic energy ratio
   zvthe,  &! Thermal velocity of electrons
   zvthi,  &! Thermal velocity of thermal ions
   zsound, &! Speed of sound
   zsound_axis,&! Speed of sound at the magnetic center
   zlog,   &! Coloumb logorithm
   zcf,    &! A factor in collision frequency
   zlare,  &! Electron Larmor radius
   zlari,  &! Ion Larmor radius
   zlarpo, &! Poloidal ion Larmor radius
   zrhos,  &! Ion larmor radius with Te
   zwn,    &! 0.3/rhos
   zgmax,  &! Upper bound a gradients
   znuei,  &! Electron collision frequency
   znude,  &! Electron magnetic drift frequency
            ! Also used as the normalization factor in W20 and DBM
   znuhat   ! Normalized electron collision frequency

Real(R8) :: &
   eta_general ! General resistivity, calculated internally or given by eta

!.. Weiland variables{{{
Real(R8) :: &
   znueff ! Effective electron collision frequency

Real(R8), Dimension(MMM_NCH) :: &
   zchieff, &! Diffusivitiy return values
   zvconv    ! Effective convective velocity return values

Real(R8) :: &
   zthte,  &! Ti/Te
   zbetae, &! Thermal/magnetic energy ratio for electrons
   ztzte,  &! Tz/Te
   zfnzne, &! n_z/n_e
   zfnsne, &! n_f/n_e, Nf is the density of superthermal ions
   zftrap, &! Fraction of trapped particles
   zkyrho, &! k_y
   znormd, &! Factor to convert normalized diffusivities
   znormv, &! Factor to convert normalized convective velocities
   shear    ! Magnetic shear

Real(R8) :: & ! Copies of gne, gnh, gte and gth
   zgne, zgnh, zgte, zgth !}}}

!.. Horton ETG variables {{{
Real(R8) :: &
   cees,  &! Scaling factor for electrostatic case
   ceem,  &! Scaling factor for electromagnetic case
   cl,    &! 3/2*sqrt(pi/2)
   zwpe,  &! Electron plasma frequency
   zlce,  &! Mixing length of ETG mode
   zdeltae ! Collisionless skin depth

Real(R8) :: &
   zgtec,           &! Critical electron temperature gradient
   zgte_thr_horton, &! in Horton's original model
   zgte_thr_jenko    ! in Jenko's modified model

Real(R8), Dimension(npoints) :: &
   zepa,   &! rmin/rmaj
   dk_deps,&! rate of change of elongation wrt epsilon
   theeg,  &! Electron thermal diffusivity
   thees,  &! Electron thermal diffusivity for electrostatic cases
   theem    ! Electron thermal diffusivity for electromagnetic cases

Logical :: nlthr ! Switch for using the Jenko model
!}}}

!.. DRIBM variables {{{
Integer :: nmodes ! Number of DRIBM unstable modes

Real(R8) :: &
   wexb_drbm,&! ExB shearing rate for DRIBM
   difftemp   ! Temporary storage
!}}}

!.. Paleoclassical variables {{{

REAL(R8), Dimension(npoints):: &
   ychepa ! Electron thermal diffusivity by PLC

Real(R8), Dimension(npoints) :: &
   z_L_ratio  ! L / ( pi * R * q )

Real(R8) :: &
   zlambda_e, &! v_{the} / \nu_e
   z_l,       &! min ( zlambda_e, z_pi R q n_max )
   z_n_max,   &! n_max
   zln_lambda  ! ln( \Lambda / 17)

Real(R8), Parameter :: &
   z_scale = 1D0 ! scale factor in front of transport

Real(R8) :: &
   z_deltae,  &! normalized electromagnetic skin depth
   z_q_prime, &! dq/dr
   z_lmax,    &! length over which magnetic field diffuse radially
   z_omegape   ! plasma frequency
!}}}
!}}}

!=== EXECUTION BODY ====================================================

!.. Initialize arrays {{{

thiig  = 0D0
thdig  = 0D0
theig  = 0D0
thzig  = 0D0
thtig  = 0D0
thttig = 0D0
theeg  = 0D0
thees  = 0D0
theem  = 0D0
velthi = 0D0
vflux  = 0D0
zchieff= 0D0
zvconv = 0D0

nerr = 0

!.. Assign adjustable parameters accoding to optional arguments or
!   default values if the optional argument is omitted
If ( Present( cmodel ) ) Then
   cscal(1:MMM_NMODE) = cmodel(1:MMM_NMODE)
Else
   cscal(1:MMM_NMODE) = 1D0
End If
!
csw = 0D0
If ( Present( cswitch ) ) Then
   jz = min( MAXNOPT, size(cswitch,1) )
   csw(1:jz,1:MMM_NMODE) = cswitch(1:jz,1:MMM_NMODE)
Else
   csw(1,KW20) = 1D0
   csw(1,KETG) = 1D0
End If
!
lsw = 0
If ( Present( lswitch ) ) Then
   jz = min( MAXNOPT, size(lswitch,1) )
   lsw(1:jz,1:MMM_NMODE) = lswitch(1:jz,1:MMM_NMODE)
End If

!.. Diffusivities components
If ( Present( xkiW20 ) ) xkiW20 = 0D0
If ( Present( xdiW20 ) ) xdiW20 = 0D0
If ( Present( xkeW20 ) ) xkeW20 = 0D0
If ( Present( xkiDBM ) ) xkiDBM = 0D0
If ( Present( xkhDBM ) ) xkhDBM = 0D0
If ( Present( xkeDBM ) ) xkeDBM = 0D0
If ( Present( xkeETG ) ) xkeETG = 0D0
If ( Present( xkePLC ) ) xkePLC = 0D0
If ( Present( yvelpa ) ) yvelpa = 0D0

!.. Calculate some ETG variables only if ETG is included
If ( Present( dkdeps ) ) Then
   dk_deps(1:npoints) = dkdeps(1:npoints)
Else
   If ( cscal(KETG) > 0D0 .and. npoints >= 3 ) Then
      zepa(1:npoints) = rmin(1:npoints) / rmaj(1:npoints)
      do jz = 2, npoints-1
         dk_deps(jz) = ( elong(jz+1) - elong(jz-1) ) / &
            ( 1.0E-6 + zepa(jz+1) - zepa(jz-1) )
      end do
      dk_deps(1)=dk_deps(2)
      dk_deps(npoints)=dk_deps(npoints-1)
   End If
End If
!}}}

!.. Start the main do-loop over the radial index "jz"
MAINLOOP: &
Do jz = 1, npoints

   ! Avoid computation of fluxes at rmin(jz) < 1e-4 * rmaj(jz)
   If (rmin(jz) < (1D-4 * rmaj(jz))) cycle

   ! Calculation of physical quantities {{{

   shear = gq(jz) * rmin(jz) / rmaj(jz)

   !  compute inverse aspect ratio
   zep    = max( rmin(jz)/rmaj(jz), zepslon )

   zgyrfe = zce * btor(jz) / zcme  ! electron plasma frequency
   zgyrfi = zce * btor(jz) / (zcmp * aimass(jz))
   zbeta  = (2D0* zcmu0 * zckb / btor(jz)**2) * (ne(jz) * te(jz) + ni(jz) * ti(jz))
   zvthe  = sqrt(2D0* zckb * te(jz) / zcme)
   zvthi  = sqrt(2D0* zckb * ti(jz) / (zcmp * aimass(jz)))
   zsound = sqrt(zckb * te(jz) / (zcmp * aimass(jz)))
   If ( Present(csnd0) ) Then
      zsound_axis = csnd0
   Else
      zsound_axis = sqrt(zckb * te(1) / (zcmp * aimass(jz)))
   End If

   zrhos  = zsound / zgyrfi

   zlog   = 37.8D0-log(sqrt(ne(jz)) / te(jz))
   zcf    = (4D0* sqrt(zpi) / 3D0) * &
             (zce / (4D0* zpi * zceps0))**2 * &
             (zce / zckb) * sqrt( (zce/zcme) * (zce/zckb) )
   znuei  = zcf*sqrt(2D0)* ne(jz) * zlog * xzeff(jz) / (te(jz) * sqrt(te(jz)))

   ! General resistivity, used in PLC
   If ( Present(eta) ) Then
      eta_general = eta(jz)
   Else
      eta_general = znuei*0.51D0*zcme/(ne(jz)*zce**2)
   End If

   zlari  = zvthi / zgyrfi
   zlare  = zvthe / zgyrfe
   zlarpo = max(zlari * q(jz) / zep, zepslon)

   zwn    = 0.3D0 / zrhos
   znude  = 2D0 * zwn * zrhos * zsound / rmaj(jz)

   zgmax = rmaj(jz) / zlarpo

   !.. Hydrogen species
   zthte  = ti(jz) / te(jz)

   zbetae = 2D0 * zcmu0 * zckb * ne(jz) * te(jz) / btor(jz)**2

   !.. Impurity species (use only impurity species 1 for now)
   !  assume T_Z = T_H throughout the plasma here

   !Ratio of impurity temperature to electron temperature
   ztzte  = ti(jz) / te(jz)

   !Ratio of impurity density to electron density
   zfnzne = nz(jz) / ne(jz)

   !Fraction of superthermal ions (beam, RF minority) to electrons
   zfnsne = nf(jz) / ne(jz)

   !Fraction of trapped particles
   zftrap = sqrt( 2D0* rmin(jz) / ( rmaj(jz) * ( 1D0+ rmin(jz) / rmaj(jz) )))

   !Default k_y rho for Weiland model
   zkyrho = 0.316D0

   ! End calculation of physical quantities!}}}

   ! Weiland20 {{{
   If ( cscal(KW20) > 0D0 ) Then

      call w20main ( lprint, &
           te(jz), ne(jz), vtorin(jz), vpolin(jz), vparin(jz), &
           btor(jz), rmin(npoints), rmin(jz), rmaj(jz), rmin(npoints)/rmaj(1), &
           aimp(jz), ahyd(jz), zimp(jz), &
           gte(jz), gti(jz), gti(jz), &
           gne(jz), gni(jz), gnz(jz), &
           gvrin(jz), gvpin(jz), gvparin(jz), &
           zkyrho, zthte, ztzte, zfnzne, zfnsne, zftrap, &
           q(jz), shear, elong(jz), wexbs(jz), &
           zsound_axis, &
           zchieff, zvconv )

      !.. If nerr not equal to 0 an error has occured
      If (nerr /= 0) Return

      !  Conversion factors for diffusion and convective velocity
      znormd = 2D0 * zsound * zrhos**2 / ( rmaj(jz) * zkyrho )
      znormv = 2D0 * zsound * zrhos**2 / ( rmaj(jz)**2 * zkyrho )

      !.. Weiland's own diffusivities
      If ( Present( xkiW20 ) ) xkiW20(jz) = zchieff(1)
      If ( Present( xdiW20 ) ) xdiW20(jz) = zchieff(2)
      If ( Present( xkeW20 ) ) xkeW20(jz) = zchieff(3)

      !  Store effective diffusivities from Weiland19 model
      thiig(jz) = cscal(KW20) * zchieff(1)
      thdig(jz) = cscal(KW20) * zchieff(2)
      theig(jz) = cscal(KW20) * zchieff(3)
      thzig(jz) = cscal(KW20) * zchieff(4)
      thtig(jz) = cscal(KW20) * zchieff(5)
      thttig(jz)= cscal(KW20) * zchieff(6)

      !  start computing the fluxes
      vflux(1,jz) = cscal(KW20) * thiig(jz) * gti(jz) / rmaj(jz)
      vflux(2,jz) = cscal(KW20) * thdig(jz) * gnh(jz) / rmaj(jz)
      vflux(3,jz) = cscal(KW20) * theig(jz) * gte(jz) / rmaj(jz)
      vflux(4,jz) = cscal(KW20) * thzig(jz) * gnz(jz) / rmaj(jz)

      !  Add momentum transport convective velocities from Reynold's source term
      !  in Weiland model and other pinches
      velthi(1,jz) = cscal(KW20)*zvconv(1)
      velthi(2,jz) = cscal(KW20)*zvconv(2)
      velthi(3,jz) = cscal(KW20)*zvconv(3)
      velthi(4,jz) = cscal(KW20)*csw(1,KW20)*zvconv(4)
      velthi(5,jz) = cscal(KW20)*csw(1,KW20)*zvconv(5)
      velthi(6,jz) = cscal(KW20)*csw(1,KW20)*zvconv(6)

      !Option 1: Use Weiland's 2006 EPS model
      If ( lsw(1,KW20)==1 ) Then
         velthi(5,jz) = velthi(4,jz)
      Else
         velthi(5,jz) = velthi(4,jz) + velthi(5,jz)
      End If

   End If ! End of Weiland model!}}}

   ! DRIBM {{{
   If ( cscal(KDBM) > 0D0 ) Then

      If ( lsw(1,KDBM) > 0 ) Then
         zgmax = rmaj(jz) / zlarpo
         zgne = sign( min( abs( gne(jz) ), zgmax ), gne(jz) )
         zgnh = sign( min( abs( gnh(jz) ), zgmax ), gnh(jz) )
         zgte = sign( min( abs( gte(jz) ), zgmax ), gte(jz) )
         zgth = sign( min( abs( gti(jz) ), zgmax ), gti(jz) )
      Else
         zgne = gne(jz)
         zgnh = gnh(jz)
         zgth = gti(jz)
         zgte = gte(jz)
      End If

      zkyrho = 0.2D0
      zwn    = 0.2D0 / zrhos
      znude  = 2D0 * zwn * zrhos * zsound / rmaj(jz)
      znuhat = znuei / znude

      If ( lsw(2,KDBM)>0 ) Then
         wexb_drbm=wexbs(jz) ! ExB effect included
      Else
         wexb_drbm=0 ! Ignore ExB effect
      End If

      call drbm ( zgne,   zgnh,   zgte,   zgth, &
         zthte,   zbetae, znuhat, q(jz),  zkyrho, wexb_drbm, &
         zchieff, nmodes, nerr )

      !.. Conversion factors for diffusion and convective velocity
      znormd = 2D0 * zsound * zrhos**2 / ( rmaj(jz) * zkyrho )
      znormv = 2D0 * zsound * zrhos**2 / ( rmaj(jz)**2 * zkyrho )

      difftemp = cscal(KDBM)*znormd* max(0D0,zchieff(1))
      If ( Present( xkiDBM ) ) xkiDBM(jz) = difftemp
      thiig(jz) = thiig(jz) + difftemp

      difftemp = cscal(KDBM)*znormd* max(0D0,zchieff(2))
      If ( Present( xkhDBM ) ) xkhDBM(jz) = difftemp
      thdig(jz) = thdig(jz) + difftemp

      difftemp = cscal(KDBM)*znormd* max(0D0,zchieff(3))
      If ( Present( xkeDBM ) ) xkeDBM(jz) = difftemp
      theig(jz) = theig(jz) + difftemp

      !.. Start computing the fluxes
      vflux(1,jz) = vflux(1,jz) + thiig(jz) * zgth / rmaj(jz)
      vflux(2,jz) = vflux(2,jz) + thdig(jz) * zgnh / rmaj(jz)
      vflux(3,jz) = vflux(3,jz) + theig(jz) * zgte / rmaj(jz)

   End If!}}}

   ! Horton ETG{{{
   If ( cscal(KETG) > 0D0 ) Then
      !Adjustable constants
      cees  = 0.06D0 * csw(1,KETG) !0.06
      ceem  = 0.06D0 * csw(1,KETG) !0.06
      nlthr = lsw(1,KETG) > 0 ! [SWITCH]

      cl = 1.88D0

      !Thermal velocity of electrons [ m / s ]
      zvthe = sqrt( zckb * te(jz) / zcme )

      !Electron gyrofrequency [ 1 / s ]
      zgyrfe = zce * btor(jz) / zcme

      !Electron gyroradius [ m ]
      zlare = zvthe / zgyrfe

      !Electron plasma frequency [ 1 / s ]
      zwpe = sqrt( (ne(jz)*zce**2) / (zcme * zceps0) )

      !Mixing length of ETG mode
      zlce = max( 0D0, q(jz) * zlare * gte(jz))

      !Collisionless skin depth
      zdeltae = zcc / zwpe

      !Critical electron temperature gradient in Horton's paper
      zgte_thr_horton = cl * &
           abs( shear ) / q(jz) * &
           (1. + xzeff(jz) * te(jz) / ti(jz) )
      !Critical electron temperature gradient by Jenko ETG model
      zgte_thr_jenko = max( ( 1D0 + xzeff(jz) * te(jz) / ti(jz) ) * &
           ( 1.33D0+1.91D0*abs( shear ) / q(jz) ) * &
           ( 1D0-1.5D0*zep ) * ( 1D0 + 0.3D0 * zep * dk_deps(jz) ), &
           0.8D0*gne(jz) )

      If ( nlthr ) Then
         zgtec = zgte_thr_jenko
      Else
         zgtec = zgte_thr_horton
      End If

      thees(jz) = cees * (q(jz) * zlare)**2 * zvthe &
              * gte(jz)* sqrt( abs(gte(jz)) ) / rmaj(jz) * &
              max( 0D0, ( gte(jz)- zgtec ) )

      theem(jz) = ceem * zdeltae**2 * zvthe &
              * sqrt( max( 1D-6, gte(jz)) ) / rmaj(jz)

      If ( zlce > zdeltae ) then
         theeg(jz) = thees(jz)
      Else if ( zlce < zdeltae ) Then
         theeg(jz) = theem(jz)
      End If

      theeg(jz) = min(1D2, theeg(jz))

      !Enforce Jenko ETG thresold
      If ( nlthr ) Then

         If ( gte(jz)< zgte_thr_jenko ) Then
            theeg(jz) = 0D0
         Else
            If ( zlce > zdeltae ) Then
               !Electrostatic case: threshold already applied
               !theeg(jz) = theeg(jz)
            Else
               !Electromagnetic case: apply threshold
               theeg(jz) = tanh( ( gte(jz)- zgte_thr_jenko ) ) * theeg(jz)
            End If
            !Apply bound on diffusivity
            theeg(jz) = min( theeg(jz), 1D2 )
         End If
      End If

      !.. ETG's own diffusivity
      If ( Present( xkeETG ) ) xkeETG(jz) = theeg(jz)

      theig(jz) = theig(jz) + cscal(KETG) * theeg(jz)

      !thr = zgte_thr_jenko
   End If!}}}

   ! Paleoclassical {{{
   If ( cscal(KPLC) > 0D0 ) Then

      zln_lambda = zlog / 17D0
      z_q_prime=gq(jz)*q(jz)/rmaj(jz)

      zlambda_e  = 1.2D16 * ( 1.0D3 * TE(jz) )**2 / &
           ( NE(jz) * xzeff(jz) * zln_lambda )

      ! Electron plasma frequency
      z_omegape = 1.7835D11*sqrt(ne(jz)*1D-19)

      ! Normalized EM skin depth
      z_deltae = zcc / ( z_omegape*rmin(jz) )
      z_q_prime = sign( max( 0.01D0, abs(z_q_prime) ), z_q_prime )
      z_n_max = 1D0/sqrt(zpi*z_deltae*abs(z_q_prime))

      ! For now, just use z_L = zlambda_e
      !     z_L = min ( zlambda_e, zpi * rmajor * z_n_max * q(jz) )

      z_lmax = zpi * rmaj(jz) * z_n_max * q(jz)

      z_L = zlambda_e*z_lmax/(zlambda_e+z_lmax)

      z_L_ratio(jz) = z_L / ( q(jz) * zpi * rmaj(jz) )

      If ( jz == npoints ) Then
         z_L_ratio(npoints) = z_L_ratio(npoints-1)
      End If

      ! General resistivity is used here
      ychepa(jz) = z_scale * 1.5D0 * (z_L_ratio(jz)+1D0) * eta_general / zcmu0

      If ( Present( yvelpa ) ) Then
         yvelpa(jz) = z_L_ratio(jz)
         yvelpa(1) = 0D0
      End If

      !.. PLC's own diffusivity
      If ( Present( xkePLC ) ) xkePLC(jz) = ychepa(jz)

      theig(jz) = theig(jz) + cscal(KPLC) * ychepa(jz)

   End If!}}}

End Do MAINLOOP ! End of the main do-loop over the radial index, "jz"

End Subroutine mmm7_1

Subroutine set_mmm7_1_switches & !{{{
   ( cmmm, lmmm, &
     KW20_C_MOM_PINCH_SCALE, KW20_L_EPS2006, &
     KDBM_C_SHEAR_LBOUND, KDBM_L_GRND_LIMIT, KDBM_L_EXB, &
     KETG_C_CEE_SCALE, KETG_L_NLTHR )

Implicit None

Real(R8), Intent(Out), Dimension(:,:) :: cmmm
Integer, Intent(Out), Dimension(:,:) :: lmmm

Real(R8), Optional, Intent(In) :: &
   KW20_C_MOM_PINCH_SCALE, &! W20 C1
   KDBM_C_SHEAR_LBOUND,    &! DBM C1
   KETG_C_CEE_SCALE         ! ETG C1

Integer, Optional, Intent(In) :: &
   KW20_L_EPS2006,         &! W20 L1
   KDBM_L_GRND_LIMIT,      &! DBM L1
   KDBM_L_EXB,             &! DBM L2
   KETG_L_NLTHR             ! ETG L1

!------- Execution Body --------------------------------------------------

cmmm = 0D0
lmmm = 0

cmmm(1,KW20) = 1D0
If (Present(KW20_C_MOM_PINCH_SCALE)) cmmm(1,KW20) = KW20_C_MOM_PINCH_SCALE

lmmm(1,KW20) = 0
If (Present(KW20_L_EPS2006)) lmmm(1,KW20) = KW20_L_EPS2006

cmmm(1,KDBM) = 0D0
If (Present(KDBM_C_SHEAR_LBOUND)) cmmm(1,KDBM) = KDBM_C_SHEAR_LBOUND
lmmm(1,KDBM) = 0
If (Present(KDBM_L_GRND_LIMIT)) lmmm(1,KDBM) = KDBM_L_GRND_LIMIT
lmmm(2,KDBM) = 0
If (Present(KDBM_L_EXB)) lmmm(2,KDBM) = KDBM_L_EXB

cmmm(1,KETG) = 1D0
If (Present(KETG_C_CEE_SCALE)) cmmm(1,KETG) = KETG_C_CEE_SCALE

lmmm(1,KETG) = 0
If (Present(KETG_L_NLTHR)) lmmm(1,KETG) = KETG_L_NLTHR

End Subroutine set_mmm7_1_switches!}}}

Subroutine drbm & !{{{
   (  gnein,  gnhin,   gtein, gthin,             &
      tauhin, betaein, vefin, q, ekyrhoin, wexb, &
      chieff, nmodes,  nerr  )

! SUBROUTINE SPECIFICATIONS {{{
Implicit None

!----- CONSTANTS ---------------------------------------------------

Integer, Parameter :: &
   neq = 6 ! Number of equations

!----- ARGUMENTS ---------------------------------------------------

Real(R8), Intent(In) :: &
   gnein,  gnhin,  gtein, gthin, &
   tauhin, betaein, q, vefin, ekyrhoin, wexb

Integer, Intent(Out) :: &
   nmodes, & ! nmodes = number of unstable modes
   nerr      ! Error value

Real(R8), Intent(Out) :: &
   chieff(:)

!----- LOCAL VARIABLES  --------------------------------------------

Integer :: j1, j2, j

Real(R8) :: &
   vef, zflh, zgamax, zflxph, zflxnh, &
   zflxpe, zphsph, zphsnh, zphspe, ztemp1, &
   zreal, zimag, zerreal, zerimag, zerrmax

!.. Arrays for eigenvalue problem
! ( zamr(i,j), zami(i,j) ) = matrix A
! ( zbmr(i,j), zbmi(i,j) ) = matrix B
!
! Note that the eigenvalues are
!
!   zomega(j) = zalfr(j) / zbeta(j)
!   zgamma(j) = zalfi(j) / zbeta(j)
!
! and the eigenvectors are
!
!   zevec(j) = cmplx( zvr(j), zvi(j) )
! Here, zbeta(j) will be 0.0 in the case of an infinite eigenvalue
!
! Here, zbeta(j) will be 0.0 in the case of an infinite eigenvalue
Real(R8), Dimension(neq,neq) :: &
   zamr, zami, zbmr, zbmi, &
   zamrt, zamit, zbmrt, zbmit, &
   zvr, zvi

Real(R8), Dimension(neq) :: &
   zomega, zgamma, &
   zalfr, zalfi, zbeta, &
   ztempa, ztempb

Complex*16 :: &
   zevec(neq,neq) = 0D0

!  These are the thermal and particle fluxes and effective diffusivities
!  Normalized transport fluxes         Eff. Diffusivities
!       zflxm(1)   n_H T_H             zchim(1)    chi_i
!       zflxm(2)   n_H                 zchim(2)    D_i
!       zflxm(3)   n_e T_e             zchim(3)    chi_e
Real(R8) :: zflxm(3), zchim(3)

Real(R8) :: &
   zkpsh, zalp, zalf

Integer :: ifail

! END SUBROUTINE SPECIFICATIONS }}}

!.. Initialize variables  {{{
!
! CAUTION: Initial values during definition may not have the expected
! effects, because all local variables have an implicit "SAVE"
! attribute in Fortran 90.
!
zomega = 0D0
zgamma = 0D0
zalfr  = 0D0
zalfi  = 0D0
zbeta  = 0D0
zflxm  = 0D0
zchim  = 0D0
zami   = 0D0
zbmr   = 0D0
zevec  = ( 0D0, 0D0 )
zamr   = 0D0
zami   = 0D0
zbmr   = 0D0
zbmi   = 0D0
zamrt  = 0D0
zamit  = 0D0
zbmrt  = 0D0
zbmit  = 0D0
zvr    = 0D0
zvi    = 0D0
nerr = 0

vef = 1D0 * vefin ! Possible scaling

!.. End of initialization

!.. compute the rest of the dimensionless variables needed
!     zflh  = k_y^2 \rho_{sH}^2

zflh = ekyrhoin**2

If ( gnhin + gthin < zepsqrt ) then
   zalp  = 0D0
   zalf  = 0D0
   zkpsh = 0D0
Else
   zalp = 1D0
   zkpsh = 1D0 * 0.5D0 * sqrt( zalp / zflh ) / q ![eric]scaling
End If
!}}}

!.. Equations for e \phi / T_e, T_H, n_H, T_e, A// and v// {{{

!.. Hydrogen density
zamr(1,1) = - 1D0 + 0.5D0 * (gnhin - zflh * tauhin * (gnhin+gthin))
zamr(1,2) = - tauhin
zamr(1,3) = - tauhin
zamr(1,6) = zkpsh
zbmr(1,1) = zflh
zbmr(1,3) = 1D0

!.. Total momentum equation
zamr(2,2) = zkpsh*tauhin
zamr(2,3) = zkpsh*(1D0+tauhin)
zamr(2,4) = zkpsh
zamr(2,5) = -0.5D0 * (tauhin*(gthin+gnhin) + gtein+gnein)
zbmr(2,6) = 1D0

!.. Hydrogen energy
zamr(3,1) = -0.5D0 * ( gthin - (2D0/3D0)*gnhin )
zamr(3,2) = tauhin * (5D0/3D0)
zbmr(3,2) = - 1D0
zbmr(3,3) =  (2D0/3D0)

!.. Electron energy
zamr(4,1) =  0.5D0 * ( gtein - (2D0/3D0)*gnein )
zamr(4,4) =  (5D0/3D0)
zamr(4,5) = - zkpsh*0.96D0*zflh / betaein
zami(4,4) = - 3600D0*zkpsh**2 / vef
zami(4,5) =   918.0D0*zkpsh*gtein / vef
zbmr(4,3) = - (2D0/3D0)
zbmr(4,4) = 1D0

!.. Vorticity equation
zamr(5,1) = - 0.5D0*zflh*tauhin*(gnhin+gthin)
zamr(5,2) = - tauhin
zamr(5,3) = - (1D0 + tauhin)
zamr(5,4) = - 1D0
zamr(5,5) = 2D0*zkpsh*zflh / betaein
zbmr(5,1) = zflh

!.. Ohms law
zamr(6,1) = zkpsh
zamr(6,3) = - zkpsh
zamr(6,4) = - 1.71D0*zkpsh
zamr(6,5) = 0.5D0*(gnein + 1.71D0*gtein)
zami(6,5) = - 0.0054D0 * vef * zflh / betaein
zbmr(6,5) = 1D0 + 2D0*0.0054D0*zflh / betaein
zbmr(6,6) = - 0.0054D0

!.. Find the eigenvalues and eigenvectors using ACM/TOMS routine 535
Ifail = -1

!.. Update A matrix w.r.t wexb and save local copies of the matrices
Do j1=1,neq
   Do j2=1,neq
      zamr(j1,j2)  = zamr(j1,j2) + abs(wexb)*zbmi(j1,j2)
      zami(j1,j2)  = zami(j1,j2) - abs(wexb)*zbmr(j1,j2)
      zamrt(j1,j2) = zamr(j1,j2)
      zamit(j1,j2) = zami(j1,j2)
      zbmrt(j1,j2) = zbmr(j1,j2)
      zbmit(j1,j2) = zbmi(j1,j2)
   End Do
End Do
!}}}

!.. Solution of the generalized eigenvalue problem {{{
Call r8tomsqz(neq,neq,zamr,zami,zbmr,zbmi,zalfr,zalfi,zbeta,zvr,zvi,ifail)
!
If ( ifail /= 0 ) Then
   nerr = ifail
   Return
End If!}}}

!.. Storing the results -  eigenvalues and eigenvectors {{{
zgamax = 0D0
Do j=1,neq
   ztemp1 = zbeta(j)
   If ( abs(zbeta(j)) < zepsqrt ) ztemp1 = zepsqrt
   zomega(j) = zalfr(j)  / ztemp1
   zgamma(j) = zalfi(j)  / ztemp1
   zgamax = max ( zgamax, zgamma(j) )
   Do j1=1,neq
      zevec(j1,j) = DCMPLX ( zvr(j1,j), zvi(j1,j))
   End Do
End Do !}}}

!.. Calculation of diffusivities and fluxes if transport exists {{{
If_HAVETXP: &
If ( zgamax > zepsqrt ) Then ! NO -> No transport, go to finalization

!  Real and imaginary parts
   zerrmax = 0D0
   nmodes = 0

! Loop over number of possible modes
   LOOP_ALLEQ: &
   Do j=1,neq

      ! Start an error check on the eigensolution comparin RHS w LHS
      Do j1=1,neq
         ztempa(j1) = 0D0
         ztempb(j1) = 0D0

         Do j2=1,neq
            zerreal =                                           &
                  zamrt(j1,j2) * DBLE (zevec(j2,j))             &
                - zamit(j1,j2) * DIMAG(zevec(j2,j))             &
                - zbmrt(j1,j2) * DBLE (zevec(j2,j)) * zomega(j) &
                + zbmrt(j1,j2) * DIMAG(zevec(j2,j)) * zgamma(j) &
                + zbmit(j1,j2) * DBLE (zevec(j2,j)) * zgamma(j) &
                + zbmit(j1,j2) * DIMAG(zevec(j2,j)) * zomega(j)

            zerimag =                                           &
                  zamrt(j1,j2) * DIMAG(zevec(j2,j))             &
                + zamit(j1,j2) * DBLE (zevec(j2,j))             &
                - zbmrt(j1,j2) * DIMAG(zevec(j2,j)) * zomega(j) &
                - zbmrt(j1,j2) * DBLE (zevec(j2,j)) * zgamma(j) &
                + zbmit(j1,j2) * DIMAG(zevec(j2,j)) * zgamma(j) &
                - zbmit(j1,j2) * DBLE (zevec(j2,j)) * zomega(j)

            ztempa(j1) = ztempa(j1) + zerreal
            ztempb(j1) = ztempb(j1) + zerimag
         End Do

         zerrmax = max ( zerrmax, abs(ztempa(j1)), abs(ztempb(j1)) )
      End Do

      !.. Compute effective diffusivities directly from eigenvalues
      ! Assuming eigenvectors are arranged in the order of
      !    e\Phi/Te, Ti, nH, Te, A//, v//

      !.. Calculate norm of solution in terms of (e\phi/T_e) **2
      ztemp1 = DBLE(zevec(1,j)) *  DBLE(zevec(1,j)) &
            + DIMAG(zevec(1,j)) * DIMAG(zevec(1,j))

      !.. Check if current mode j is unstable
      If ( zgamma(j)>zepsqrt .and. ztemp1>zepsqrt ) Then

         nmodes = nmodes + 1

         !.. Define fluxes : Thermal hydrogen flux
         zreal =  DBLE(zevec(2,j)) +  DBLE(zevec(3,j))
         zimag = DIMAG(zevec(2,j)) + DIMAG(zevec(3,j))
         !
         !.. phase difference
         zphsph = - ( zimag * DBLE(zevec(1,j) ) &
                  -   zreal *DIMAG(zevec(1,j) ) ) / ztemp1
         !
         !.. flux from j:th mode
         zflxph = 2D0* zphsph * zgamma(j) * zgamma(j)
         !
         !.. flux summed over all unstable modes
         zflxm(1) = zflxm(1) + zflxph
         !
         !.. Define hydrogen density flux - phase shift
         zphsnh = - ( DIMAG(zevec(3,j)) * DBLE(zevec(1,j)) &
                - DBLE(zevec(3,j)) * DIMAG(zevec(1,j)) ) / ztemp1
         !
         !.. Flux from j:th mode
         zflxnh = 2D0* zphsnh * zgamma(j) * zgamma(j)
         !
         !.. Flux summed over all unstable modes
         zflxm(2) = zflxm(2) + zflxnh
         !
         !.. Define thermal electron flux
         zreal =  DBLE(zevec(4,j)) +  DBLE(zevec(3,j))
         zimag =  DIMAG(zevec(4,j)) +  DIMAG(zevec(3,j))
         !
         zphspe = - ( zimag * DBLE(zevec(1,j)) &
                  -  zreal * DIMAG(zevec(1,j)) ) / ztemp1
         !
         zflxpe = 2D0* zphspe * zgamma(j)* zgamma(j)
         zflxm(3) = zflxm(3) + zflxpe

      End If

   End Do LOOP_ALLEQ !.. End of flux and diffusivity definitions loop

   !..compute effective total diffusivities
   zchim(1) = zflxm(1) / sign( max( abs( gthin ), zepsqrt ),  gthin )
   zchim(2) = zflxm(2) / sign( max( abs( gnhin ), zepsqrt ),  gnhin )
   zchim(3) = zflxm(3) / sign( max( abs( gtein ), zepsqrt ),  gtein )

   !..save effective diffusivities and fluxes
   chieff(1:3) = zchim(1:3)

End If IF_HAVETXP! No transport }}}

End Subroutine drbm!}}}

End Module modmmm7_1
