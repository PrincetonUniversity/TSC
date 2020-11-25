      subroutine mmm95 (                                                 &  
     &   rminor,  rmajor,   elong                                        &  
     & , dense,   densh,    densimp,  densfe                             &  
     & , xzeff,   tekev,    tikev,    q,       btor                      &  
     & , avezimp, amassimp, amasshyd, aimass,  wexbs                     &  
     & , grdne,   grdni,    grdnh,    grdnz,   grdte,   grdti,  grdq     &  
     & , thiig,   thdig,    theig,    thzig                              &  
     & , thirb,   thdrb,    therb,    thzrb                              &  
     & , thikb,   thdkb,    thekb,    thzkb                              &  
     & , gamma,   omega,    difthi,   velthi,  vflux                     &  
     & , matdim,  npoints,  nprout,   lprint,  nerr                      &  
     & , lsuper,  lreset,   lswitch,  cswitch, fig,    frb,     fkb)
!
 
 
 
! 29mar1999 fgtok -s cgg.table "dmc:  rename eispack routines"
! 29mar1999 fgtok -s rr.table "dmc:  rename intrinsic REAL -> REAL"
!
! dmc -- Cray/workstation portable real*8<-->complex*16 conversion routines
 
!      include "f77_dcomplx.h"
 
!
!glf2d.f 11-nov-99 Kinsey
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
 
!
!    All the following 1-D arrays are assumed to be defined on flux
!    surfaces called zone boundaries where the transport fluxes are
!    to be computed.  The number of flux surfaces is given by npoints
!    (see below).  For example, if you want to compute the transport
!    on only one flux surface, set npoints = 1.
!
!    Note, if the mmm95 module is used in another code, the minimum
!    dimension required for cswitch is 23, the minimum dimension for
!    lswitch is 5 and fig, frb, and fkb must all be dimensioned to a
!    minimum of 4.
 
!  Input arrays:
!  -------------
!
!  rminor(jz)   = minor radius (half-width) of zone boundary [m]
!  rmajor(jz)   = major radius to geometric center of zone bndry [m]
!  elong(jz)    = local elongation of zone boundary
!
!  dense(jz)    = electron density [m^-3]
!  densh(jz)    = sum over thermal hydrogenic ion densities [m^-3]
!  densimp(jz)  = sum over impurity ion densities [m^-3]
!  densfe(jz)   = electron density from fast (non-thermal) ions [m^-3]
!
!  xzeff(jz)    = Z_eff
!  tekev(jz)    = T_e (electron temperature) [keV]
!  tikev(jz)    = T_i (temperature of thermal ions) [keV]
!  q(jz)        = magnetic q-value
!  btor(jz)     = ( R B_tor ) / rmajor(jz)  [tesla]
!
!  avezimp(jz)  = average density weighted charge of impurities
!               = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
!                 sum_imp = sum over impurity ions with charge state Z_imp
!
!  amassimp(jz) = average density weighted atomic mass of impurities
!               = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where
!                 sum_imp = sum over impurity ions, each with mass M_imp
!
!  amasshyd(jz) = average density weighted atomic mass of hydrogen ions
!               = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
!                 sum_hyd = sum over hydrogenic ions, each with mass M_hyd
!
!  aimass(jz)   = mean atomic mass of thermal ions [AMU]
!               = ( sum_i n_i M_i ) / ( sum_i n_i ) where
!                 sum_i = sum over all ions, each with mass M_i
!
!  wexbs(jz)    = ExB shearing rate in [rad/s].  See  K.H. Burrell,
!                 "Effects of {ExB} velocity shear and magnetic shear
!                 on turbulence and transport in magnetic confinement
!                 devices", Phys. of Plasmas, 4, 1499 (1997).
!
!    All of the following normalized gradients are at zone boundaries.
!    r = half-width, R = major radius to center of flux surface
!
!  grdne(jz) = -R ( d n_e / d r ) / n_e
!  grdni(jz) = -R ( d n_i / d r ) / n_i
!  grdnh(jz) = -R ( d n_h / d r ) / n_h
!  grdnz(jz) = -R ( d Z n_Z / d r ) / ( Z n_Z )
!  grdte(jz) = -R ( d T_e / d r ) / T_e
!  grdti(jz) = -R ( d T_i / d r ) / T_i
!  grdq (jz) =  R ( d q   / d r ) / q    related to magnetic shear
!
!  where:
!    n_i     = thermal ion density (sum over hydrogenic and impurity)
!    n_h     = thermal hydrogenic density (sum over hydrogenic species)
!    n_Z     = thermal impurity density,  Z = average impurity charge
!                      sumed over all impurities
!
!  Output:
!  -------
!
!    The following effective diffusivities represent contributions
!    to the total diffusivity matrix (difthi and velthi given below)
!    from each of the models that contribute to the Multi-Mode model.
!    Generally, these arrays are used for diagnostic output only.
!
!  thiig(jz) = ion thermal diffusivity from the Weiland model
!  thdig(jz) = hydrogenic ion diffusivity from the Weiland model
!  theig(jz) = elelctron thermal diffusivity from the Weiland model
!  thzig(jz) = impurity ion diffusivity from the Weiland model
!
!  thirb(jz) = ion thermal diffusivity from resistive ballooning modes
!  thdrb(jz) = hydrogenic ion diffusivity from resistive ballooning modes
!  therb(jz) = elelctron thermal diffusivity from resistive ballooning modes
!  thzrb(jz) = impurity ion diffusivity from resistive ballooning modes
!
!  thikb(jz) = ion thermal diffusivity from kinetic ballooning modes
!  thdkb(jz) = hydrogenic ion diffusivity from kinetic ballooning modes
!  thekb(jz) = elelctron thermal diffusivity from kinetic ballooning modes
!  thzkb(jz) = impurity ion diffusivity from kinetic ballooning modes
!
!    The following are growth rates and mode frequencies from the
!    Weiland model for drift modes such as ITG and TEM.
!    These arrays are intended for diagnostic output.
!
!  gamma(jm,jz) = growth rate for mode jm at point jz ( 1/sec )
!  omega(jm,jz) = frequency for mode jm at point jz ( radians/sec )
!
!    All of the transport coefficients are given in the following two
!    matricies for diffusion difthi and convection velthi in MKS units.
!    See the LaTeX documentation for difthi and velthi just below.
!
!    NOTE:  difthi and velthi include all of the anomalous transport.
!    There are no additional contributions to the heat fluxs from
!    charged particle convection.
!
!  difthi(j1,j2,jz) = full matrix of anomalous transport diffusivities
!  velthi(j1,jz)    = convective velocities
!  vflux(j1,jz)     = flux matrix
!
!  Input integers:
!  ---------------
!
!  matdim  = first and second dimension of transport matricies
!            difthi(j1,j2,jz) and the first dimension of
!            velthi(j1,jz), vflux(j1,jz), gamma(j1,jz), and omega(j1,jz).
!            matdim must be at least 5
!
!  npoints = number of values of jz in all of the above arrays
!
!  nprout  = output unit number for long printout
!
!
!  Input switches
!  --------------
!
!  lprint      controls the amount of printout (0 => no printout)
!              higher values yield more diagnostic output
!
!  lsuper   = 0 for simulations of all other discharges
!           > 0 for supershot simulations; substantially reduces
!               contribution from kinetic ballooning mode
!
!
!
!
!  lreset  = 0 to use only internal settings for lswitch, cswitch
!              and for the coefficients fig, frb, and fkb that control
!              the contributions form the various instability modes
!
!    Note that when lreset = 0, the values of the switches and
!    coefficients in the argument list are ignored and all the
!    switches and coefficients are set internally.
!
!    WARNING:  use lreset > 0 only if you want to pass all the switches
!              lswitch, cswitch, fig, frb, and fkb through the
!              argument list.
!
!    WARNING:  NTCC users desiring to use anything other than lreset = 0
!              should consult with the mmm95 code authors first.
!
!
!  Output Integer
!  --------------
!
!  nerr        status code returned; 0 = OK; .ne. 0 indicates error
!
!
!  Internal control variables:
!  ---------------------------
!
!  lswitch(j), j=1,8   integer control variables:
!
!  cswitch(j), j=1,25   general control variables:
!
!  lswitch(1)  controls which version of the Weiland model is used
!                  Default lswitch(1) = 10
!             = 2  2 eqn  Weiland model Hydrogen \eta_i mode only
!             = 4  4 eqn  Weiland model with Hydrogen and trapped electrons
!             = 5  5 eqn  Weiland model with trapped electrons, FLR effects,
!                         and parallel ion motion
!             = 6  6 eqn  Weiland model Hydrogen, trapped electrons,
!                    and one impurity species
!             = 7  7 eqn   Weiland model Hydrogen, trapped electrons,
!                  one impurity species, and collisions
!             = 8  8 eqn  Weiland model Hydrogen, trapped electrons,
!                  one impurity species, collisions, and parallel
!                  ion (hydrogenic) motion
!             = 9  9 eqn  Weiland model Hydrogen, trapped electrons,
!                  one impurity species, collisions, and finite beta
!             = 10 10 eqn Weiland model Hydrogen, trapped electrons,
!                  one impurity species, collisions, parallel
!                  ion (hydrogenic) motion, and finite beta
!             = 11 11 eqn Weiland model Hydrogen, trapped electrons,
!                  one impurity species, collisions, parallel
!                  ion (hydrogenic, impurity) motion, and finite beta
!
!  lswitch(2) = 0  full matrix representation for difthi and velthi
!                  Default lswitch(2) = 2
!             = 1  set diagonal matrix elements difthi and velthi
!             = 2  set diagonal matrix elements = effective diffusivities
!
!  lswitch(3)  controls \kappa scaling
!                  Default lswitch(3) = 0
!             = 0  use \kappa scaling raised to
!                  exponents (cswitch(3) - cswitch(5))
!             = 1  use (1+\kappa^2)/2 instead of \kappa scaling
!
!  lswitch(4) > 0  to replace negative diffusivity with velocity
!                  Default lswitch(4) = 1
!
!  lswitch(5) = 1  to limit magnitude of all normalized gradients
!                  to ( major radius ) / ( ion Larmor radius )
!                  Default lswitch(5) = 1
!
!  cswitch(1)   0.5  minimum value of shear
!  cswitch(2)   3.5  coeff in exponential (fbeta-th) of kinetic ballooning model
!  cswitch(3)  -4.0  exponent of local elongation multiplying drift waves
!  cswitch(4)  -4.0  exponent of local elongation multiplying resistive
!                     ballooning modes
!  cswitch(5)  -4.0  exponent of local elongation multiplying
!                     kinetic balllooning modes
!  cswitch(6)   0.0  k_y \rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
!  cswitch(8)   1.0  coeff of beta_prime_1 in kinetic ballooning mode
!  cswitch(9)  0.15  alpha in diamagnetic stabil. in kinetic ballooning model
!  cswitch(10)  0.0  rel fract of ion thermal diffusivity given to convection
!  cswitch(11)  0.0  rel fract of hydrogen particle diffusivity given to convect
!  cswitch(12)  0.0  rel fract of el thermal diffusivity given to convection
!  cswitch(13)  0.0  rel fract of impurity particle diffusivity given to convect
!  cswitch(14)  1.0  coef of finite beta effect in weiland14 = cetain(20)
!  cswitch(15)  0.0  min value of impurity charge state zimpz
!  cswitch(16)  0.0  coef of fast particle fraction (superthermal ions)
!                    in weiland model -- coef of densfe
!  cswitch(17)  1.0  coeff of k_\parallel (parallel ion motion) in weiland14
!                    = cetain(10)
!  cswitch(18)  0.0  coeff of nuhat (effect of collisions) in weiland14
!                    = cetain(15)
!  cswitch(19)  0.0  coeff for including v_parallel in strong ballooning limit
!                    = cetain(12); cswitch(19) = 1 for inclusion of v_par effect
!  cswitch(20)  0.0  trapping fraction used in weiland14 (when > 0.0)
!                    multiplies electron trapping fraction when < 0.0
!                    no effect when cswitch(20) = 0.0
!  cswitch(21)  1.0  multiplier for wexbs (flow shear rate) in Weiland model
!  cswitch(22)  0.0  ranges from 0.0 to 1.0 adds impurity heat flow to total
!                    ionic heat flow for the weiland model
!  cswitch(23)  0.0  controls finite diff to construct the zgm matrix
!                    = cetain(30)
!
!     contributions to vfluxes and interchanges:
!
!  fig(1)   hydrogen particle transport from ITG (eta_i) mode
!  fig(2)   electron thermal  transport from ITG (eta_i) mode
!  fig(3)   ion      thermal  transport from ITG (eta_i) mode
!  fig(4)   impurity particle transport from ITG (eta_i) mode
!
!  frb(1)   hydrogen particle transport from resistive ballooning mode
!  frb(2)   electron thermal  transport from resistive ballooning mode
!  frb(3)   ion      thermal  transport from resistive ballooning mode
!  frb(4)   impurity particle transport from resistive ballooning mode
!
!  fkb(1)   hydrogen particle transport from kinetic ballooning mode
!  fkb(2)   electron thermal  transport from kinetic ballooning mode
!  fkb(3)   ion      thermal  transport from kinetic ballooning mode
!  fkb(4)   impurity particle transport from kinetic ballooning mode
!
!
!***********************************************************************
!
!-----------------------------------------------------------------------
!
!  Compile this routine and routines that it calls with a compiler
!  option, such as -r8, to convert real to double precision when used on
!  workstations.
!
!-----------------------------------------------------------------------
!
!  External dependencies:
!
!  Call tree: MMM95 calls the following routines
!
!  WEILAND14       - Computes diffusion matrix and convect velocities
!                        for the Weiland transport model
!    WEILAND14FLUX - Calculates fluxes and effective diffusivities
!      TOMSQZ      - Wrapper for QZ algorithm solving Ax = lambda Bx
!         CQZHES   - First step in QZ algorithm
!         CQZVAL   - Second and third step in QZ algorithm
!         CQZVEC   - Fourth step in QZ algorithm
!
!-----------------------------------------------------------------------
 
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer km, klswitch, kcswitch
!
      integer  matdim,  npoints, nprout,  lprint,   nerr
!
      integer  lsuper,  lreset,  lswitch(*)
!
      parameter ( km = 12, klswitch = 8, kcswitch = 25 )
!
      REAL*8                                                             &  
     &   rminor(*),  rmajor(*),   elong(*)                               &  
     & , dense(*),   densh(*),    densimp(*),  densfe(*)                 &  
     & , xzeff(*),   tekev(*),    tikev(*),    q(*),       btor(*)       &  
     & , avezimp(*), amassimp(*), amasshyd(*), aimass(*),  wexbs(*)      &  
     & , grdne(*),   grdni(*),    grdnh(*),    grdnz(*)                  &  
     & , grdte(*),   grdti(*),    grdq(*)
!
      REAL*8                                                             &  
     &   thiig(*),   thdig(*),    theig(*),    thzig(*)                  &  
     & , thirb(*),   thdrb(*),    therb(*),    thzrb(*)                  &  
     & , thikb(*),   thdkb(*),    thekb(*),    thzkb(*)                  &  
     & , omega(matdim,*),         gamma(matdim,*)                        &  
     & , difthi(matdim,matdim,*), velthi(matdim,*)                       &  
     & , vflux(matdim,*)
!
      REAL*8     cswitch(*)
!
      REAL*8     fig(*),  fkb(*),  frb(*)
!
!..physical constants
!
      REAL*8 zpi,  zcc,  zcmu0,  zceps0,  zckb,  zcme,  zcmp,  zce
!
!  zpi     = pi
!  zcc     = speed of light                  [m/sec]
!  zcmu0   = vacuum magnetic permeability    [henrys/m]
!  zceps0  = vacuum electrical permittivity  [farads/m]
!  zckb    = energy conversion factor        [Joule/keV]
!  zcme    = electron mass                   [kg]
!  zcmp    = proton mass                     [kg]
!  zce     = electron charge                 [Coulomb]
!
!..computer constants
!
      REAL*8  zepslon, zlgeps
!
!  zepslon = machine epsilon [smallest number so that 1.0+zepslon>1.0]
!  zlgeps  = ln ( zepslon )
!
!
!..local variables
!
      integer  jz, j1, j2, jm
 
      REAL*8  zelong,  zelonf,  zai,     zne,     zni,    zte,    zti    &  
     & ,    zq,      zeff,    zgne,    zgni,    zgnh,   zgnz,   zgte     &  
     & ,    zgth,    zgtz,    zshear,  zrmin,   zrmaj,  zbtor,  zep      &  
     & ,    zgyrfi,  zbeta,   zvthe,   zvthi,   zsound, zlog,   zcf      &  
     & ,    znuei,   znueff,  zlari,   zlarpo,  zrhos,  zwn,    znude    &  
     & ,    znuhat,  zgpr,    zscyl,   zsmin,   zshat,  zgmax
!
!.. variables for Weiland model
!
!  iletai(j1) and cetain(j1) are control variables
!
      integer        iletai(32)
!
      REAL*8  cetain(32), zomega(km), zgamma(km), zchieff(km)
!
      REAL*8           zdfthi(km,km),    zvlthi(km),     zflux(km)
!
      integer        idim,    ieq,     imodes
!
      REAL*8  zthte,   zbetae,  znz,     zmass,  zimpz,  ztzte           &  
     & ,    zfnzne,  zmzmh,   zfnsne,  zftrap, zkyrho, zomegde           &  
     & ,    zwexb,   znormd,  znormv
!
!  zexb    = local copy of ExB shearing rate
!  znormd  = factor to convert normalized diffusivities
!  znormv  = factor to convert normalized convective velocities
!
!..local variables for resistive ballooning modes
!
      REAL*8  zgyrfe, zlare,   zgddia, zgdp
!
!..local variables for kinetic ballooning modes
!
      REAL*8  zbprim, zbcoef1, zbc1,   zelfkb,  zfbthn, zdk
!
!-----------------------------------------------------------------------
!
!..initialize imodes
      imodes  = 0
!..physical constants
!
        zpi     = atan2 ( 0.0_R8, -1.0_R8)
        zcc     = 2.997925E+8_R8
        zcmu0   = 4.0E-7_R8* zpi
        zceps0  = 1.0_R8/ ( zcc**2 * zcmu0 )
        zckb    = 1.60210E-16_R8
        zcme    = 9.1091E-31_R8
        zcmp    = 1.67252E-27_R8
        zce     = 1.60210E-19_R8
!
!..computer constants
!
        zepslon = 1.0E-34_R8
        zlgeps  = log ( zepslon )
!
!
!..initialize arrays
!
      do jz = 1, npoints
        thiig(jz)  = 0._R8
        thdig(jz)  = 0._R8
        theig(jz)  = 0._R8
        thzig(jz)  = 0._R8
        therb(jz)  = 0._R8
        thirb(jz)  = 0._R8
        thdkb(jz)  = 0._R8
        thekb(jz)  = 0._R8
        thzkb(jz)  = 0._R8
        thikb(jz)  = 0._R8
        thdkb(jz)  = 0._R8
        thekb(jz)  = 0._R8
        thzkb(jz)  = 0._R8
      enddo
!
      do jz = 1, npoints
        do j1 = 1, matdim
          velthi(j1,jz) = 0.0_R8
          vflux(j1,jz) = 0.0_R8
          gamma(j1,jz) = 0.0_R8
          omega(j1,jz) = 0.0_R8
          do j2 = 1, matdim
            difthi(j1,j2,jz) = 0.0_R8
          enddo
        enddo
      enddo
!
      nerr = 0
!
!..if lreset < 1, use internal settings for switches and coefficients
!  otherwise, use values passed through the argument list above
!
      if ( lreset .lt. 1 ) then
!
!..initialize switches
!
      do j1=1,kcswitch
        cswitch(j1) = 0.0_R8
      enddo
!
      do j1=1,klswitch
        lswitch(j1) = 0
      enddo
!
!
!  Multi Mode Model in sbrtn THEORY version MMM95
!  for use in the BALDUR transport code
!
      lswitch(1) = 10 ! Weiland ITG model weiland14 (10 eqns, no collisions)
      lswitch(2) = 2  ! use effective diffusivities
      lswitch(3) = 0  ! use kappa instead of (1+\kappa^2)/2
      lswitch(4) = 1  ! replace -ve diffusivity with convective velocity
      lswitch(5) = 1  ! limit gradients by major radius / ion Larmor radius
!
!  misc. parameters for subroutine mmm95
!
      cswitch(1)  =  0.5_R8  ! minimum value of shear
      cswitch(2)  =  3.5_R8  ! coef in exponential (fbeta-th) in kinetic ballooning
      cswitch(3)  = -4.0_R8  ! elongation scaling for drift wave mode
      cswitch(4)  = -4.0_R8  ! elongation scaling for resistive ballooning mode
      cswitch(5)  = -4.0_R8  ! elongation scaling for kinetic ballooning mode
      cswitch(6)  =  0.0_R8  ! k_y \rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
      cswitch(8)  =  1.0_R8  ! coeff of beta_prime_1 in kinetic ballooning mode
      cswitch(9)  = 0.15_R8  ! alpha in diamagnetic stabil. in kinetic balloon mode
      cswitch(10) =  0.0_R8  ! relative fraction of ion thermal diffusivity
                          ! given to convection
      cswitch(11) =  0.0_R8  ! relative fract of hydrogen particle diffusivity
                          ! given to convection
      cswitch(12) =  0.0_R8  ! relative fraction of el thermal diffusivity
                          ! given to convection
      cswitch(13) =  0.0_R8  ! relative fract of impurity particle diffusivity
                          ! given to convection
      cswitch(14) =  1.0_R8  ! coef of finite beta effect in weiland14 = cetain(20)
      cswitch(15) =  0.0_R8  ! min value of impurity charge state zimpz
      cswitch(16) =  1.0_R8  ! coef of fast particle fraction (superthermal ions)
                          ! in weiland model -- coef of densfe
      cswitch(17) =  1.0_R8  ! coeff of k_\parallel (parallel ion motion) in
                          ! weiland14 = cetain(10)
      cswitch(18) =  0.0_R8  ! coeff of nuhat in weiland14 = cetain(15)
      cswitch(19) =  0.0_R8  ! coeff for including v_parallel in strong ballooning
                          ! limit = cetain(12); cswitch(19) = 1 for inclusion
                          ! of v_par effect
      cswitch(20) =  0.0_R8  ! trapping fraction used in weiland14 (when > 0.0)
                          ! multiplies electron trapping fraction when < 0.0
                          ! no effect when cswitch(20) = 0.0
      cswitch(21) =  1.0_R8  ! multiplier for wexbs (flow shear rate)
                          ! in Weiland model
      cswitch(22) =  0.0_R8  ! multiplier to impurity heat flux
      cswitch(23) =  0.0_R8  ! controls finite diff to construct the
                          ! zgm matrix = cetain(30)
 
!  contributions to hydrogenic particle, elec-energy, ion-energy,
!    and impurity ion fluxes
 
        fig(1) = 0.80_R8
        fig(2) = 0.80_R8
        fig(3) = 0.80_R8
        fig(4) = 0.80_R8
!
        fkb(1) = 1.00_R8
         fkb(2) = 0.65_R8
        fkb(3) = 0.65_R8
        fkb(4) = 1.00_R8
!
        if ( lsuper .gt. 0 ) then
          fkb(1) = 0.045_R8
          fkb(2) = 0.010_R8
          fkb(3) = 0.010_R8
          fkb(4) = 0.045_R8
        endif
!
        frb(1) = 1.00_R8
        frb(2) = 1.00_R8
        frb(3) = 1.00_R8
        frb(4) = 1.00_R8
!
      endif
 
!
!.. start the main do-loop over the radial index "jz"..........
!
!
      do 300 jz = 1, npoints
!
!  transfer common to local variables to compact the notation
!
      zelong = max (zepslon,elong(jz))
      if ( lswitch(3) .eq. 1 ) then
        zelonf = ( 1._R8+ zelong**2 ) / 2._R8
      else
        zelonf = zelong
      endif
!
      zai    = aimass(jz)
      zne    = dense(jz)
      zni    = densh(jz) + densimp(jz)
      znz    = densimp(jz)
      zte    = tekev(jz)
      zti    = tikev(jz)
      zq     = q(jz)
      zeff   = xzeff(jz)
!
!  normalized gradients
!
      zgne   = grdne(jz)
      zgni   = grdni(jz)
      zgnh   = grdnh(jz)
      zgnz   = grdnz(jz)
      zgte   = grdte(jz)
      zgth   = grdti(jz)
      zgtz   = grdti(jz)
 
      zrmin  = max( rminor(jz), zepslon )
      zrmaj  = rmajor(jz)
      zshear = grdq(jz) * zrmin / zrmaj
      zbtor  = btor(jz)
!
!  compute inverse aspect ratio
!
      zep    = max( zrmin/zrmaj, zepslon )
!
!
      zgyrfi = zce * zbtor / (zcmp * zai)
      zbeta  = (2._R8* zcmu0 * zckb / zbtor**2) * (zne * zte + zni *     &  
     & zti)
      zvthe  = sqrt(2._R8* zckb * zte / zcme)
      zvthi  = sqrt(2._R8* zckb * zti / (zcmp * zai))
      zsound = sqrt(zckb * zte / (zcmp * zai))
      zlog   = 37.8_R8-log(sqrt(zne) / zte)
      zcf    = (4._R8* sqrt(zpi) / 3._R8)
      zcf    = zcf * (zce / (4._R8* zpi * zceps0))**2
      zcf    = zcf * (zce / zckb) * sqrt( (zce/zcme) * (zce/zckb) )
      znuei  = zcf * sqrt(2._R8) * zne * zlog * zeff / (zte * sqrt(zte))    
!
      znueff = znuei / zep
      zlari  = zvthi / zgyrfi
      zlarpo = max(zlari * zq / zep, zepslon)
      zrhos  = zsound / zgyrfi
      zwn    = 0.3_R8/ zrhos
      znude  = 2 * zwn * zrhos * zsound / zrmaj
      znuhat = znueff / znude
!
!..if lswitch(5) = 1, limit magnitude of normalized gradients
!                    to ( major radius ) / ( ion Larmor radius )
!
      zgmax = zrmaj / zlarpo
!
      if ( lswitch(5) .eq. 1 ) then
!
        zgne = sign ( min ( abs ( zgne ), zgmax ), zgne )
        zgni = sign ( min ( abs ( zgni ), zgmax ), zgni )
        zgnh = sign ( min ( abs ( zgnh ), zgmax ), zgnh )
        zgnz = sign ( min ( abs ( zgnz ), zgmax ), zgnz )
        zgte = sign ( min ( abs ( zgte ), zgmax ), zgte )
        zgth = sign ( min ( abs ( zgth ), zgmax ), zgth )
        zgtz = sign ( min ( abs ( zgtz ), zgmax ), zgtz )
!
      endif
!
!  zgpr = -R ( d p   / d r ) / p    for thermal pressure
!
!  Compute the pressure scale length using smoothed and bounded
!  density and temperature
!
      zgpr = ( zne * zte * ( zgne + zgte )                               &  
     &         + zni * zti * ( zgni + zgth ) )                           &  
     &         / ( zne * zte + zni * zti )
!
      if ( lswitch(5) .eq. 1 )                                           &  
     &  zgpr = sign ( min ( abs ( zgpr ), zgmax ), zgpr )
!
!
!
      zscyl=max(abs(zshear),zepslon)
      zsmin=max(cswitch(1),zepslon)
      zshat=max(zsmin,zscyl)
!
!
        do j1=1,32
          iletai(j1) = 0
          cetain(j1) = 0.0_R8
        enddo
!
        thiig(jz) = 0.0_R8
        theig(jz) = 0.0_R8
        thdig(jz) = 0.0_R8
        thzig(jz) = 0.0_R8
!
!..set the number of equations to use in the Weiland model
!
        if ( (lswitch(1) .lt. 2) .or. (lswitch(1) .gt. 11 )) then
          nerr = -10
          return
        elseif (lswitch(1) .eq. 3) then
        nerr = -10
          return
        else
          ieq = lswitch(1)
        endif
!
        cetain(11) = 1.0_R8
!
!.. coefficient of k_parallel for parallel ion motion
!.. cswitch(19) for v_parallel in strong ballooning limit
!.. in 9 eqn model
!
        cetain(10) = cswitch(17)
        cetain(12) = cswitch(19)
        cetain(15) = cswitch(18)
        cetain(20) = cswitch(14)
!
        iletai(10) = 0
!
        idim   = km
!
!  Hydrogen species
!
        zthte  = zti / zte
 
        zbetae = 2._R8* zcmu0 * zckb * zne * zte / zbtor**2
!
!  Impurity species (use only impurity species 1 for now)
!  assume T_Z = T_H throughout the plasma here
!
        znz    = densimp(jz)
        zmass  = amassimp(jz)
        zimpz  = avezimp(jz)
        zimpz  = max ( zimpz, cswitch(15) )
!
        ztzte  = zti / zte
        zfnzne = znz / zne
        zmzmh  = zmass / amasshyd(jz)
!
!  superthermal ions
!
!  zfnsne = ratio of superthermal ions to electrons
!  L_ns   = gradient length of superthermal ions
!
        zfnsne = max ( cswitch(16) * densfe(jz) / dense(jz) , 0.0_R8)
!
        zftrap = sqrt ( 2._R8* zrmin / ( zrmaj * ( 1._R8+ zrmin / zrmaj   &  
     & )))
        if ( cswitch(20) .gt. zepslon ) zftrap = cswitch(20)
        if ( cswitch(20) .lt. -zepslon )                                 &  
     &       zftrap = abs(cswitch(20))*zftrap
!
        if ( abs(cswitch(6)) .lt. zepslon ) then
          zkyrho = 0.316_R8
        else
          zkyrho = cswitch(6)
        endif
!
!
!...Define a local copy of normalized ExB shearing rate : pis
!
        zomegde = 2.0_R8* zkyrho * zsound / zrmaj
!
        zwexb = cswitch(21) * wexbs(jz) / zomegde
!
!
        cetain(30) = cswitch(23)
        iletai(6)  = 0
        if ( lswitch(2) .lt. 1 ) iletai(7) = 1
!
!  if lswitch(2) .lt. 1, compute only the effective diffusivities
!
        iletai(9) = lswitch(2)
!
        call weiland14 (                                                 &  
     &   iletai,   cetain,   lprint,   ieq,      nprout,   zgne          &  
     & , zgnh,     zgnz,     zgte,     zgth,     zgtz,     zthte         &  
     & , ztzte,    zfnzne,   zimpz,    zmzmh,    zfnsne,   zbetae        &  
     & , zftrap,   znuhat,   zq,       zshat,    zkyrho,   zwexb         &  
     & , idim,     zomega,   zgamma,   zdfthi,   zvlthi,   zchieff       &  
     & , zflux,    imodes,   nerr )
!
!  If nerr not equal to 0 an error has occured
!
      if (nerr .ne. 0) return
!
!
!  Growth rates for diagnostic output
!    Note that all frequencies are normalized by \omega_{De}
!      consequently, trapped electron modes rotate in the positive
!      direction (zomega > 0) while eta_i modes have zomega < 0.
!
        jm = 0
        do j1=1,imodes
          if ( zgamma(j1) .gt. zepslon ) then
            jm = jm + 1
            gamma(jm,jz) = zgamma(j1) / zomegde
            omega(jm,jz) = zomega(j1) / zomegde
          endif
        enddo
!
!  conversion factors for diffusion and convective velocity
!
        znormd = zelonf**cswitch(3) *                                    &  
     &    2.0_R8* zsound * zrhos**2 / ( zrmaj * zkyrho )
        znormv = zelonf**cswitch(3) *                                    &  
     &    2.0_R8* zsound * zrhos**2 / ( zrmaj**2 * zkyrho )
!
!  compute effective diffusivites for diagnostic purposes only
!
        thdig(jz) = fig(1) * znormd * zchieff(2)
        theig(jz) = fig(2) * znormd * zchieff(3)
        thiig(jz) = fig(3) * znormd * zchieff(1)                         &  
     &  + fig(3) * znormd * zchieff(5) * cswitch(22) * znz / zni
        thzig(jz) = fig(4) * znormd * zchieff(4)
!
!  start computing the fluxes
!
        vflux(1,jz) = vflux(1,jz) + thiig(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdig(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + theig(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzig(jz) * zgnz / zrmaj
!
!  compute diffusivity matrix
!
        do j1=1,matdim
          velthi(j1,jz) = 0.0_R8
          vflux(j1,jz) = 0.0_R8
          do j2=1,matdim
            difthi(j1,j2,jz) = 0.0_R8
          enddo
        enddo
!
!..set difthi and velthi
!
        if ( lswitch(2) .gt. 1 ) then
!
!  diagonal elements of matrix = effective diffusivities
!
          difthi(1,1,jz) = difthi(1,1,jz) + thiig(jz)
          difthi(2,2,jz) = difthi(2,2,jz) + thdig(jz)
          difthi(3,3,jz) = difthi(3,3,jz) + theig(jz)
          difthi(4,4,jz) = difthi(4,4,jz) + thzig(jz)
!
        else
!
!..full matrix form of model
!
          if ( ieq .eq. 2 ) then
            difthi(1,1,jz) = difthi(1,1,jz) +                            &  
     &        fig(3) * znormd * zdfthi(1,1)
            velthi(1,jz)   = velthi(1,jz) +                              &  
     &        fig(3) * znormv * zvlthi(1)
          elseif ( ieq .eq. 4 ) then
            do j2=1,3
              difthi(1,j2,jz) = difthi(1,j2,jz) +                        &  
     &          fig(3) * znormd * zdfthi(1,j2)
              difthi(2,j2,jz) = difthi(2,j2,jz) +                        &  
     &          fig(1) * znormd * zdfthi(2,j2)
              difthi(3,j2,jz) = difthi(3,j2,jz) +                        &  
     &          fig(2) * znormd * zdfthi(3,j2)
              difthi(4,j2,jz) = difthi(4,j2,jz) +                        &  
     &          fig(4) * znormd * zdfthi(2,j2)
            enddo
              velthi(1,jz)    = velthi(1,jz) +                           &  
     &          fig(3) * znormv * zvlthi(1)
              velthi(2,jz)    = velthi(2,jz) +                           &  
     &          fig(1) * znormv * zvlthi(2)
              velthi(3,jz)    = velthi(3,jz) +                           &  
     &          fig(2) * znormv * zvlthi(3)
              velthi(4,jz)    = velthi(4,jz) +                           &  
     &          fig(4) * znormv * zvlthi(2)
          else
            do j2=1,4
              difthi(1,j2,jz) = difthi(1,j2,jz) +                        &  
     &          fig(3) * znormd * zdfthi(1,j2)
              difthi(2,j2,jz) = difthi(2,j2,jz) +                        &  
     &          fig(1) * znormd * zdfthi(2,j2)
              difthi(3,j2,jz) = difthi(3,j2,jz) +                        &  
     &          fig(2) * znormd * zdfthi(3,j2)
              difthi(4,j2,jz) = difthi(4,j2,jz) +                        &  
     &          fig(4) * znormd * zdfthi(4,j2)
            enddo
              velthi(1,jz)    = velthi(1,jz) +                           &  
     &          fig(3) * znormv * zvlthi(1)
              velthi(2,jz)    = velthi(2,jz) +                           &  
     &          fig(1) * znormv * zvlthi(2)
              velthi(3,jz)    = velthi(3,jz) +                           &  
     &          fig(2) * znormv * zvlthi(3)
              velthi(4,jz)    = velthi(4,jz) +                           &  
     &          fig(4) * znormv * zvlthi(4)
          endif
!
        endif
!
!
!..transfer from diffusivity to convective velocity
!
        if ( lswitch(4) .gt. 0 ) then
!
          if ( thiig(jz) .lt. 0.0_R8) then
            velthi(1,jz) = velthi(1,jz) - thiig(jz) * zgth / zrmaj
            thiig(jz) = 0.0_R8
            do j2=1,4
              difthi(1,j2,jz) = 0.0_R8
            enddo
          endif
!
          if ( thdig(jz) .lt. 0.0_R8) then
            velthi(2,jz) = velthi(2,jz) - thdig(jz) * zgnh / zrmaj
            thdig(jz) = 0.0_R8
            do j2=1,4
              difthi(2,j2,jz) = 0.0_R8
            enddo
          endif
!
          if ( theig(jz) .lt. 0.0_R8) then
            velthi(3,jz) = velthi(3,jz) - theig(jz) * zgte / zrmaj
            theig(jz) = 0.0_R8
            do j2=1,4
              difthi(3,j2,jz) = 0.0_R8
            enddo
          endif
!
          if ( thzig(jz) .lt. 0.0_R8) then
            velthi(4,jz) = velthi(4,jz) - thzig(jz) * zgnz / zrmaj
            thzig(jz) = 0.0_R8
            do j2=1,4
              difthi(4,j2,jz) = 0.0_R8
            enddo
          endif
!
        else
!
!..shift from diffusion to convective velocity
!
          if ( abs(cswitch(10)) + abs(cswitch(11)) + abs(cswitch(12))    &  
     &       + abs(cswitch(13)) .gt. zepslon ) then
!
            velthi(1,jz) = velthi(1,jz)                                  &  
     &       + cswitch(10) * thiig(jz) * zgth / zrmaj
            velthi(2,jz) = velthi(2,jz)                                  &  
     &       + cswitch(11) * thdig(jz) * zgnh / zrmaj
            velthi(3,jz) = velthi(3,jz)                                  &  
     &       + cswitch(12) * theig(jz) * zgte / zrmaj
            velthi(4,jz) = velthi(4,jz)                                  &  
     &       + cswitch(13) * thzig(jz) * zgnz / zrmaj
!
!..alter the effective diffusivities
!  if they are used for more than diagnostic purposes
!
            thiig(jz) = ( 1.0_R8- cswitch(10) ) * thiig(jz)
            thdig(jz) = ( 1.0_R8- cswitch(11) ) * thdig(jz)
            theig(jz) = ( 1.0_R8- cswitch(12) ) * theig(jz)
            thzig(jz) = ( 1.0_R8- cswitch(13) ) * thzig(jz)
!
            do j2=1,4
              difthi(1,j2,jz) = ( 1.0_R8- cswitch(10) ) * difthi(1,j2,   &  
     & jz)
              difthi(2,j2,jz) = ( 1.0_R8- cswitch(11) ) * difthi(2,j2,   &  
     & jz)
              difthi(3,j2,jz) = ( 1.0_R8- cswitch(12) ) * difthi(3,j2,   &  
     & jz)
              difthi(4,j2,jz) = ( 1.0_R8- cswitch(13) ) * difthi(4,j2,   &  
     & jz)
            enddo
!
          endif
!
        endif
!
!..end of Weiland model
!
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
!
!..Guzdar-Drake theory (Phys Fluids B 5 (1993) 3712
!..L_p used instead of L_n
!
        zgyrfe = zce * zbtor / zcme  ! electron plasma frequency
        zlare = zvthe / zgyrfe    ! electron Larmor radius
!
!..   Diamagnetic stabilization
!
          zgddia = cswitch(9)
!
!..   Diffusivities
!
        zgdp = 2._R8* zpi * ((zq * zlare)**2._R8) * znuei                &  
     &    * zgpr * 100._R8* zgddia
 
        thdrb(jz) = frb(1) * zgdp * zelonf**cswitch(4)
        therb(jz) = frb(2) * zgdp * zelonf**cswitch(4)
        thirb(jz) = frb(3) * zgdp * zelonf**cswitch(4)
        thzrb(jz) = frb(4) * zgdp * zelonf**cswitch(4)
!
!  add to the fluxes
!
        vflux(1,jz) = vflux(1,jz) + thirb(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdrb(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + therb(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzrb(jz) * zgnz / zrmaj
!
        difthi(1,1,jz) = difthi(1,1,jz) + thirb(jz)
        difthi(2,2,jz) = difthi(2,2,jz) + thdrb(jz)
        difthi(3,3,jz) = difthi(3,3,jz) + therb(jz)
        difthi(4,4,jz) = difthi(4,4,jz) + thzrb(jz)
!
 
! ..................................
! .  the kinetic ballooning model  .
! ..................................
!
!       zbprim and zbc1 computed above under drift model
!
      if (  abs(cswitch(2)) .gt. zepslon                                 &  
     &   .and.  zgpr .gt. 0.0_R8) then
!
      zbprim = abs( zbeta * zgpr / zrmaj )
      zbcoef1 = 1.0_R8
      if ( abs(cswitch(8)) .gt. zepslon ) zbcoef1 = cswitch(8)
      zbc1   = zbcoef1 * abs(zshat)/(1.7_R8*zq**2*zrmaj)
      zelfkb = zelonf**cswitch(5)
!
        zfbthn = exp( min(abs(zlgeps),                                   &  
     &     max(-abs(zlgeps),cswitch(2)*(zbprim/zbc1 - 1._R8))) )
!
        zdk = abs( zsound * zrhos**2 * zfbthn * zgpr / zrmaj )
!
        thdkb(jz) = zdk*fkb(1)*zelfkb
        thekb(jz) = zdk*fkb(2)*zelfkb
        thikb(jz) = zdk*fkb(3)*zelfkb
        thzkb(jz) = zdk*fkb(4)*zelfkb
!
!  add to the fluxes
!
        vflux(1,jz) = vflux(1,jz) + thikb(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdkb(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + thekb(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzkb(jz) * zgnz / zrmaj
!
        difthi(1,1,jz) = difthi(1,1,jz) + thikb(jz)
        difthi(2,2,jz) = difthi(2,2,jz) + thdkb(jz)
        difthi(3,3,jz) = difthi(3,3,jz) + thekb(jz)
        difthi(4,4,jz) = difthi(4,4,jz) + thzkb(jz)
!
      endif
!
 300  continue
!
!
!   end of the main do-loop over the radial index, "jz"----------
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
