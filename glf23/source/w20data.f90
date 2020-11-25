module w20data
  implicit none

  !Precision constants
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(10)
  INTEGER, PARAMETER :: CKIND = RKIND

  !Order of the dispersion relation
  integer, parameter :: neq = 9

  !Flags used for searching ion or electron drift modes
  integer :: searchmode
  integer, parameter :: S_ION = 2
  integer, parameter :: S_ELC = 3
  integer, parameter :: S_ANY = 0
  integer, parameter :: N_NONE = 0
  real, parameter :: grdmin = 0.001

  !Printing flag
  integer :: lprint

  !Species mass
  real(rkind) :: aimp    ! Atomic mass a_z (a.m.u.)
  real(rkind) :: ahyd    ! Atomic mass a_i (a.m.u.)

  !Species charge
  real(rkind) :: zimp    ! Impurity charge Z_z
  real(rkind) :: zhyd    ! Hydrogenic charge Z_i
  real(rkind) :: zfion   ! Fast ion charge Z_f
  real(rkind) :: zeff    ! Effective charge

  !Density and temperature gradients
  !Defined as -(R/U)(dU/dr), derivative respect to minor radius r (not rho)
  real(rkind) :: gte     ! Electron temperature gradient 
  real(rkind) :: gti     ! Ion temperature gradient
  real(rkind) :: gtz     ! Impurity temperature gradient
  real(rkind) :: gne     ! Electron density gradient
  real(rkind) :: gni     ! Ion density gradient
  real(rkind) :: gnz     ! Impurity density gradient
  real(rkind) :: gvt     ! Toroidal velocity gradient
  real(rkind) :: gvp     ! Poloidal velocity gradient
  real(rkind) :: gvpar   ! Parallel velocity gradient

  !Correlation lengths
  real(rkind) :: kyrho   ! k_y rho = 0.316, ( k_y rho )^2 = 0.1
  real(rkind) :: flh     ! Hydrogenic modes correlation length ( k_y rho )**2~0.10
  real(rkind) :: fle     ! Electron modes correlation length
  real(rkind) :: zflh    ! Atomic mass * Correlation length / hydrogenic charge squared
  real(rkind) :: zflz    ! Same, for impurities
  real(rkind) :: epsncorr ! epsilon_n correlation factor ( \epsilon_n = 2 L_n / L_B = 2 / gne )

  !Electron temperature and temperature ratios
  real(rkind) :: te      ! Te
  real(rkind) :: th      ! Th (hydrogenic ion temperature)
  real(rkind) :: tz      ! Tz (impurity temperature)
  real(rkind) :: tau     ! Te/Ti
  real(rkind) :: tau_inv ! Ti/Te
  real(rkind) :: tauh    ! Ti/Te
  real(rkind) :: tauz    ! Tz/Te

  real(rkind) :: ztauh   ! Z_i Ti/Te
  real(rkind) :: ztauz   ! Z_z Tz/Te

  !Electron density and density ratios
  real(rkind) :: ne      ! Electron density n_e
  real(rkind) :: nh      ! Main hydrogenic ion density n_h
  real(rkind) :: nz      ! Impurity ion density n_z
  real(rkind) :: ns      ! Fast particle density n_s
  real(rkind) :: znz_ne  ! Z_imp nz/ne
  real(rkind) :: zns_ne  ! Z_(fast ion) ns/ne
  real(rkind) :: ne_nh   ! ne/nh
  real(rkind) :: nh_ne   ! ni/ne (fraction of main ions)
  real(rkind) :: nz_ne   ! nz/ne (fraction of impurity ions)
  real(rkind) :: ns_ne   ! ns/ne (fraction of fast ions)
  real(rkind) :: fte     ! n_e (trapped) / n_e

  !Fractions
  real(rkind) :: fft
  real(rkind) :: fzft
  real(rkind) :: ftrt
  real(rkind) :: gm
  real(rkind) :: bta
  real(rkind) :: xt

  !Collision frequency related
  real(rkind) :: vei     ! Electron-ion collision rate
  real(rkind) :: vef     ! nu effective
  real(rkind) :: bt
  real(rkind) :: bt1

  !Geometry and equilibrium quantities
  real(rkind) :: rhos    ! Gyroradius
  real(rkind) :: rmaj    ! Local major radius
  real(rkind) :: rmaj0   ! Major radius at magnetic axis
  real(rkind) :: rmin    ! Local minor radius
  real(rkind) :: amin    ! Plasma minor radius a
  real(rkind) :: eps0    ! Global aspect ratio a/R_0
  real(rkind) :: eps     ! Local aspect ratio r/R
  real(rkind) :: btor    ! Local toroidal magnetic field
  real(rkind) :: q       ! Magnetic safety factor
  real(rkind) :: shear   ! Magnetic shear (r/q)(dq/dr)
  real(rkind) :: shat    ! Effective magnetic shear modified for elongation
  real(rkind) :: kappa   ! Plasma elongation
  real(rkind) :: alpha_MHD   ! MHD alpha

  !Plasma beta
  real(rkind) :: beta    ! Plasma beta
  real(rkind) :: betae   ! Electron beta

  !ExB shear rates
  real(rkind) :: wexb    ! Shear rate w_ExB

  !Saved internal variables used in dispersion relation

  real(rkind) :: str = 7.0/3.0
  real(rkind) :: ftr = 5.0/3.0   
  real(rkind) :: tvr = 2.0/3.0
  real(rkind) :: em2 = 1.0
  real(rkind) :: em1 = 1.0
  real(rkind) :: em  = 1.0
  real(rkind) :: kpc = 1.0

  !Plasma velocity normalized by the local ion soundspeed Cs = sqrt ( Te/mi )
  real(rkind) :: vtor       ! Toroidal velocity
  real(rkind) :: vpol       ! Poloidal velocity
  real(rkind) :: vpar       ! Parallel velocity

  !Some other useful quantities
  real(rkind) :: wde        ! 2.0 * kyrho**2 * csound / R
  real(rkind) :: Csound     ! speed of sound
  real(rkind) :: Csound0    ! speed of sound at magnetic axis
  real(rkind) :: diffBohm   ! Bohm diffusion coefficient

  contains

    subroutine w20init( i_print, &
         z_te,   z_ne,    z_vtor, z_vpol,   z_vpar, &
         z_btor, z_amin,  z_rmin, z_rmaj,   z_eps0, &
         z_aimp, z_ahyd,  z_zimp, &
         z_gte,  z_gti,   z_gtz , &
         z_gne,  z_gni,   z_gnz, &
         z_gvt,  z_gvp,   z_gvpar, &
         z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
         z_q,    z_shear, z_kappa, z_wexb , z_Csound0 )

      implicit none

      real(rkind), intent(in) :: &
           z_te,   z_ne,   z_vtor, z_vpol, z_vpar,  &
           z_btor, z_amin, z_rmin, z_rmaj, z_eps0,  &
           z_aimp, z_ahyd, z_zimp, &
           z_gte,  z_gti,  z_gtz , &
           z_gne,  z_gni,  z_gnz,  &
           z_gvt,  z_gvp,  z_gvpar, &
           z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
           z_q,    z_shear, z_kappa, z_wexb, z_Csound0

      integer, intent(in) :: i_print
      
      lprint = i_print
      te     = z_te
      ne     = z_ne*1.0E-19
      btor   = z_btor
      amin   = z_amin
      rmin   = z_rmin
      rmaj   = z_rmaj
      eps0   = z_eps0
      aimp   = z_aimp
      ahyd   = z_ahyd
      zimp   = z_zimp
      kyrho  = z_kyrho
      tauh   = z_tauh
      tauz   = z_tauz
      nz_ne  = z_nz_ne
      ns_ne  = z_ns_ne
      fte    = z_fte
      q      = z_q
      shear  = z_shear
      kappa  = z_kappa
      Csound0= z_Csound0 
      
      !Initialize some additional variables
      zhyd  = 1.0
      zfion = 1.0

      !Aspect ratio and major radius at axis

      rmaj0 = amin / eps0
      eps   = max( 0.01, rmin / rmaj0 )
!      rmaj  = rmaj0 * ( 1.0+1.084*eps )

      !Gradients
      gte    = sign( max( abs(z_gte), grdmin ), z_gte ) * rmaj / rmaj0
      gti    = sign( max( abs(z_gti), grdmin ), z_gti ) * rmaj / rmaj0
      gtz    = sign( max( abs(z_gtz), grdmin ), z_gtz ) * rmaj / rmaj0
      gne    = sign( max( abs(z_gne), grdmin ), z_gne ) * rmaj / rmaj0
      gni    = sign( max( abs(z_gni), grdmin ), z_gni ) * rmaj / rmaj0
      gnz    = sign( max( abs(z_gnz), grdmin ), z_gnz ) * rmaj / rmaj0
      gvp    = sign( max( abs(z_gvp), grdmin ), z_gvp ) * rmaj / rmaj0
      gvt    = sign( max( abs(z_gvt), grdmin ), z_gvt ) * rmaj / rmaj0
      gvpar  = sign( max( abs(z_gvpar), grdmin ), z_gvpar ) * rmaj / rmaj0

!      gte    = max(z_gte, grdmin ) * rmaj / rmaj0
!      gti    = max(z_gti, grdmin ) * rmaj / rmaj0
!      gtz    = max(z_gtz, grdmin ) * rmaj / rmaj0
!      gne    = max(z_gne, grdmin ) * rmaj / rmaj0
!      gni    = max(z_gni, grdmin ) * rmaj / rmaj0
!      gnz    = max(z_gnz, grdmin ) * rmaj / rmaj0
!      gvp    = max(z_gvp, grdmin ) * rmaj / rmaj0
!      gvt    = max(z_gvt, grdmin ) * rmaj / rmaj0
!      gvpar  = max(z_gvpar, grdmin ) * rmaj / rmaj0


      !Effective magnetic shear
      shat = sqrt(max(2*shear-1.0+(kappa*(shear-1))**2,0.0))

      !Temperatures and temperature ratios
      th = te*tauh
      tz = te*tauz

      tau      = te  / max( th, 0.001 )    ! Te/Ti
      tau_inv  = th  / max( te, 0.001 )    ! Te/Ti

      ztauh = tauh / zhyd               ! T_h / ( Z_h Te )
      ztauz = tauz / zimp               ! T_z / ( Z_z Te )

      !Correlation lengths (k_y rho )**2
      epsncorr = min( 2.0_rkind / abs(gne), 4.0 )**0.17_rkind
      flh   = (0.7+2.4/(7.14*q*shat+0.1))*kyrho**2.0*2.0/(1.0+tau_inv)
!      zflh  = flh
      kyrho = (dsqrt(flh)-0.1*fte)*epsncorr
      zflh  = kyrho**2.0 / zhyd**2.0
      zflz  = aimp * kyrho**2.0 / zimp**2.0    ! A_z * ( k_y rho ) / Z_z**2.0

      znz_ne = z_zimp*nz_ne
      zns_ne = zfion*ns_ne

      !Densities and density ratios
      nh_ne = (1.0 - znz_ne - zns_ne )
      ne_nh = 1.0 / nh_ne
      nh    = ne * nh_ne
      nz    = ne * nz_ne
      ns    = ne * ns_ne

      !Effective charge Zeff
      zeff = (nh_ne*zhyd**2 + nz_ne*zimp**2 + ns_ne*zfion**2 ) 

      fft   = (1.0 - znz_ne - zns_ne ) / ( 1.0 - fte )
      fzft  = znz_ne / ( 1.0 - fte )
      ftrt  = ( 5.0/3.0 ) * tauh

      !Collision related quantities
      bt    = 1.5
      bt1   = bt - 2.5

      !Plasma beta
      beta = 4.027E-22 * ( z_ne*te*(1.0 + nh_ne*tauh + nz_ne*tauz ) ) / (btor**2)

      !Electron beta
!      betae = 0.01D0*0.4092D0*ne*te/(btor**2)
      betae = 4.027E-22 * z_ne*te / (btor**2)      

      !MHD alpha
      alpha_MHD = em*(q**2.0)*betae*( gne + gte + tauh*(gni + gti) )

      !Speed of sound
      Csound = sqrt( te * 1000.0 * 1.602E-19 / ( 1.6725E-27 * ahyd ) )
!      Csound = 0.311*dsqrt(te)*10**6.0/dsqrt(ahyd)

      !Ion gyroradius
      rhos = Csound / ( 1.D8*0.957*btor )
!      rhos = (4.57E-3 * sqrt(ahyd * te*tauh)) / btor 

      !Temporary: use Csound at axis for momentum transport stuff
!      Csound = Csound0

      !Toroidal velocity normalized by speed of sound
      vtor = z_vtor / Csound0

      !Poloidal velocity normalized by speed of sound
      vpol = z_vpol / Csound0

      !Parallel velocity normalized by speed of sound
      vpar = z_vpar / Csound0

      !Magnetic drift frequency with fixed k_y rho
      wde = 2.0 * dsqrt(flh) * Csound / rmaj
!      wde = 2.0 * kyrho * Csound / rmaj
!      wde = 2.0 * kyrho * Csound / (rmaj*(1.0+1.084*eps))

      !ExB shear rates
      wexb   = abs(z_wexb / wde)

      !Bohm diffusion coefficient
      diffBohm = rhos * Csound

      !Electron ion collision rate (main ions only)
      vei = 0.09*1.0E4 * nh * ( 15.95- dlog( dsqrt(ne)/te) ) / te**1.5
!      vei = 0.09174*1.0E4 * nh * (15.2 - 0.5*dlog(ne/10) + dlog(te)) / te**1.5
!      vei = (4.0/3.0*sqrt(2.0*3.14159265))*nh*(37.8-dlog(dsqrt(ne)/te))/te**1.5

      !Effective collision rate
      vef = vei / ( eps * wde )
!      vef = vei * q / ( eps * wde )

      !Electron correlation length
      fle = max(0.005,0.024-0.0064*(vei/10000.0-0.61)+0.4*betae)

      gm   = 1.0 / (1.0-fte)
      bta  = ftrt*(gm+tauh)
      xt   = 1.0 / (1.0+tauh)

!      lprint = 9
      lprint = i_print

    end subroutine w20init

    !This function returns the position of the mode with
    !the largest growthrate
    !
    !INPUT: zz(nmod): eigenvalues of drift mode dispersion relation )
    !       nmode   : dimension of the vector
    !       S_SRCH  : direction S_ION or S_ELC or S_ANY
    !
    !OUTPUT: nmod_ptr : position of the largest growing mode in
    !                   direction S_SRCH in array zz
    !                   will return 0 if no 
    function w20wmunst ( zz, nmod, S_SRCH ) result ( nmod_ptr )
      implicit none

      integer, intent(in) :: nmod
      complex(ckind), intent(in), dimension(nmod) :: zz
      integer, intent(in) :: S_SRCH

      integer :: nmod_ptr, i, S_FND 
      real(rkind) :: gamma_max

      gamma_max = -10000.0
      S_FND = 0
      nmod_ptr = N_NONE

      do i=1,nmod

         if ( S_SRCH .ne. S_ANY ) then
            S_FND = w20s_fnd( zz(i) )
         end if
         
         if ( dimag(zz(i)) > gamma_max .and. S_FND .eq. S_SRCH ) then
            gamma_max = dimag(zz(i))
            nmod_ptr = i
         end if
      end do

    end function w20wmunst

    !This function takes in a complex frequency and
    !evaluates whether the mode is an electron or ion mode
    !
    ! INPUT: zz = complex frequency ( w_r, gamma )
    ! OUTPUT: S_FND = S_ELC or S_ION ( integer flag for electron or ion mode)
    !
    function w20s_fnd ( zz ) result ( S_FND )
      implicit none

      complex(ckind), intent(in) :: zz

      integer :: S_FND

      if ( dreal(zz) > 0.0 ) then
         S_FND = S_ELC
      else
         S_FND = S_ION
      end if

    end function w20s_fnd

    !This function takes in a set of complex frequencies
    !and counts the number of unstable modes in any direction
    !
    ! INPUT: zz(nmod) = complex frequency
    !        nmod     = dimension of zz vector
    !        S_SRCH   = direction of search S_ANY, S_ION, S_ELC
    !
    ! OUTPUT: iunst   = number of unstable modes in S_SRCH direction

    function w20nunst ( zz, nmod, S_SRCH ) result ( iunst )
      implicit none
      
      integer, intent(in) :: nmod, S_SRCH
      complex(ckind), intent(in), dimension(nmod) :: zz
      
      integer :: iunst, i, S_FND
      
      iunst = N_NONE
      S_FND = S_ANY
      
      do i=1,nmod
         
         if ( S_SRCH .ne. S_ANY ) then
            S_FND = w20s_fnd( zz(i) )
         end if
         
         if ( dimag(zz(i)) > 0.0 .and. S_SRCH .eq. S_FND ) then
            iunst = iunst + 1
         end if
         
      end do
     
    end function w20nunst

    !This function returns the growthrate of the drift mode
    !at the position nmod_ptr. If nmod_ptr is not between
    !1 and nmod, then the result is zero. This is useful
    !when no unstable modes are found in one of the functions
    !above
    !
    !INPUT: zz(nmod) = complex eigenfrequencies
    !       nmode    = dimension of zz vector
    !       nmod_ptr = position of the growthrate
    !
    function w20gamma ( zz, nmod, nmod_ptr ) result ( gamma )
      implicit none
      
      integer, intent(in) :: nmod, nmod_ptr
      complex(ckind), intent(in), dimension(nmod) :: zz

      real(rkind) :: gamma
      integer :: i

      if ( nmod_ptr .eq. N_NONE ) then
         gamma = 0.0
      else if ( nmod_ptr .ge. 1 .and. nmod_ptr .le. nmod ) then
         gamma = dimag(zz(nmod_ptr))
      else
         write (*,*) "Invalid mode number: ",nmod_ptr
      end if

    end function w20gamma

    !This function takes in the real and imaginary parts of the
    !eigenvalues and folds them into a complex frequency
    !The frequencies are then normalized

    function w20omg ( zomega, zgamma, zbeta, nmod ) result ( zz )
      implicit none

      integer, intent(in) :: nmod
      real(rkind) , intent(in), dimension(nmod) :: zomega, zgamma, zbeta

      complex(ckind), dimension(nmod) :: zz
      real, dimension(nmod) :: zb

      zz = (0.0,0.0)
      zb = max(1.0E-4,zbeta)
      
      zz(1:nmod) = ( zomega(1:nmod)+(0.0,1.0)*zgamma(1:nmod) ) / zb(1:nmod)

    end function w20omg

    !This function computes the relative difference between two modes

    function w20delwz ( wz1, wz2 ) result ( delwz )
      implicit none
      
      complex(ckind) :: wz1, wz2

      real(rkind) :: delwz,denom
      real(rkind) :: eps = 0.001
      real(rkind) :: wr1, wr2, wi1, wi2

      eps = 0.001
      wr1 = dreal(wz1)
      wr1 = sign( max( abs(wr1), eps ), wr1)
      wi1 = dimag(wz1)
      wi1 = sign( max( abs(wi1), eps ), wi1)

      wr2 = dreal(wz2)
      wi2 = dimag(wz2)

      denom = ( max( abs(wz1), eps ) )**2
      denom = 0.5 * ( abs(wz1) + abs(wz2) )

      delwz = dsqrt( ( wr1-wr2 )**2 + ( wi1-wi2 )**2 ) / denom
!      delwz = ( (1.0-wr2/wr1)**2 + (1.0-wi2/wi1)**2 ) ! / denom 



    end function w20delwz

    !This function returns the mode closest to the input mode
    function w20wmatch ( wz, zz, nmod ) result ( nmod_ptr )
      implicit none

      integer, intent(in) :: nmod
      complex(ckind), intent(in) :: wz
      complex(ckind), intent(in), dimension(nmod):: zz

      integer :: nmod_ptr,i
      real :: delwz, delwzmin

      nmod_ptr = N_NONE
      delwzmin = 10000000.0
      delwz = 0.0

      do i=1,nmod
         delwz = w20delwz( wz, zz(i) )
         if (  delwz .lt. delwzmin ) then
            delwzmin = delwz
            nmod_ptr = i
         end if
      end do

    end function w20wmatch

end module w20data
