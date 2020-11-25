!For definition of input parameters look at file w20_README.txt
!or at the variable definitions in w20data.f90

Subroutine w20main ( i_print, &
     z_te,   z_ne,    z_vtor, z_vpol, z_vpar, &
     z_btor, z_amin,  z_rmin, z_rmaj, z_eps0, &
     z_aimp, z_ahyd,  z_zimp, &
     z_gte,  z_gti,   z_gtz , &
     z_gne,  z_gni,   z_gnz, &
     z_gvt,  z_gvp,   z_gvpar, &
     z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
     z_q,    z_shear, z_kappa, z_wexb, z_Csound0, &
     diffs,  vconv )

  use w20data
  implicit none

  !Variable declation
  integer, intent(in) :: i_print
  
  real(rkind), intent(in) :: &
       z_te,   z_ne,    z_vtor,  z_vpol,  z_vpar, &
       z_btor, z_amin,  z_rmin,  z_rmaj,  z_eps0, &
       z_aimp, z_ahyd,  z_zimp,  &
       z_gte,  z_gti,   z_gtz ,  &
       z_gne,  z_gni,   z_gnz,   &
       z_gvt,  z_gvp,   z_gvpar, & 
       z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
       z_q,    z_shear, z_kappa, z_wexb,  z_Csound0

  integer :: ierr, iiu, ieu, iimax, iemax, eiu, eeu, eimax, eemax, j

  !Resulting fluxes
  real(rkind), intent(inout), dimension(6) :: diffs, vconv

  !Individual fluxes from ion and electron call
  real(rkind), dimension(6) :: idiffs, ivconv, ediffs, evconv

  !Eigenvectors resulting from solved system of equations
  real(rkind), dimension(neq,neq) :: izvr, izvi, ezvr, ezvi, zvr, zvi

  !Complex frequency of the mode 
  complex(ckind), dimension(neq) :: izz, ezz, zz

  !Geometric average factor multiplying W_de, used to compute momentum thermoelectric pinch
  real(rkind) :: G_ave_i, G_ave_e

  !Parallel wavenumber of drift mode, used to compute momentum pinches
  real(rkind) :: kps_i, kps_e

  !Logical variable decides if separate correlation length needed
  logical :: l_elc

  !Initialize all internal variables
  call w20init( i_print, &
     z_te,   z_ne,    z_vtor,  z_vpol,  z_vpar, &
     z_btor, z_amin,  z_rmin,  z_rmaj,  z_eps0, &
     z_aimp, z_ahyd,  z_zimp, &
     z_gte,  z_gti,   z_gtz,  &
     z_gne,  z_gni,   z_gnz,  &
     z_gvt,  z_gvp,   z_gvpar,&
     z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
     z_q,    z_shear, z_kappa, z_wexb,  z_Csound0 )


  ! Algorithm starts here

  !Initialize geometric factor to 1.0 (circular)
  G_ave_e = 1.0
  G_ave_i = 1.0

  !Search for most unstable ion mode
  searchmode = S_ION

  !Variable correlation length for ion modes

  call w20solv( izz, izvr, izvi, iiu, ieu, iimax, iemax, G_ave_i, kps_i, ierr )

  !Search for most unstable electron mode using varying correlation length

  searchmode = S_ELC

  !Different correlation length for electrons
  zflh = fle
!  zflh = ahyd * fle / zhyd**2.0
  zflz = aimp * fle / zimp**2.0
  kyrho = sqrt(fle)
  wde  = 2.0 * kyrho * csound / rmaj

  call w20solv( ezz, ezvr, ezvi, eiu, eeu, eimax, eemax, G_ave_e, kps_e, ierr )

  !Compute transport coefficients

  diffs = 0.0
  vconv = 0.0

  idiffs = 0.0
  ediffs = 0.0

  ivconv = 0.0
  evconv = 0.0

  !G_ave_i = 1.0

  !No unstable modes, then no transport
  if ( iiu .eq. 0 .and. ieu .eq. 0 .and. eiu .eq. 0 .and. eeu .eq. 0 ) then
     return
  end if

  !Consider whether a separate correlation length for electrons is needed
  !Conditions for yes:
  ! Fastest growing mode in electron w20solv call gives more transport than
  ! fastest growing mode in ion w20solv call
  ! (note that unstable modes must exist for the comparison to occur)
  ! No electron modes in ion w20solv call
  !Conditions for no
  ! All others
  if ( eeu > 0 ) then

     if ( ieu > 0 ) then
        if ( dimag(izz(iemax))/sqrt(flh) .lt. dimag(ezz(eemax))/sqrt(fle) ) then
           l_elc=.true.
        else
           l_elc=.false.
        end if
     else
        l_elc=.true.
     end if

  else if ( eeu .eq. 0 ) then
     l_elc=.false.
  end if

  if ( l_elc ) then
    
     !Set up call to flux computation using only
     ! -Eigenvector associated with fastest growing electron mode in electron w20solv call
     ! -Eigenvalue associated with fastest growing electron mode in electron w20solv call
     zz = ezz
     ezz= (0.0,0.0)
     ezz(1) = zz(eemax)
     
     ezvr(1:neq,1) = ezvr(1:neq,eemax) 
     ezvi(1:neq,1) = ezvi(1:neq,eemax) 

     !Use different correlation length for electron modes
     zflh = fle
!     zflh = ahyd * fle / zhyd**2.0
     zflz = aimp * fle / zimp**2.0
     kyrho = sqrt(fle)
     wde  = 2.0 * kyrho * csound / rmaj

     !Compute fluxes
     call w20diff( ezz, ezvr, ezvi, G_ave_i, kps_e, ediffs, evconv )

     !Set up call to flux computation with ion modes
     ! -Remove fastest growing electron mode in ion w20solv call
     if (iemax>0) izz(iemax)=(0.0,0.0)

  end if

  !Use correlation length for ion modes
  kyrho = (dsqrt(flh)-0.1*fte)*epsncorr
  zflh = kyrho**2.0 / zhyd**2.0
!  zflh = ahyd * flh / zhyd**2.0
  zflz = aimp * kyrho**2.0 / zimp**2.0

  wde  = 2.0 * kyrho * csound / rmaj

  !Compute fluxes
  call w20diff( izz, izvr, izvi, G_ave_i, kps_i, idiffs, ivconv )

  !Add transport from two calls
  diffs = idiffs + ediffs
  vconv = ivconv + evconv

End Subroutine w20main
