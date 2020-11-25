      subroutine weiland14 (                                             &  
     &   letain,   cetain,   lprint,   neq,      nout,     gnein         &  
     & , gnhin,    gnzin,    gtein,    gthin,    gtzin,    tauhin        &  
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein       &  
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb          &  
     & , ndim,     omega,    gamma,    difthi,   velthi,   chieff        &  
     & , vflux,    nmodes,   nerr )
 
 
 
!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt
!| \topmargin 0pt  \textheight 9.0in
!| \oddsidemargin 0pt \textwidth 6.5in
!|
!| \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!| \newcommand{\jacobian}{{\cal J}}
!|
!| \begin{document}
!|
!| \begin{center}
!| {\bf
!| Weiland Model for Transport Driven by
!| Toroidal Ion Temperature Gradient and \\
!| Trapped Electron Modes  \\
!| {\tt weiland14.tex} \\
!| \vspace{1pc}
!| Glenn Bateman and Arnold Kritz \\
!| Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem, PA 18015 \\
!| \vspace{1pc}
!| Jan Weiland, Hans Nordman and P{\"a}r Strand\\
!| Department of Electromagnetics \\
!| Chalmers University of Technology \\
!| S-412 96 G\"{o}teborg, Sweden \\
!| \vspace{1pc}
!| Jon Kinsey \\
!| General Atomics \\
!| P.O. Box 85608, San Diego, CA 92186} \\
!| \vspace{1pc}
!| \end{center}
!| This subroutine evaluates the transport matrix for $\eta_i$ and trapped
!| electron modes derived by Jan Weiland, H. Nordman and their group in
!| G\"{o}teborg Sweden \cite{bate98a} -- \cite{jarm87a}.
!| The equations in this routine include fast Hydrogenic ions, impurities,
!| trapped electron and finite Larmor radius effects. New options include
!| parallel ion motion, finite beta and collisional effects together with an
!| approximate treatment of ${\bf E}\times {\bf B}$ flow shear reduction of
!| transport.
!|
!| The dimensionless diffusivity matrix $ D_{j_1,j_2} = {\tt difthi(j1,j2)}$
!| and convective velocity array $ v_{j_1} = {\tt velthi(j1)} $
!| are given as:
!| $$ \frac{\partial}{\partial t}
!|  \left( \begin{array}{c} n_H T_H  \\ n_H \\ n_e T_e \\
!|     n_Z \\ n_Z T_Z \\ \vdots
!|     \end{array} \right)
!|  = - \nabla \cdot
!| \left( \begin{array}{l} {\rm vFlux}_1 n_H T_H \\
!|  {\rm vFlux}_2 n_H \\
!|  {\rm vFlux}_3 n_e T_e \\
!|  {\rm vFlux}_4 n_Z \\
!|  {\rm vFlux}_5 n_Z T_Z \\
!|  \vdots \end{array} \right)
!|  + \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  = \nabla \cdot
!| \left( \begin{array}{llll}
!| D_{1,1} n_H & D_{1,2} T_H & D_{1,3} n_H T_H / T_e \\
!| D_{2,1} n_H / T_H & D_{2,2} & D_{2,3} n_H / T_e \\
!| D_{3,1} n_e T_e / T_H & D_{3,2} n_e T_e / n_H & D_{3,3} n_e & \vdots \\
!| D_{4,1} n_Z / T_H & D_{4,2} n_Z / n_H & D_{4,3} n_Z / T_e \\
!| D_{5,1} n_Z T_Z / T_H & D_{5,2} n_Z T_Z / n_H &
!|         D_{5,3} n_Z T_Z / T_e \\
!|  & \ldots & & \ddots
!| \end{array} \right)
!|  \nabla
!|  \left( \begin{array}{c}  T_H \\ n_H \\  T_e \\
!|    n_Z \\  T_Z \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  + \nabla \cdot
!| \left( \begin{array}{l} {\bf v}_1 n_H T_H \\ {\bf v}_2 n_H \\
!|    {\bf v}_3 n_e T_e \\
!|    {\bf v}_4 n_Z \\ {\bf v}_5 n_Z T_Z \\ \vdots \end{array} \right) +
!|  \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right) $$
!| Note that all the diffusivities in this routine are normalized by
!| $ \omega_{De} / k_y^2 $,
!| convective velocities are normalized by $ \omega_{De} / k_y $,
!| and all the frequencies are normalized by $ \omega_{De} $.
!| \newpage
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Names of variables in the argument list:}
!| \begin{tabular}{lllp{3.0in}}
!| variable & status & symbol & meaning \\
!| {\tt letain(j)} & input & & Integer control variables
!|                             (see table below).\\
!| {\tt cetain(j)} & input & & Real-valued control variables
!|                             (see table below).\\
!| {\tt lprint}    & input & & Controls printout.
!|  Higher values produce more. \\
!| {\tt neq} & input & & number of equations \\
!| {\tt nout} & input & & output unit for error messages \\
!| {\tt gnein} & input & $  $ & $ - R\hat{r} \cdot \nabla n_e / n_e $ \\
!| {\tt gnhin} & input & $  $ & $ - R\hat{r} \cdot \nabla n_H / n_H $ \\
!| {\tt gnzin} & input & $  $ & $ - R\hat{r} \cdot \nabla n_Z / n_Z $ \\
!| {\tt gtein} & input & $  $ & $ - R\hat{r} \cdot \nabla T_e / T_e $ \\
!| {\tt gthin} & input & $  $ & $ - R\hat{r} \cdot \nabla T_H / T_H $ \\
!| {\tt gtzin} & input & $  $ & $ - R\hat{r} \cdot \nabla T_Z / T_Z $ \\
!| {\tt tauhin} & input & $\tau_H$ & $ \tau_H = T_H / T_e $ \\
!| {\tt tauzin} & input & $\tau_Z$ & $ \tau_Z = T_Z / T_e $ \\
!| {\tt fnzin} & input & $ f_{nZ} $ & $f_{nZ} = n_Z / n_e $ \\
!| {\tt czin}  & input & $ Z $ & impurity charge number \\
!| {\tt azin}  & input & $ m_Z / m_H $
!|      & impurity mass to hydrogen isotope mass. \\
!|  & & & Note that ``hydrogen'' may include a deuterium or tritium mix. \\
!| {\tt fnsin}  & input & $ f_s $
!|   & $ f_s = n_s / n_e $ fraction of superthermal hydrogenic ions \\
!| {\tt betaein} & input & $\beta_e$ &
!|      $ = n_e T_e / ( B^2 / 2 \mu_0 ) $ \\
!| {\tt ftrapein}  & input & $f_{trap} $ &
!|     fraction of trapped electrons \\
!| {\tt vef}       & input & $ \nu_{th} / \omega_{De} $ &
!|      thermal collision frequency, normalized \\
!| {\tt q}         & input & $ q $ & magnetic q-value \\
!| {\tt shear}     & input & $ s $ & $ d \ln q / d \ln r $ \\
!| {\tt ekyrhoin}  & input & $ k_y \rho_s $ & normalized poloidal
!|                     wave number \\
!| {\tt wexb}      & input & $\omega_{E\times B}$& $E\times B$ shearing rate
!| (normalized with $\omega_{D_e}$) \\
!| {\tt ndim} & input & & first dimension of the 2-D array difthi
!|                and the maximum number of unstable modes allowed \\
!| {\tt omega(j)}  & output & $\omega / \omega_{De} $ &
!|      real part of the frequencies normalized by $ \omega_{De} $ \\
!| {\tt gamma(j)}  & output & $\gamma / \omega_{De} $ &
!|      growth rates normalized by $ \omega_{De} $ \\
!| {\tt difthi(i,j)}      & output & $ D \omega_{De} / k_y^2 $
!|       & diffusivity matrix normalized by $ k_y^2 / \omega_{De} $ \\
!| {\tt velthi(j)}      & output & $ v \omega_{De} / k_y $
!|       & convective velocities normalized by $ k_y / \omega_{De} $ \\
!| {\tt chieff(j)} & output & $ \chi_{\rm eff} \omega_{De} / k_y^2 $
!|       & effective total diffusivities
!|         for $ n_H T_H $, $ n_H $, $ n_e T_e $,
!|         $ n_Z $, $ n_Z T_Z $, \ldots
!|         normalized by $ k_y^2 / \omega_{De} $ \\
!| {\tt vflux(j)}   & output & $ ({F}/{nT})({R k_y^2}/{\omega_{De}}) $
!|    & array containing computed fluxes \\
!| {\tt nmodes} & output & & number of unstable modes \\
!| {\tt nerr}       & output & & nerr $\neq 0 \rightarrow$ error\\
!| \end{tabular}
!| \end{center}
!| \newpage
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Integer control variables in the argument list:}
!| \begin{tabular}{lp{6.0in}}
!| {\bf variable} & {\bf meaning} \\
!| {\tt letain(2)} & = number of elements computed for transport matrix
!|                     (only when $> 0$) \\
!| {\tt letain(7)} &$ > 0 \rightarrow$ rescale transport matrix with velthi(j) = 0 \\
!| {\tt letain(9)} &$ > 1 \rightarrow$ do not produce transport matrix, only
!|                        effective diffusivities needed \\
!|                 & $ = 1 \rightarrow$ only diagonal elements of diffusion
!|                           matrix \\
!|                 & $ < 1 \rightarrow$ full diffusion matrix \\
!| {\tt letain(29)} & $ > 0 $ to print frequencies and fluxes mode by mode \\
!| \end{tabular}
!| \end{center}
!|
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Real-valued control variables in the argument list:}
!| \begin{tabular}{lp{6.0in}}
!| {\bf variable} & {\bf meaning} \\
!| {\tt cetain(10)} & = coefficient of $ k_\parallel $ [default to 1.0 in mmm95 model cswitch(17)] \\
!| {\tt cetain(11)}  & = 1./kpc, used to normalize $A_\parallel$ and $K$ [default to 1.0 in mmm95 model] \\
!| {\tt cetain(12)} & = coefficient in expression for $H$ [default to 0.0 in mmm95 model cswitch(19)] \\
!| {\tt cetain(15)} & = coefficient of $ \hat{\nu} $ [default to 0.0 in mmm95 model cswitch(18)] \\
!| {\tt cetain(20)} & = coefficient of $ \beta_{e,h,z} $ [default to 1.0 in mmm95 model cswitch(14)] \\
!| {\tt cetain(29)} & = radius used in printouts in weiland14flux\\
!| {\tt cetain(30)} & = finite difference used to construct
!|                    transport matrix \\
!| \end{tabular}
!| \end{center}
!|
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Effects included with different numbers of equations:}
!| \begin{tabular}{lp{5.0in}}
!| {\bf neq} & {\bf effects} \\
!| 2 & Hydrogenic Ion Temperature Gradient mode only \\
!| 4 & Hydrogenic ITG + trapped electron modes \\
!| 5 & Hydrogenic ITG with parallel ion motion + TEM \\
!| 6 & Hydrogenic + impurity ITG + TEM \\
!| 7 & Hydrogenic + impurity ITG + TEM + collisions \\
!| 8 & Hydrogenic ITG with parallel hydrogen ion motion
!|   + impurity ITG (without parallel impurity ion motion)
!|   + TEM + collisions\\
!| 9 & Hydrogenic ITG + ipurity ITG (without any parallel ion motion)
!|   + TEM + electromagnetic (finite $\beta$) effects \\
!| 10 & Hydrogenic ITG with parallel ion motion
!|   + impurity ITG without parallel ion motion
!|   + TEM + electromagnetic (finite $\beta$) effects + collisions \\
!| 11 & Hydrogenic ITG with parallel ion motion
!|   + impurity ITG with parallel ion motion
!|   + TEM + electromagnetic (finite $\beta$) effects + collisions \\
!|
!| \end{tabular}
!| \end{center}
!|
!| This routine {\tt weiland14} calls the routine {\tt weiland14flux}
!| to compute transport fluxes {\tt vflux} and the effective diffusivities
!| {\tt chieff}.
!| If a full transport matrix is requested ($ {\tt letain(9)} < 1 $),
!| derivatives of the transport fluxes are taken with respect to the
!| normalized gradients to compute {\tt difthi(i,j)}.
!| These diffusive fluxes are then subtracted from the transport fluxes
!| to compute the convective velocities {\tt velthi(j)}.
!|
!-----------------------------------------------------------------------
!
!  External dependencies:
!
!  Call tree: WEILAND14 calls the following routines
!
!    WEILAND14FLUX - Calculates fluxes and effective diffusivities
!      TOMSQZ      - Wrapper for QZ algorithm solving Ax = lambda Bx
!        CQZHES    - First step in QZ algorithm
!        CQZVAL    - Second and third step in QZ algorithm
!        CQZVEC    - Fourth step in QZ algorithm
!
!-----------------------------------------------------------------------
 
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer idp
      parameter (idp = 15)
 
      integer                                                            &  
     &   letain,   lprint, neq,    ndim,    nmodes,   nerr,   nout       &  
     & , imatrx,   j1,     j2,     j
!
! ndim  = first dimension of the 2-D array difthi
!         and the maximum number of unstable modes allowed
! nmodes = number of unstable modes
! imatrx= the number of elements computed along each row and
!          column of the transport matrix
 
      logical zdiffinit
      data zdiffinit /.true./
!
      dimension                                                          &  
     &   letain(32),        cetain(32),     omega(*),         gamma(*)   &  
     & , chieff(*),         vflux(*),       difthi(ndim,*),   velthi(*)
!
      REAL*8                                                             &  
     &   cetain,   gnein,  gnhin,  gnzin,   gtein,    gthin,  gtzin      &  
     & , tauhin,   tauzin, fnzin,  czin,    azin,     fnsin,  betaein    &  
     & , ftrapein, vef,    q,      shear,   ekyrhoin, wexb,   omega      &  
     & , gamma,    chieff, vflux,  difthi,  velthi
!
!
      REAL*8                                                             &  
     &   zgne,     zgnh,   zgnz,    zgns,    zgte,      zgth, zgtz       &  
     & , zone,     zdg,    ztemp,   zepsqrt, zepsmach,  zdgflux(idp)     &  
     & , ztwo,     zhalf,  chidum(idp),      zfluxmax
 
      save                                                               &  
     &   zone,     ztwo,   zhalf,   zepsqrt, zepsmach, zdiffinit
 
      if (zdiffinit) then
        zone  = 1.0_R8
        ztwo  = 2.0_R8
        zhalf = zone/ztwo
        zepsmach = zhalf
  2     if ( zhalf * zepsmach + zone .gt. zone ) then
          zepsmach = zhalf * zepsmach
          go to 2
        endif
        zepsqrt = sqrt(zepsmach)
        zdiffinit = .false.
      endif
 
!
!..initialize arrays
!
      do j1=1,ndim
        omega(j1)  = 0.0_R8
        gamma(j1)  = 0.0_R8
        chieff(j1) = 0.0_R8
        vflux(j1)  = 0.0_R8
        velthi(j1) = 0.0_R8
        do j2=1,ndim
          difthi(j1,j2) = 0.0_R8
        enddo
      enddo
!
      do j1=1,idp
        zdgflux(j1) = 0.0_R8
        chidum(j1)  = 0.0_R8
      enddo
!
!...Define Transport matrix sizes
!
      imatrx = min ( neq - 1, 4 )
      if ( letain(2) .gt. 0 ) imatrx = min ( imatrx, letain(2)-1 )
 
!
!...copy gradients to local variables
!
 
      zgne = gnein
      zgnh = gnhin
      zgnz = gnzin
      zgth = gthin
      zgte = gtein
      zgtz = gtzin
      zgns = zgne-zgnh*( zone-czin*fnzin-fnsin ) - zgnz*czin*fnzin
!
!... Establish fluxes
!
 
      call weiland14flux (                                               &  
     &   letain,   cetain,   lprint,   neq,      nout,     zgne          &  
     & , gnhin,    gnzin,    gtein,    gthin,    gtzin,    tauhin        &  
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein       &  
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb          &  
     & , ndim,     omega,    gamma,    chieff,   vflux,    nmodes        &  
     & , nerr )
!
!  return on error condition
!
      if ( nerr .ne. 0 ) then
         nerr = -7
         return
      endif
!
!  return if only chieff array is needed
!
      if (letain(9) .gt. 1) return
!
!  return if there is no flux
!
      zfluxmax = 0.0_R8
      do j = 1, imatrx + 1
        zfluxmax = zfluxmax + abs( vflux(j) )
      enddo
      if ( zfluxmax .lt. zepsqrt ) return
!
!... Define forward differencing step
!
      zdg = cetain(30)
      if ( abs(zdg) .lt. zepsqrt ) zdg = 0.01_R8
!
!... Take the derivative w.r.t gth
!
 
      zgth = zgth + sign(zdg,zgth)
      call weiland14flux (                                               &  
     &   letain,   cetain,   lprint,   neq,      nout,     zgne          &  
     & , gnhin,    gnzin,    gtein,    zgth,     gtzin,    tauhin        &  
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein       &  
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb          &  
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes        &  
     & , nerr )
!
!  return on error condition
!
      if ( nerr .ne. 0 ) then
         nerr = -7
         return
      endif
 
!
      do j = 1, imatrx +1
         difthi(j,1) = (zdgflux(j) - vflux(j))/sign(zdg,zgth)
      enddo
      zgth = gthin
 
!
!... Take derivatives w.r.t gnh
!
      zgnh = zgnh + sign(zdg,zgnh)
      zgne = zgnh*( zone-czin*fnzin-fnsin )+gnzin*czin*fnzin + zgns
 
      call weiland14flux (                                               &  
     &   letain,   cetain,   lprint,   neq,      nout,     zgne          &  
     & , zgnh,     gnzin,    gtein,    gthin,    gtzin,    tauhin        &  
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein       &  
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb          &  
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes        &  
     & , nerr )
!
!  return on error condition
!
      if ( nerr .ne. 0 ) then
         nerr = -7
         return
      endif
!
      do j = 1, imatrx +1
         difthi(j,2) = (zdgflux(j) - vflux(j))/sign(zdg,zgnh)
      enddo
      zgnh = gnhin
!
!... Take derivatives w.r.t gte
!
      zgte = zgte + sign(zdg,zgte)
      zgne = zgnh*( zone-czin*fnzin-fnsin )+gnzin*czin*fnzin + zgns
 
      call weiland14flux (                                               &  
     &   letain,   cetain,   lprint,   neq,      nout,     zgne          &  
     & , gnhin,    gnzin,    zgte,     gthin,    gtzin,    tauhin        &  
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein       &  
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb          &  
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes        &  
     & , nerr )
!
!  return on error condition
!
      if ( nerr .ne. 0 ) then
         nerr = -7
         return
      endif
!
      do j = 1, imatrx +1
         difthi(j,3) = (zdgflux(j) - vflux(j))/sign(zdg,zgte)
      enddo
      zgte = gtein
!
!... Take derivatives w.r.t gnz
!
 
      if (neq .gt. 4) then
 
        zgnz = zgnz + sign(zdg,zgnz)
        zgne = zgnh*( zone-czin*fnzin-fnsin ) + zgnz*czin*fnzin + zgns
 
        call weiland14flux (                                             &  
     &   letain,   cetain,   lprint,   neq,      nout,     zgne          &  
     & , gnhin,    zgnz,     gtein,    gthin,    gtzin,    tauhin        &  
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein       &  
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb          &  
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes        &  
     & , nerr )
!
!  return on error condition
!
      if ( nerr .ne. 0 ) then
         nerr = -7
         return
      endif
!
        do j = 1, imatrx + 1
           difthi(j,4) = (zdgflux(j) - vflux(j))/sign(zdg,zgnz)
        enddo
        zgnz = gnzin
 
!
!... Take derivatives w.r.t gtz
!
        zgtz = zgtz + sign(zdg,zgtz)
        zgne = zgnh*( zone-czin*fnzin-fnsin ) + zgnz*czin*fnzin + zgns
 
        call weiland14flux (                                             &  
     &   letain,   cetain,   lprint,   neq,      nout,     zgne          &  
     & , gnhin,    gnzin,    gtein,    gthin,    zgtz,     tauhin        &  
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein       &  
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb          &  
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes        &  
     & , nerr )
!
!  return on error condition
!
      if ( nerr .ne. 0 ) then
         nerr = -7
         return
      endif
!
        do j = 1, imatrx + 1
          difthi(j,5) = (zdgflux(j) - vflux(j))/sign(zdg,zgtz)
        enddo
        zgtz = gtzin
      else
        do j = 1,imatrx+1
          difthi(j,4) = 0.0_R8
          difthi(j,5) = 0.0_R8
        enddo
      endif
 
 
!
!... Convective velocities
!
      velthi(1) = vflux(1) - difthi(1,1) * gthin - difthi(1,2) * gnhin   &  
     &                     - difthi(1,3) * gtein - difthi(1,4) * gnzin
!
      velthi(2) = vflux(2) - difthi(2,1) * gthin - difthi(2,2) * gnhin   &  
     &                     - difthi(2,3) * gtein - difthi(2,4) * gnzin
!
      velthi(3) = vflux(3) - difthi(3,1) * gthin - difthi(3,2) * gnhin   &  
     &                     - difthi(3,3) * gtein - difthi(3,4) * gnzin
!
      if (neq .gt. 4) then
!
      velthi(4) = vflux(4) - difthi(4,1) * gthin - difthi(4,2) * gnhin   &  
     &                     - difthi(4,3) * gtein - difthi(4,4) * gnzin
!
      velthi(5) = vflux(5) - difthi(5,1) * gthin - difthi(5,2) * gnhin   &  
     &                     - difthi(5,3) * gtein - difthi(5,4) * gnzin
      else
         velthi(4) = 0.0_R8
         velthi(5) = 0.0_R8
      endif
!
!... Temporary lines Needed for correspondance with etaw14a and earlier
!
      do j = 1, imatrx+1
         difthi(j,5) = 0.0_R8
         difthi(5,j) = 0.0_R8
      enddo
      velthi(5) = 0.0_R8
!
!
!... Rescale diffusivity matrix and set velthi = 0 if that is requested
!    through letain(7) > 0
!
      if (letain(7) .gt. 0) then
         do j = 1, imatrx +1
            ztemp = vflux(j) - velthi(j)
            do j1 = 1, imatrx +1
                difthi(j,j1) = difthi(j,j1) * vflux(j) / ztemp
            enddo
            velthi(j) = 0.0_R8
         enddo
      endif
!
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
