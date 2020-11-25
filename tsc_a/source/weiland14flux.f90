!#include "f77_dcomplx.h"
      subroutine weiland14flux (                                         &  
     &   letain,   cetain,   lprint,   neq,      nout,     gnein         &  
     & , gnhin,    gnzin,    gtein,    gthin,    gtzin,    tauhin        &  
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein       &  
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb          &  
     & , ndim,     omega,    gamma,    chieff,   vfluxout, nmodes        &  
     & , nerr )
!@weiland14flux
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
!..See description of the argument list at the beginning of this file
!
!|
!| This subroutine evaluates the eigenvalue problem for $\eta_i$
!| and trapped electron modes derived by
!| Jan Weiland, H. Nordman and their group in G\"{o}teborg Sweden
!| \cite{bate98a} -- \cite{jarm87a}.
!| Transport fluxes are then computed using a quasi-linear estimate.
!| The equations in this routine include fast Hydrogenic ions,
!| impurities, trapped electron and finite Larmor radius effects.
!| New options include parallel ion motion, finite beta and collisional effects
!| together with an approximate treatment of ${\bf E}\times {\bf B}$ flow
!| shear reduction of transport.
!|
!| The essential idea is to linearize the fluid equations
!| with magnetic drifts for a given Fourier harmonic
!| \[ u = u_0(x) + \tilde{u} \exp\left[ -i (\omega t
!|    + \vec{k} \cdot \vec{x} ) \right] \]
!| Compute the eigenvalues and eigenvectors from these equations and then
!| compute the quasi-linear heat and particle fluxes.
!|
!|
!| The fundamental equations used in this model are the fluid equations.
!| For ions (hydrogen or impurity), the equation of continuity is
!| \[ \Partial{n_i}{t} + \nabla \cdot ( n_i v_i ) = 0 \]
!| where
!| \[ v_i = v_E + v_{*i} + v_{Pi} + v_{\pi i}
!|    + \hat{B} v_{\parallel i} \]
!| are the fluid flows where
!| \[ v_E = E \times \hat{B} / B = \hat{B} \times \nabla \phi / B
!|  = - i k_y \phi \hat{x} / B \]
!| is the E cross B drift;
!| \[ v_{*i} = \frac{ \hat{B} \times \nabla ( n_i T_i )}{Z_i e n_i B}
!|  = i ( 1 + \eta_i ) \omega_{*i} / k_y \qquad{\rm where}\qquad
!|  \omega_{*i} = \frac{ - k_y T_i g_{ni} }{ Z_i e B R } \]
!| is the diamagnetic drift ($ Z_i = 1 $ for hydrogen isotopes and
!| $ Z_i = -1 $ for electrons);
!| \[ g_{ni} \equiv - R \hat{x} \cdot \nabla n_i / n_i \]
!| is the density gradient scale length
!| \[ \eta_i \equiv g_{Ti} / g_{ni}
!|  = \frac{n_i \hat{x} \cdot \nabla T_i}{T_i \hat{x} \cdot \nabla n_i}, \]
!| \[ v_{Pi} = \frac{ d E }{ dt } / ( B \Omega_i ) \]
!| is the polarization drift with
!| \[ \Omega_i = \frac{ Z_i e B }{ m_i}; \]
!| and
!| \[ v_{\pi i} = \frac{ \hat{B} \times \nabla \cdot \pi_i}{
!|   Z_i e n_i B } \]
!| is the drift due to the off-diagonal elements of the stress tensor, $\pi_i$.
!|
!| The parallel ion motion $ v_{\parallel i} $ is determined by the parallel
!| ion momentum equation driven by electromagnetic forces and ion pressure
!| gradient along the field lines\cite{weil92a}:
!| \begin{equation}
!|   m_i n_i \Partial{v_{\parallel i}}{t}
!|  =  - e Z_i n_i \left[ \nabla_{\parallel} \phi + \frac{1}{c} \left(
!|     \Partial{A_\parallel}{t}
!|      - ( {\mathbf v}_{*i} \times {\mathbf B}_{\perp} ) \cdot \hat{{\mathbf B}}
!|          \right) \right]  - \nabla_{\parallel} ( n_i T_i )
!| \end{equation}
!| where $A_\parallel$ is the parallel component of the vector potential and
!| $ {\mathbf B}_{\perp} = {\mathbf \nabla} \times
!|     ( A_\parallel \hat{{\mathbf B}} )
!|    \approx {\mathbf \nabla} A_\parallel \times  \hat{{\mathbf B}} $.
!|
!| \noindent
!| The ion energy balance equation is
!| \[ \frac{3}{2} n_i \left( \Partial{ }{t} + \vec{v}_i \cdot \nabla \right) T_i
!|   + n_i T_i \nabla \cdot \vec{v}_i = - \nabla \cdot \vec{q}_{*i}
!|   = \frac{5}{2} n_i ( \vec{v}_{*i} - \vec{v}_{Di} ) \cdot \nabla T_i \]
!| where $ \vec{q}_{*i} $ is the diamagnetic ion heat flow and
!| \[ \vec{v}_{Di} = \frac{T_i}{m_i \Omega_i} \hat{B} \times
!|   \left( \frac{\nabla B}{B} + \vec{\kappa} \right) \]
!| is the drift due to $\nabla |B|$ and magnetic curvature
!| $ \vec{\kappa} = \hat{B} \cdot \nabla \hat{B} $.  Note
!| \[ \vec{k} \cdot \vec{v}_{Di} = \omega_{Di} = \frac{-2 k_y T_i}{Z_i e B R} \]
!| is the diamagnetic drift frequency.
!|
!| The equations for trapped electron continuity and energy flow have
!| the same form except that the polarization drift and the drift
!| due to the stress tensor
!| can be neglected and the electron fluid velocity can be written
!| \[ v_e = v_E +  v_{*e}\]
!| where
!| \[ v_{*e} = -\frac{ \hat{B} \times \nabla ( n_e T_e )}{ e n_e B}
!|  = i ( 1 + \eta_e ) \omega_{*e} / k_y \qquad{\rm where}\qquad
!|  \omega_{*e} = \frac{  k_y T_e g_{ne} }{  e B R } \]
!|
!| Let $f_t$ be the fraction of trapped electrons, $ f_t = n_{et} / n_e $,
!| where $n_{et}$ is the density of trapped electrons and
!| $n_{ef}$ and the density of free (circulating) electrons.  Then
!| \[ n_e = n_{et} + n_{ef}. \]
!| The perturbed density can be written
!| \[ \frac{\tilde{n}_e}{n_e}
!|  =  \frac{\tilde{n}_{et}}{n_e} + \frac{\tilde{n}_{ef}}{n_e}
!|  =  f_t \frac{\tilde{n}_{et}}{n_{et}}
!|        + ( 1 - f_t ) \frac{\tilde{n}_{ef}}{n_{ef}} \]
!| %  =  f_t \frac{\tilde{n}_{et}}{n_{et}}
!| %        + ( 1 - f_t ) \frac{ e \tilde{\phi} }{T_e}. \]
!|
!| Let $ f_Z $ be the fraction of impurity ions with charge $Z$
!| relative to the electron density $ f_Z = n_Z / n_e $.
!| Charge neutrality gives
!| \[ n_e = n_H + Z n_Z + n_s \]
!| where $n_H$ is the density of the hydrogen isotopes and
!| $n_Z$ is the density of the impurity species and $n_s$ is the density of
!| superthermal hydrogenic ions.  Then
!| \[ \frac{\tilde{n}_e}{n_e}
!|  =  \frac{\tilde{n}_H}{n_e} + \frac{Z \tilde{n}_Z}{n_e}
!|  =  ( 1 - Z f_Z - f_s ) \frac{\tilde{n}_H}{n_H}
!|        + Z f_Z \frac{\tilde{n}_Z}{n_Z}  \]
!| where $ f_s \equiv n_s / n_e $.
!| For now, it is assumed that the superthermal ions
!| do not take part in the perturbation $ \tilde{n}_s = 0.$
!|
!| A relation between the perturbed free (circulating) electron density
!| $ n_{ef} $ and the perturbed electric and magnetic potentials
!| can be obtained from the momentum equation for free electrons parallel
!| to the unperturbed magnetic field
!| \begin{equation}
!|  m_e \frac{\partial v_{\parallel e}}{\partial t}  =
!|    e \left( \hat{{\mathbf B}} \cdot {\mathbf \nabla}\phi
!|    + \frac{1}{c} \frac{\partial A_{\parallel}}{\partial t} \right)
!|    - \frac{e}{c} \left( {\mathbf v}_{*ep} \times \delta {\mathbf B}_{\perp}
!|       \right) \cdot \hat{{\mathbf B}}
!|    - \frac{1}{n} \hat{{\mathbf e}}_{\parallel} \cdot {\mathbf \nabla} p
!| \end{equation}
!| where
!| the unperturbed perpendicular free electron velocity is the
!| diamagnetic velocity (including the full electron pressure gradient)
!| $ {\mathbf v}_{*ep} = - \hat{{\mathbf B}} \times {\mathbf \nabla}
!|    ( n_{ef} T_e ) / ( e n_{ef} B ) $.
!| The inertial term on the left of this momentum equation can be neglected
!| if $ k_\parallel v_{the} \gg \omega $.
!|
!| Assuming the electron temperature is nearly uniform along
!| the perturbed magnetic field
!| $ {\mathbf B} \cdot {\mathbf \nabla} T_e \approx 0 $,
!| if follows that perturbed electron temperature is approximately
!| \begin{equation}
!|  \hat{T}_e =
!|    \frac{ \omega_{*ep} }{ k_{\parallel} c }
!|    \frac{ g_{Te} }{  g_{ne} }
!|       \hat{A}_{\parallel}
!| \qquad{\rm where}\qquad  \hat{A}_\parallel \equiv e A_\parallel / T_e .
!| \end{equation}
!|
!| The perturbed free electron density then follows by neglecting electron
!| inertia in Eq. (2)
!| \begin{equation}
!| {{\tilde{n}_{ef}}\over n_e}\equiv \hat{n}_{ef} =
!|    \hat{\phi} - ( \omega - \omega_{*e} ) \hat{A}_\parallel
!|     / ( k_\parallel c ).
!| \end{equation}
!|
!| In order to relate $\phi$ and $A_{\parallel}$, use the
!| electron continuity equation for the free electrons,
!| \begin{equation}
!| \frac{\partial n_{ef}}{\partial t}  +
!|   {\mathbf v}_{E}\cdot{\mathbf \nabla} n_{ef}
!|   + n_{ef}{\mathbf \nabla}\cdot{\mathbf v}_{E}
!|    +  {\mathbf \nabla}\cdot \left( n_{ef} {\mathbf v}_{*e} \right)
!|   - \frac{1}{e}{\mathbf \nabla}\cdot
!|    \left( J_{\parallel ef}\hat{{\mathbf B}} \right)
!|    = 0
!| \end{equation}
!| Only the free electrons contribute to the parallel electron current.
!| When parallel ion motion is included in the current,
!| Ampere's law implies,
!| \begin{equation}
!| J_{\parallel} = J_{\parallel ef} + J_{\parallel i} = - \frac{c}{4 \pi}
!|    \Delta A_{\parallel}
!| \end{equation}
!| Combining the last four equations,
!| eliminating $\tilde{n}_{ef}$, using the curvature relations
!| ${\mathbf \nabla} \cdot {\mathbf v}_E$ and
!|  ${\mathbf \nabla} \cdot \left( n_e {\mathbf v}_{*e} \right)$,
!| and assuming the same unperturbed temperature for the free and trapped
!| electrons, we obtain
!| \begin{equation}
!| \hat{A}_{\parallel}
!|  = \frac{k_{\parallel} c \left( \omega_{*e} - \omega \right) \hat{\phi}}
!| { \omega \left( \omega_{*e} - \omega \right)
!|    + \left( k_{\perp} \rho k_{\parallel} v_{A} \right)^2
!|    + \omega_{De} \left( \omega - \omega_{*eT} \right) }
!| \end{equation}
!| where $\omega_{*eT} = \omega_{De} ( g_{Te} + g_{ne} ) / 2 $.
!|
!| Justification for the use of all of the above equations
!| and more complete derivations
!| can be found in the references.\cite{weil92a,nord90a}
!|
!| Note that
!| \[ n_i \nabla \cdot \vec{v}_i = - \Partial{n_i}{t} - \vec{v}_i \cdot \vec{n}_i \]
!| from the equation of continuity can be used to eliminate
!| $ \nabla \cdot \vec{v}_i $ in the energy equation (divided by 3/2):
!| \[ n_i \left( \Partial{}{t} + \vec{v}_i \cdot \nabla \right) T_i
!| -\frac{2}{3} T_i \Partial{n_i}{t}
!|  - \frac{2}{3} T_i \vec{v}_i \cdot \nabla n_i
!|  = \frac{5}{3} n_i ( \vec{v}_{*i} - \vec{v}_{Di} ) \cdot \nabla T_i \]
!| and a similar equation for the trapped electrons.
!| The $\vec{v}_{*i}$ (and corresponding $\vec{v}_{*e}$)
!| contributions drop out because
!| \[ n_i \vec{v}_{*i} \cdot \nabla T_i
!|  - \frac{2}{3} T_i \vec{v}_{*i} \cdot \nabla n_i
!|  - \frac{5}{3} n_i \vec{v}_{*i} \cdot \nabla T_i
!|  = \frac{2}{3} \vec{v}_{*i} \cdot \nabla ( n_i T_i ) = 0 \]
!| since the diamagnetic drift
!| $ \vec{v}_{*i} = \hat{B} \times \nabla ( n_i T_i ) / ZeB n_i $
!| is perpendicular to the gradient of the pressure
!| $ \nabla ( n_i T_i ) $.
!|
!| See the derivations in the notes by Weiland\cite{weil92a} for
!| the relations
!| \[ \nabla \cdot \delta ( {n_i v_{*i}} )
!|  = \vec{v}_{Di} \cdot \nabla \delta ( {n_i T_i} ) / T_i \]
!| (Ref. [2], page 132 Eq. (4.20)),
!| \[ \nabla \cdot \tilde{v}_E
!|    = \frac{Z e}{T_i} \vec{v}_{Di} \cdot \nabla \tilde{\phi} \]
!| (Ref. [2], page 132 Eq. (4.21)), and
!| \[ \nabla \cdot [ n_i ( \vec{v}_{pi} + \vec{v}_{\pi i} ) ]
!|  \approx - i n_i k_y^2 \rho_{si}^2
!|  [ \omega - \omega_{*i} ( 1 + \eta_i ) ] \frac{e \tilde{\phi}}{T_e} \]
!| (Ref. [2], page 25 Eq. (1.28)) where
!| \[ \rho_{si}^2 = \frac{T_e}{m_i \Omega_i^2} \]
!| with $ v_{thi}^2 = 2 T_i / m_i $ and $ \Omega_i = Z e B / m_i $.
!|
!| Using these relations in the continuity, momentum, and energy equations,
!| we obtain the following perturbed relations,
!| \[ (-\omega + \omega_{Di}) \hat{n}_i + \omega_{Di} \hat{T}_i
!|  + [ ( \omega_{Di} - \omega_{*i} ) Z_i T_e/T_i
!|  - k_y^2 \rho_{si}^2 ( \omega - \omega_{*i} ( 1 + \eta_i ) )
!|  ] \hat{\phi} + k_\parallel v_{\parallel i} = 0 \]
!| \[ - \omega m_i v_{\parallel i} + k_\parallel Z_i T_e \hat{\phi}
!|  + k_\parallel T_i ( \hat{n}_i + \hat{T}_i ) = 0 \]
!| \[ (-\omega + \frac{5}{3} \omega_{Di} ) \hat{T}_i
!|    + \frac{2}{3} \omega \hat{n}_i
!|    + \omega_{*e} \frac{g_{ni}}{g_{ne}}
!|    ( \eta_i - \frac{2}{3} ) \hat{\phi} = 0 \]
!| with corresponding equations for trapped electrons ($i=et$)
!| and impurities ($i=Z$).
!| Here $\hat{n} \equiv \tilde{n} / n$,
!| $\hat{T} \equiv \tilde{T} / T$, and
!| $\hat{\phi} \equiv e \tilde{\phi} / T_e$ are dimensionless forms
!| of the perturbation.
!|
!| Normalizing all frequencies by
!| $\omega_{De} = \vec{k} \cdot \vec{v}_{De} = 2 k_y T_e /  e B R $,
!| and normalizing the parallel ion velocity by the speed of sound
!| in that ion,
!| $ \hat{v}_{\parallel i} \equiv v_{\parallel i} / c_{si} $
!| where $ c_{si} \equiv \sqrt{T_e / m_i} $,
!| we obtain
!| \[ ( \hat{\omega} + \frac{T_i}{Z_i T_e} ) \hat{n}_i
!|   + \frac{T_i}{Z_i T_e} \hat{T}_i
!|   + \left[ 1 - \frac{g_{ni}}{2}
!|   + k_y^2 \rho_{si}^2 \left( \hat{\omega}
!|     + \frac{T_i}{Z_i T_e} \frac{g_{ni}+g_{Ti}}{2}
!|     \right) \right] \hat{\phi}
!|   - \frac{k_\parallel c_{si}}{\omega_{De}} \hat{v}_{\parallel i}
!|   = 0 \]
!| \[ \hat{\omega} \hat{v}_{\parallel i}
!|  - \frac{k_\parallel c_{si}}{\omega_{De}} Z_i \left[ \hat{\phi}
!|    + \frac{T_i}{Z_i T_e} ( \hat{n}_i + \hat{T}_i ) \right] = 0 \]
!| \[ \left( - \hat{\omega} - \frac{5}{3} \frac{T_i}{Z_i T_e}
!|    \right) \hat{T}_i + \frac{2}{3} \hat{\omega} \hat{n}_i
!|    + \frac{1}{2} \left( g_{Ti} - \frac{2}{3} g_{ni}
!|      \right) \hat{\phi} = 0 \]
!| where $ \hat{\omega} \equiv \omega / \omega_{De} $
!| and
!| \[ g_{ni} \equiv - R \hat{x} \cdot \nabla n_i / n_i \]
!| is the normalized gradient.
!|
!| Now use
!| \[ f_t \hat{n}_{et} = ( 1 - Z f_Z - f_s ) \hat{n}_H
!|    + Z f_Z \hat{n}_Z - ( 1 - f_t ) \hat{\phi} \]
!| to eliminate the perturbed trapped electron density $\hat{n}_{et}$
!| in favor of the perturbed ion densities and perturbed potential.
!| The perturbed hydrogen and impurity equations remain as above
!| while the trapped electron density and energy equations become:
!| \[ (1-\hat{\omega}) ( 1 - Z f_Z - f_s ) \hat{n}_H
!|  + (1-\hat{\omega}) Z f_Z \hat{n}_Z  + f_t \hat{T}_{et}
!|  - [ 1 - f_t g_{ne} / 2 - ( 1 - f_t ) \hat{\omega} ] \hat{\phi}
!|  = 0 \]
!| \[ \frac{2}{3} \hat{\omega} ( 1 - Z f_Z - f_s ) \hat{n}_H
!|  + \frac{2}{3} \hat{\omega} Z f_Z \hat{n}_H
!|  + \left( \frac{5}{3} - \hat{\omega} \right) f_t \hat{T}_{et}
!|  + \left[ f_t \frac{g_{Te} - (2/3) g_{ne}}{2}
!|    - \frac{2}{3} \hat{\omega} ( 1 - f_t ) \right] \hat{\phi}
!|  = 0  \]
!|
!| The parallel wavenumber $ k_\parallel = \sqrt{\alpha} / ( q R ) $
!| is estimated in reference \cite{weil92a}
!| in terms of the ballooning parameter $\alpha$ given by
!| \begin{equation}
!| \alpha = \left( k_1 + \sqrt{k_1 + S^2 k_2} \right) / 2
!| \end{equation}
!| with
!| \begin{equation}
!| k_1 = \frac{q^2}{4} \left( k_{\perp}\rho \right)^2
!|    \sqrt{ \frac{T_i}{T_e} \frac{ ( g_{Ti} + g_{ni} ) }{ 2 ( 1 - f_t ) } }
!| \qquad{\rm and}\qquad
!| k_2 = q^2 \left( k_{\perp}\rho \right)^4
!|     \frac{T_i}{T_e} \frac{ ( g_{Ti} + g_{ni} ) }{ 2 ( 1 - f_t ) }
!| \end{equation}
!| where $ S = d \ln q / d \ln r $ is the magnetic shear.
!|
!| In all of these expressions we use the notation
!| $$ \omega_{*e} = k_y T_e g_{ne} / ( e B R ) $$
!| $$ \omega_{*i} = - k_y T_i g_{ni} / ( Z_i e B R )
!|     = - T_i \omega_{*e} / ( Z_i T_e ) $$
!| $$ \omega_{De} = 2 k_y T_e / e B R = 2 \omega_{*e} / ( g_{ne} ) $$
!| $$ \omega_{Di} = 2 k_y T_i / Z_i e B R = 2 \omega_{*i} / ( g_{ni} )
!|    = - 2 T_i \omega_{*e} / (  g_{ni} T_e ) $$
!| $$ g_{ni} = - R \hat{x} \cdot \nabla n_i / n_i $$
!| $$ \eta_i = g_{Ti} / g_{ni}  $$
!| $$ \eta_e = g_{Te} / g_{ne} $$
!| and
!| $$ f_t \approx \sqrt{ \frac{ 2 r/R }{ 1 + r/R } } $$
!| is the trapped electron fraction.
!| Note that here I define
!| \[ \epsilon_{ni} \equiv 1 / ( g_{ni} ) \]
!| while $ \epsilon_{ni} $ in the Weiland papers is defined to be
!| twice this value.
!|
!| The basic technique is to set up the generalized eigenvalue problem
!| \[ A v = \lambda B v \]
!| where $ \lambda = \hat{\omega} + i \hat{\gamma} $ is the eigenvalue
!| and $ v $ is the corresponding eigenvector.
!| Hence the eigenvalues give the frequency and growth rates of the
!| modes while the eigenvectors give the phase of the perturbed
!| variables relative to one another.
!|
!| In order to approximate the reduction of transport with $E\times B$
!| velocity shear the $E\times B $ shearing rate is subtracted from the
!| linear growth rates obtained as above. Currently this is implemented as
!| follows: The generalized eigenvalue problem is redefined using
!|  \[ (A - {\it i} \omega_{E\times B} B) v' = {\lambda}' B v' \]
!| and fluxes are calculated rom the new growth rates   ${\lambda}'$
!| and eigenvectors $v'$. This method do give the same results as a
!| direct reduction of the growth rates only would give but is more
!| integrated with the current framework of the model.
!|
!|
!| In this routine, the perturbed variables are always computed in
!| the following order:
!| \[ \hat{\phi}, \hat{T}_H, \hat{n}_H, \hat{T}_{et},
!|     \hat{n}_Z, \hat{T}_Z, \ldots \]
!|
!
!.. the table of argument input and output variable names and
!   definitions is givin in tex documentation
!
!  Note that external subroutines need to be provided:
!
!  tomsqz   Code wrapper for ACM/TOMS routine 535 implementing the QZ algorithm
!           complex generailized eigenvalue problem (eigenvalues and
!           eigenvectors)
!           The algorithm itself consists of the three routines
!           cqzhes, cqzval and cqzvec
!
!
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer idp
      parameter ( idp = 15 )
!
      integer                                                            &  
     &   letain,   lprint,  neq,     ndim,   nmodes, nerr,      nout     &  
     & , ieq,      idim,    imatrx,  j1,     j2,     j
!
! ieq    = number of equations
!
! imatrx = the number of elements computed along each row and
!          column of the transport matrix
!
      logical inital
      data inital /.true./
!
      dimension                                                          &  
     &   letain(32),        cetain(32),      omega(*)                    &  
     & ,  gamma(*),         vfluxout(*),     chieff(*)
!
      REAL*8                                                             &  
     &   cetain,   gnein,   gnhin,   gnzin,  gtein,  gthin,     gtzin    &  
     & , tauhin,   tauzin,  fnzin,   czin,   azin,   fnsin,     betaein  &  
     & , ftrapein, vef,     q,       shear,  h,      ekyrhoin,  wexb     &  
     & , omega,    gamma,   chieff,  vfluxout
!
!  Note:  Normally the transport matrix is  ieq-1 by ieq-1
!    however, if there are 6 equations, compute only a 4 by 4 matrix
!    since the impurity temperature equation is not used by most
!    transport codes including BALDUR
!
      REAL*8                                                             &  
     &   zamr(idp,idp),       zami(idp,idp),     zbmr(idp,idp)           &  
     & , zbmi(idp,idp),       zamrt(idp,idp),    zamit(idp,idp)          &  
     & , zbmrt(idp,idp),      zbmit(idp,idp),    zvr(idp,idp)            &  
     & , zvi(idp,idp),        zomega(idp),       zgamma(idp)             &  
     & , zalfr(idp),          zalfi(idp),        zbeta(idp)
      REAL*8 AREAL
!
      integer ifail
!
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
!
! Here, zbeta(j) will be 0.0 in the case of an infinite eigenvalue
!
      COMPLEX*16  zevec(idp,idp)
!
      REAL*8                                                             &  
     &   zepsmach, zepsqrt,   zone,  ztwo,   zthree,   zfour,   zfive    &  
     & , zhalf,    zquarter,  ztvr,  zftr,   ztwohlf,  zgne,    zgnh     &  
     & , zgnz,     zgns,      zgte,  zgth,   zgtz,     ztauh,   ztauz    &  
     & , zft,      zimp,      zfnz,  zmass,  zfns,     zflh,    zflz     &  
     & , zgamax,   zetae,     zetah
!
! zepsmach = machine epsilon
! zepsqrt  = sqrt ( machine epsilon )
! zone     = 1.0
! ztwo     = 2.0
! zthree   = 3.0
! zfour    = 4.0
! zfive    = 5.0
! zhalf    = 0.5
! zquarter = 0.25
! ztwohlf  = 2.5
! ztvr     = 2. / 3.
! zftr     = 5. / 3.
!
      REAL*8                                                             &  
     &   zflxph,   zflxnh,   zflxpe, zflxnz, zflxpz,   zphsph,  zphsnh   &  
     & , zphspe,   zphsnz,   zphspz, ztemp1, zreal,    zimag             &  
     & , ztempa(idp),        ztempb(idp),    zerreal,  zerimag, zerrmax
!
!  These are the thermal and particle fluxes and effective diffusivities
!
      REAL*8 zflxm(idp), zchim(idp)
!
!  Normalized transport fluxes         Eff. Diffusivities
!       zflxm(1)   n_H T_H             zchim(1)    chi_i
!       zflxm(2)   n_H                 zchim(2)    D_i
!       zflxm(3)   n_e T_e             zchim(3)    chi_e
!       zflxm(4)   n_Z                 zchim(4)    D_q
!       zflxm(5)   n_Z T_Z             zchim(5)    chi_q
!
      REAL*8                                                             &  
     &   bt,        bt1,      zeni,  k1,      k2,      zkpsh,   zkpsz    &  
     & , zalp,      zalf,     zanorm,         zbetae
!
!|
!| Definitions of some of the internal variables:
!|
!| \renewcommand{\arraystretch}{1.4}
!| \begin{center}
!| \begin{tabular}{lp{4.0in}}
!| variable & meaning \\
!| \\
!| {\tt zgth} & $ - \frac{R}{T_H} \Partial{T_H}{r} $ \\
!| {\tt zgte} & $ - \frac{R}{T_e} \Partial{T_e}{r} $ \\
!| {\tt zgtz} & $ - \frac{R}{T_Z} \Partial{T_Z}{r} $ \\
!| {\tt zgnh} & $ - \frac{R}{n_H} \Partial{n_H}{r} $ \\
!| {\tt zgnz} & $ - \frac{R}{n_Z} \Partial{n_Z}{r} $ \\
!| \end{tabular}  \end{center}
!|
!
      save                                                               &  
     &   idim,  zepsmach, zepsqrt,  zone,   ztwo,   zthree,   zfour      &  
     & , zfive, zhalf,    zquarter, ztvr,   zftr,   ztwohlf, inital
!
!..initialize variables
!
!
      ieq = max ( 2, neq )
!
      vef = cetain(15) * vef
      imatrx = min ( ieq - 1, 4 )
      if ( letain(2) .gt. 0 ) imatrx = min ( imatrx, letain(2)-1 )
!
      if ( inital ) then
!
        idim = idp
!
        zone   = 1.0_R8
        ztwo   = 2.0_R8
        zthree = 3.0_R8
        zfour  = 4.0_R8
        zfive  = 5.0_R8
!
        zhalf    = zone  / ztwo
        ztwohlf  = ztwo  + zhalf
        zquarter = zone  / zfour
        ztvr     = ztwo  / zthree
        zftr     = zfive / zthree
!
        zepsmach = zhalf
  2     if ( zhalf * zepsmach + zone .gt. zone ) then
          zepsmach = zhalf * zepsmach
          go to 2
        endif
!
        zepsqrt = sqrt ( zepsmach )
!
        inital = .false.
!
      endif
!
!..end of initialization
!
      nerr = 0
!
!..print header
!
      if ( lprint .gt. 0 ) then
!
        write (nout,*)
        write (nout,*)                                                   &  
     & 'Weiland-Nordman eigenvalue equations, subroutine weiland14flux'
        write (nout,*) '(all frequencies normalized by omega_{De})'
        write (nout,*) '(all diffusivities normalized by '               &  
     &    ,'omega_{De} / k_y^2'
!
        write (nout,108) 'gnein', gnein
        write (nout,108) 'gnhin', gnhin
        write (nout,108) 'gnzin', gnzin
        write (nout,108) 'gtein', gtein
        write (nout,108) 'gthin', gthin
        write (nout,108) 'gtzin', gtzin
        write (nout,108) 'tauhin', tauhin
        write (nout,108) 'tauzin', tauzin
        write (nout,108) 'ftrapein', ftrapein
        write (nout,108) 'fnzin', fnzin
        write (nout,108) 'czin', czin
        write (nout,108) 'azin', azin
        write (nout,108) 'fnsin', fnsin
        write (nout,108) 'betaein', betaein
        write (nout,108) 'vefin', vef
        write (nout,108) 'qin', q
        write (nout,108) 'shearin', shear
        write (nout,108) 'ekyrhoin', ekyrhoin
        write (nout,108) 'wexb', wexb
 108    format (1x,a8,' = ',1pe14.6,',')
!
      endif
!
!...copy to local variables
!
      zgne = gnein
      zgnh = gnhin
      zgnz = gnzin
      zgth = gthin
      zgte = gtein
      zgtz = gtzin
!
!..check validity of input data
!
      if ( ndim .gt. idim ) then
         nerr = -6
         return
      endif
!
!..initialize arrays
!
      do j1=1,ndim
        omega(j1)    = 0.0_R8
        gamma(j1)    = 0.0_R8
        vfluxout(j1) = 0.0_R8
        chieff(j1)   = 0.0_R8
      enddo
!
      do j1=1,idp
        zomega(j1)   = 0.0_R8
        zgamma(j1)   = 0.0_R8
        zalfr (j1)   = 0.0_R8
        zalfi (j1)   = 0.0_R8
        zbeta (j1)   = 0.0_R8
        zchim(j1)    = 0.0_R8
        zflxm(j1)    = 0.0_R8
        do j2=1,idp
          zevec(j1,j2)  = ( 0.0_R8, 0.0_R8)
          zamr (j1,j2)  = 0.0_R8
          zami (j1,j2)  = 0.0_R8
          zbmr (j1,j2)  = 0.0_R8
          zbmi (j1,j2)  = 0.0_R8
          zamrt(j1,j2)  = 0.0_R8
          zamit(j1,j2)  = 0.0_R8
          zbmrt(j1,j2)  = 0.0_R8
          zbmit(j1,j2)  = 0.0_R8
          zvr  (j1,j2)  = 0.0_R8
          zvi  (j1,j2)  = 0.0_R8
        enddo
      enddo
!
!..set up initial gradients
!
      zimp   = czin
      zmass  = azin
      zfnz   = fnzin * zimp
      zfns   = fnsin
!
!..Save zgns = - R ( d n_s / d r ) / n_e
!    where n_s is the fast ion density
!
      zgns = zgne - zgnh * ( zone - zfnz - zfns ) - zgnz * zfnz
!
!..compute the rest of the dimensionless variables needed
!
!     ztauz = T_Z / ( Z T_e ),     zfnz  = Z n_Z / n_e
!     zimp  = Z,                   zflh  = k_y^2 \rho_{sH}^2
!     zmass = m_Z / m_H,           zflz  = k_y^2 \rho_{sZ}^2
!
 
      ztauh  = tauhin
      zeni = zhalf * zgne
      zbetae = max ( cetain(20) * betaein, zepsqrt )
      zimp   = max ( czin, zone )
      zmass  = max ( azin, zone )
      zfnz   = fnzin * zimp
      zft    = ftrapein
      zflh   = ekyrhoin**2
!
      if ( neq .gt. 4 ) then
        ztauz  = tauzin / czin
        zflz   = zmass * zflh / zimp**2
      else
        ztauz  = zone
        zflz   = 0.0_R8
      endif
      zetae = zgte / zgne
      zetah = zgth / zgne
!
!..diagnostic output
!
      if ( lprint .gt. 2 ) then
        write (nout,*)
        write (nout,*) '--------------------------------------'
        write (nout,*)
        write (nout,*)
        write (nout,*) zgnh,' = zgnh'
        write (nout,*) zgne,' = zgne'
        write (nout,*) zgnz,' = zgnz'
        write (nout,*) zgns,' = zgns'
        write (nout,*) zgth,' = zgth'
        write (nout,*) zgte,' = zgte'
        write (nout,*) zgtz,' = zgtz'
        write (nout,*) zetae,' = zetae'
        write (nout,*) zetah,' = zetah'
        write (nout,*)
!
        write (nout,*) ztauh,' = ztauh'
        write (nout,*) ztauz,' = ztauz'
        write (nout,*)
        write (nout,*) zft,' = zft'
        write (nout,*) zimp,' = zimp'
        write (nout,*) zmass,' = zmass'
        write (nout,*) fnzin,' = fnz'
        write (nout,*) zfnz,' = zfnz'
        write (nout,*) zfns,' = zfns'
        write (nout,*) zflh,' = zflh'
        write (nout,*) zflz,' = zflz'
        write (nout,*)
        write (nout,*) zepsqrt,' = zepsqrt'
        write (nout,*) zepsmach,' = zepsmach'
        write (nout,109) 'gnsin', zgns
 109  format (1x,a8,' = ',1pe14.6,', computed from quasi-neutrality')
 
      endif
 
      if ( ieq .gt. 4 ) then
!
!..constants from Nilsson and Weiland, NF 34 (1994) 803  section 2
!  and from Weiland and Hirose NF 32 (1992) 151
!
!  ** Note:  zanorm = 1./kpc is used to normalize A_\parallel and K
!
        zanorm = cetain(11)
!
        bt  = 1.5_R8
        bt1 = bt - ztwohlf
!
!  Expressions from Weiland and Hirose, NF 32 (1992) 151  Eq. (6)
!  Note that k1 rather than k1**2 appears under sqrt in zalp.
!  This is the result of a numerical fit rather than analytic solution.
!  Also, only free electrons are used from zbetae in zalf
!  Hence, zbetae -> zbetae * (zone-zft).
!
        if ( zgnh + zgth .lt. zepsqrt ) then
!
          zalp  = 0.0_R8
          zalf  = 0.0_R8
          zkpsh = 0.0_R8
          zkpsz = 0.0_R8
!
        else
!
          k1 = zquarter * q * q * zflh *                                 &  
     &         sqrt( zhalf * (zgnh + zgth) * ztauh / (zone-zft) )
          k2 = q * q * zflh * zflh * zhalf * (zgnh + zgth) *             &  
     &         ztauh / ( zone-zft)
 
          zalp = zhalf * ( k1 + sqrt( k1 + shear * shear * k2) )
          zalf = zalp / (ztwo * zflh * q * q * zbetae * (zone-zft))
 
          zkpsh = cetain(10) * zhalf * sqrt( zalp / zflh ) / q
 
          if ( zmass .gt. zepsqrt ) then
            zkpsz = zkpsh / sqrt ( zmass )
          else
            zkpsz = zkpsh
          endif
!
        endif
!
        if (lprint .gt. 5 ) then
          write (nout,*)
          write (nout,*) k1,'  = k1'
          write (nout,*) k2,'  = k2'
          write (nout,*) zalp,' = zalp'
          write (nout,*) zalf,' = zalf'
          write (nout,*) zkpsh,' = zkpsh'
          write (nout,*) zkpsz,' = zkpsz'
          write (nout,*)
          write (nout,*) bt1,' = bt1'
          write (nout,*) zanorm,' = zanorm'
        endif
!
      endif
!
!..NOTE! Below only the section, appropriate for the number of equations
!        in the version of the Weiland model that is chosen, is used.
!        In each section the A and the B matrices are constructed.
!
!  For the MMM95 transport model used by the Lehigh Group, ieq = 10, and
!  collisions are turned off.
!
      if (ieq .eq. 11) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..Eleven equations with impurities, trapped electrons, FLR,
!  collisional effects, parallel hydrogenic and impurity ion motion,
!  electromagnetic (finite beta) effects
!
!  equations for
!  e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vph, Av, K, Vpz
!
      if (lprint .gt. 5 ) then
      write(nout,*)                                                      &  
     &  'Eleven eqns for '                                               &  
     &  ,'e phi / T_e, T_H, n_H, T_et, n_Z, T_Z, F, Vph, Av, K, Vpz'
      endif
!
!
!  hydrogenic density
!
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))    
      zamr(1,2) = - ztauh
      zamr(1,3) = - ztauh
      zamr(1,8) = zkpsh
!
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
!
!  hydrogenic energy
!
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
!
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
!
!  trapped electron density expressed through quasineutrality
!
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,9) =  - ( zone - zft ) * zeni * zanorm
      zami(3,9) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,10) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,10) = - vef * ( zone - zft ) * zanorm
!
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,10) = ( zone - zft ) * zanorm
!
!  trapped electron energy
!
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - 2.5_R8* (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
!
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
!
!      impurity density
!
      zamr(5,1) = - zone + zhalf*(zgnz - zflz * ztauz * (zgnz + zgtz))
      zamr(5,5) = - ztauz
      zamr(5,6) = - ztauz
      zamr(5,11) = zkpsz
!
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
!
!      impurity energy
!
      zamr(6,1) = zhalf * (zgtz - ztvr * zgnz)
      zamr(6,6) = - ztauz * zftr
!
      zbmr(6,5) = - ztvr
      zbmr(6,6) = zone
!
!  variable F
!
      zamr(7,1) = zhalf * zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = - vef
!
      zbmr(7,1) = - zone
      zbmr(7,7) = zone
!
!  Parallel hydrogenic ion motion Vpi/Cs
!
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh * ztauh
      zamr(8,3) = zkpsh * ztauh
!
      zbmr(8,8) = zone
!
!  electromagnetic parallel vector potential Av = e A_par /Te
!
      zamr(9,1) = zeni
      zamr(9,8) = zkpsh * ( zone - zfnz - zfns ) / ( zone - zft )
      zamr(9,9) = ( zhalf * (zgne + zgte) - zalf * zflh ) * zanorm
      zamr(9,10) = - zeni * zanorm
      zamr(9,11) = zkpsz * zfnz / ( zone - zft )
!
      zbmr(9,1) = zone
      zbmr(9,9) = zanorm
      zbmr(9,10) = - zanorm
!
!  time derivative of Av
!
      zamr(10,10) = zone
!
      zbmr(10,9) = zone
!
!..Impurity parallel ion motion
!
      zamr(11,1) = zimp * zkpsz
      zamr(11,5) = zimp * zkpsz * ztauz
      zamr(11,6) = zimp * zkpsz * ztauz
!
      zbmr(11,11) = zone
!
!
      elseif (ieq .eq. 10) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..Ten equations with impurities, trapped electrons, FLR,
!  collisional effects, parallel ion motion, and electromagnetic effects
!
!  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vp, Av, K
!
      if (lprint .gt. 5 ) then
      write(nout,*)                                                      &  
     &'Ten eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Av, K'
      endif
!
!
!  hydrogenic density
!
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))    
      zamr(1,2) = - ztauh
      zamr(1,3) = - ztauh
      zamr(1,8) = zkpsh
!
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
!
!  hydrogenic energy
!
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
!
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
!
!  trapped electron density expressed through quasineutrality
!
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,9) =  - ( zone - zft ) * zeni * zanorm
      zami(3,9) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,10) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,10) = - vef * ( zone - zft ) * zanorm
!
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,10) = ( zone - zft ) * zanorm
!
!  trapped electron energy
!
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - 2.5_R8* (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
!
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
!
!      impurity density
!
      zamr(5,1) = - zone + zhalf*(zgnz - zflz * ztauz * (zgnz + zgtz))
      zamr(5,5) = - ztauz
      zamr(5,6) = - ztauz
!
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
!
!      impurity energy
!
      zamr(6,1) = zhalf * (zgtz - ztvr * zgnz)
      zamr(6,6) = - ztauz * zftr
!
      zbmr(6,5) = - ztvr
      zbmr(6,6) = zone
!
!  variable F
!
      zamr(7,1) = zhalf * zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = - vef
!
      zbmr(7,1) = - zone
      zbmr(7,7) = zone
!
!  Parallel ion motion Vpi/Cs
!
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh * ztauh
      zamr(8,3) = zkpsh * ztauh
!
      zbmr(8,8) = zone
!
!  electromagnetic parallel vector potential Av = e A_par /Te
!
      zamr(9,1) = zeni
      zamr(9,8) = zkpsh * ( zone - zfnz - zfns ) / ( zone - zft )
      zamr(9,9) = ( zhalf * (zgne + zgte) - zalf * zflh ) * zanorm
      zamr(9,10) = - zeni * zanorm
!
      zbmr(9,1) = zone
      zbmr(9,9) = zanorm
      zbmr(9,10) = - zanorm
!
!  time derivative of Av
!
      zamr(10,10) = zone
!
      zbmr(10,9) = zone
!
!
      elseif (ieq .eq. 9) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..Nine equations with impurities, trapped electrons, FLR,
!  collisional effects, and electromagnetic effects
!  and parallel ion motion in strong ballooning approximation
!
!  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Av, K
!
      if (lprint .gt. 5 ) then
      write(nout, *)                                                     &  
     &'Nine eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Av, K'
      endif
!
!
      H = cetain(12) * 0.5_R8* abs( shear ) / max ( q, zepsqrt )
!
!     ion continuity
!
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))    
      zami(1,1) = - H
      zamr(1,2) = - ztauh
      zami(1,2) = - ztauh*H
      zamr(1,3) = - ztauh
      zami(1,3) = - ztauh*H
!
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
!
!      ion energy
!
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
!
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
!
! total electron density expressed through quasineutrality
!
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,8) = - ( zone - zft ) * zeni * zanorm
      zami(3,8) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,9) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,9) = - vef * ( zone - zft ) * zanorm
!
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,9) = ( zone - zft ) * zanorm
!
!      trapped electron energy
!
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - 2.5_R8* (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
!
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
!
!      impurity density
!
      zamr(5,1) = -zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
      zamr(5,5) = -ztauz
      zamr(5,6) = -ztauz
!
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
!
!      impurity energy
!
      zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
      zamr(6,6) = -ztauz*zftr
!
      zbmr(6,5) = -ztvr
      zbmr(6,6) = zone
!
!  variable F
!
      zamr(7,1) = zhalf*zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = -vef
!
      zbmr(7,1) = -zone
      zbmr(7,7) = zone
!
!  electromagnetic parallel vector potential Av = e A_par /Te
!
      zamr(8,1) = zeni
      zamr(8,8) = ( zhalf*(zgne+zgte)-zalf*zflh ) * zanorm
      zamr(8,9) = - zeni * zanorm
!
      zbmr(8,1) = zone
      zbmr(8,8) = zanorm
      zbmr(8,9) = - zanorm
!
!   time derivative of Av
!
      zamr(9,9) = zone
!
      zbmr(9,8) = zone
!
!
      elseif (ieq .eq. 8) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..Eight equations with impurities, trapped electrons, FLR,
!  collisional effects, and parallel ion motion
!
!  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vp
!
      if (lprint .gt. 5 ) then
      write(nout,*)                                                      &  
     &'Eight eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Vp'
      endif
!
!
!     ion continuity
!
      zamr(1,1) = -zone + zhalf*(zgnh-zflh*ztauh*(zgnh+zgth))
      zamr(1,2) = -ztauh
      zamr(1,3) = -ztauh
      zamr(1,8) = zkpsh
!
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
!
!      ion energy
!
      zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
      zamr(2,2) = -ztauh * zftr
!
      zbmr(2,2) = zone
      zbmr(2,3) = -ztvr
!
! total electron density expressed through quasineutrality
!
      zamr(3,1) = zft*zeni - zone
      zami(3,1) = vef*(zone-zft)
      zamr(3,3) = zone - zfnz- zfns
      zami(3,3) = -vef*(zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = -vef*zfnz
      zami(3,7) = vef*zft
!
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
!
!      trapped electron energy
!
      zamr(4,1) = zft*zhalf*(zgte-ztvr*zgne)
      zami(4,1) = vef*ztvr*(bt-2.5_R8*(zone-zft))
      zami(4,3) = -vef*ztvr*bt1*(zone - zfnz - zfns)
      zamr(4,4) = zft*zftr
      zami(4,5) = -vef*ztvr*bt1*zfnz
      zami(4,7) = -zftr*vef*zft
!
      zbmr(4,1) = (zone-zft)*ztvr
      zbmr(4,3) = - (zone - zfnz - zfns)*ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = -zfnz*ztvr
!
!      impurity density
!
      zamr(5,1) = -zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
      zamr(5,5) = -ztauz
      zamr(5,6) = -ztauz
!
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
!
!      impurity energy
!
      zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
      zamr(6,6) = -ztauz*zftr
!
      zbmr(6,5) = -ztvr
      zbmr(6,6) = zone
!
!  variable F
!
      zamr(7,1) = zhalf*zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = -vef
!
      zbmr(7,1) = -zone
      zbmr(7,7) = zone
!
!  Parallel ion motion Vpi/Cs
!
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh*ztauh
      zamr(8,3) = zkpsh*ztauh
!
      zbmr(8,8) = zone
!
!
      elseif ( ieq .eq. 7 ) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..Seven equations with impurities, trapped electrons, FLR,
!  and collisions
!
!  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
!  Here, F is defined as F = GM*e phi/T_e
!  where GM=1+etae/(epsn*(omega-1+i*vef))
!
      if ( lprint .gt. 5 ) then
      write (nout,*)
      write (nout,*)                                                     &  
     & 'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and F'
      endif
!
!
!  hydrogen density
!
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
!
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
!
!  hydrogen energy
!
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
!
        zbmr(2,2) = zone
        zbmr(2,3) = - ztvr
!
!  trapped electron density
!
        zamr(3,1) = - zone + zft * zeni
        zami(3,1) = vef*(zone-zft)
        zamr(3,3) = zone - zfnz - zfns
        zami(3,3) = -vef*(zone - zfnz - zfns)
        zamr(3,4) = zft
        zamr(3,5) = zfnz
        zami(3,5) = -vef*zfnz
        zami(3,7) = vef*zft
!
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
        zbmr(3,5) = zfnz
!
!  trapped electron energy
!
        zamr(4,1) = zft*zhalf*(zgte-ztvr*zgne)
        zami(4,1) = vef*ztvr*(bt-2.5_R8*(zone-zft))
        zami(4,3) = -vef*ztvr*bt1*(zone - zfnz - zfns)
        zamr(4,4) = zft * zftr
        zami(4,5) = -vef*ztvr*bt1*zfnz
        zami(4,7) = -zftr*vef*zft
!
        zbmr(4,1) = ( zone - zft ) *ztvr
        zbmr(4,3) = - ( zone - zfnz - zfns ) *ztvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * ztvr
!
!  impurity density
!
        zamr(5,1) = - zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
!
        zbmr(5,1) = zflz
        zbmr(5,5) = zone
!
!  impurity energy
!
        zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
        zamr(6,6) = - ztauz * zftr
!
        zbmr(6,5) = - ztvr
        zbmr(6,6) = zone
!
!  variable F
!
        zamr(7,1) = zhalf*zgte - zone
        zami(7,1) = vef
        zamr(7,7) = zone
        zami(7,7) = -vef
!
        zbmr(7,1) = -zone
        zbmr(7,7) = zone
!
!
      elseif ( ieq .eq. 6 ) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..Six equations with impurities, trapped electrons, and FLR
!
!  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, and T_Z
!
      if ( lprint .gt. 1 ) then
        write (nout,*)
        write (nout,*)                                                   &  
     &   'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
      endif
!
!
!  hydrogen density
!
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
!
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
!
!  hydrogen energy
!
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
!
        zbmr(2,2) = zone
        zbmr(2,3) = - ztvr
!
!  trapped electron density
!
        zamr(3,1) = zft*zeni - zone
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
        zamr(3,5) = zfnz
!
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
        zbmr(3,5) = zfnz
!
!  trapped electron energy
!
        zamr(4,1) = zft * zhalf * (zgte - ztvr*zgne)
        zamr(4,4) = zft * zftr
!
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - ( zone - zfnz - zfns ) * ztvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * ztvr
!
!  impurity density
!
        zamr(5,1) = -zone + zhalf * (zgnz - zflz*ztauz*(zgnz+zgtz))
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
!
        zbmr(5,1) = zflz
        zbmr(5,5) = zone
!
!  impurity energy
!
        zamr(6,1) = zhalf * (zgtz - ztvr*zgnz)
        zamr(6,6) = - ztauz * zftr
!
        zbmr(6,5) = - ztvr
        zbmr(6,6) = zone
!
!
      elseif ( ieq .eq. 5 ) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..5 equations with trapped electrons, FLR effects, and parallel ion motion
!
!  equations for e phi/T_e, T_H, n_i, T_e, and Vp
!
       if ( lprint .gt. 1 ) then
         write (nout,*)
         write (nout,*)                                                  &  
     &    ' Five eqns for e phi/T_e, T_H, n_H, T_e, and Vp'
       endif
!
!
!  ion continuity
!
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
        zamr(1,5) = zkpsh
!
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
!
!  ion energy
!
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
!
        zbmr(2,2) =   zone
        zbmr(2,3) = - ztvr
!
!  trapped electron continuity
!
!   Calculates the total electron density perturbation and replaces it
!   by the ion density perturbation.
!   The dilution factor 1-zfnz has now been added.
!
        zamr(3,1) = - zone + zft * zhalf * zgne
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
!
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
!
!  trapped electron energy
!
        zamr(4,1) = zft * ( zgte - ztvr * zgne ) * zhalf
        zamr(4,4) = zft * zftr
!
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - zone * ztvr
        zbmr(4,4) = zft
!
!
!  Parallel ion motion Vpi/Cs
!
        zamr(5,1) = zkpsh
        zamr(5,2) = zkpsh*ztauh
        zamr(5,3) = zkpsh*ztauh
!
        zbmr(5,5) = zone
!
      elseif ( ieq .eq. 4 ) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..4 equations with trapped electrons and FLR effects
!
!  equations for e phi/T_e, T_H, n_i, and T_e
!
       if ( lprint .gt. 1 ) then
         write (nout,*)
         write (nout,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
       endif
!
!
!  ion continuity
!
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
!
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
!
!  ion energy
!
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
!
        zbmr(2,2) =   zone
        zbmr(2,3) = - ztvr
!
!  trapped electron continuity
!
!   Calculates the total electron density perturbation and replaces it
!   The dilution factor 1-zfnz has now been added.
!
        zamr(3,1) = - zone + zft * zhalf * zgne
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
!
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
!
!  trapped electron energy
!
        zamr(4,1) = zft * ( zgte - ztvr * zgne ) * zhalf
        zamr(4,4) = zft * zftr
!
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - zone * ztvr
        zbmr(4,4) = zft
!
      elseif ( ieq .eq. 2 ) then
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!..two equations when trapped particles and FLR effects omitted
!
!
!  equations for e phi/T_e and T_H
!
       if ( lprint .gt. 1 ) then
         write (nout,*)
         write (nout,*) ' Two eqns for e phi/T_e and T_H'
       endif
!
        zamr(1,1) = zhalf * zgnh - ztauh - zone
        zamr(1,2) = - ztauh
        zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
        zamr(2,2) = - ztauh * zftr
!
        zbmr(1,1) = zone
        zbmr(1,2) = 0.0_R8
        zbmr(2,1) = - ztvr
        zbmr(2,2) = zone
!
      endif
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!..find the eigenvalues and eigenvectors using ACM/TOMS routine 535
!
      ifail = -1
!
!..Update A matrix w.r.t wexb and save local copies of the matrices
!
      do j1=1,ieq
       do j2=1,ieq
        zamr(j1,j2)  = zamr(j1,j2) + abs(wexb)*zbmi(j1,j2)
        zami(j1,j2)  = zami(j1,j2) - abs(wexb)*zbmr(j1,j2)
        zamrt(j1,j2) = zamr(j1,j2)
        zamit(j1,j2) = zami(j1,j2)
        zbmrt(j1,j2) = zbmr(j1,j2)
        zbmit(j1,j2) = zbmi(j1,j2)
       enddo
      enddo
 
!
!..diagnostic output for defining matrices
!
      if ( lprint .gt. 6 ) then
!
        write (nout,*)
        write (nout,*) ieq
!
        write (nout,*)
        write (nout,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,ieq
          write (nout,192) (zamr(j1,j2),j2=1,ieq)
        enddo
!
        write (nout,*)
        write (nout,*) ' zami(j1,j2)  j2 ->'
        do j1=1,ieq
          write (nout,192) (zami(j1,j2),j2=1,ieq)
        enddo
!
        write (nout,*)
        write (nout,*) ' zbmr(j1,j2)  j2->'
        do j1=1,ieq
          write (nout,192) (zbmr(j1,j2),j2=1,ieq)
        enddo
!
        write (nout,*)
        write (nout,*) ' zbmi(j1,j2)  j2->'
        do j1=1,ieq
          write (nout,192) (zbmi(j1,j2),j2=1,ieq)
        enddo
 192  format (1p10e12.4)
      endif
!
!... Solution of the generalized eigenvalue problem
!
      call tomsqz(idim,ieq,zamr,zami,zbmr,zbmi,ZALFR,ZALFI,ZBETA,        &  
     &            ZVR,ZVI,IFAIL)
!
!... Storing the results -  eigenvalues and eigenvectors
!
      zgamax = 0.0_R8
      do j=1,ieq
        ztemp1 = zbeta(j)
        if ( abs(zbeta(j)) .lt. zepsqrt ) ztemp1 = zepsqrt
        zomega(j) = zalfr(j)  / ztemp1
        zgamma(j) = zalfi(j)  / ztemp1
        omega(j) = zomega(j)
        gamma(j) = zgamma(j)
        zgamax = max ( zgamax, zgamma(j) )
        do j1=1,ieq
          zevec(j1,j) = CMPLX ( zvr(j1,j), zvi(j1,j),R8)
        end do
      enddo
 
      nerr = ifail
 
!
!... If no unstable roots are found set fluxes to zero and leave
!
      if ( zgamax .lt. zepsqrt ) go to 90
!
!..check the eigenfunctions
!
      if ( lprint .gt. 12 ) then
        write (nout,*)
        write (nout,*) ' Checking eigenfunctions'
      endif
!
!  Real and imaginary parts
!
      zerrmax = 0.0_R8
      nmodes = 0
!
! Loop over number of possible modes
!
      do j=1,ieq
!
! Start an error check on the eigensolution comparin RHS w LHS
!
        do j1=1,ieq
          ztempa(j1) = 0.0_R8
          ztempb(j1) = 0.0_R8
          do j2=1,ieq
            zerreal =                                                    &  
     &            zamrt(j1,j2) * REAL (zevec(j2,j))                     &  
     &          - zamit(j1,j2) * AIMAG(zevec(j2,j))                      &  
     &          - zbmrt(j1,j2) * REAL (zevec(j2,j)) * zomega(j)         &  
     &          + zbmrt(j1,j2) * AIMAG(zevec(j2,j)) * zgamma(j)          &  
     &          + zbmit(j1,j2) * REAL (zevec(j2,j)) * zgamma(j)         &  
     &          + zbmit(j1,j2) * AIMAG(zevec(j2,j)) * zomega(j)
 
            zerimag =                                                    &  
     &            zamrt(j1,j2) * AIMAG(zevec(j2,j))                      &  
     &          + zamit(j1,j2) * REAL (zevec(j2,j))                     &  
     &          - zbmrt(j1,j2) * AIMAG(zevec(j2,j)) * zomega(j)          &  
     &          - zbmrt(j1,j2) * REAL (zevec(j2,j)) * zgamma(j)         &  
     &          + zbmit(j1,j2) * AIMAG(zevec(j2,j)) * zgamma(j)          &  
     &          - zbmit(j1,j2) * REAL (zevec(j2,j)) * zomega(j)
 
            ztempa(j1) = ztempa(j1) + zerreal
            ztempb(j1) = ztempb(j1) + zerimag
 
          enddo
          zerrmax = max ( zerrmax, abs(ztempa(j1)), abs(ztempb(j1)) )
        enddo
 
        if ( lprint .gt. 12 ) then
          write (nout,*)
          write (nout,*) ' LHS - RHS for j =  ',j
          do j1=1,ieq
            write (nout,142) ztempa(j1), ztempb(j1)
          enddo
        endif
!
!|
!| The effective diffusivities can be computed directly from the
!| eigenvectors in the following way:
!| Consider, for example, the flux of hydrogen particles produced by
!| the perturbed $ E \times B $ motion of the plasma
!| \[ \Gamma_H = \tilde{n}_H \tilde{v}_E^* + c.c.
!|  = 2 ( {\rm Re} \tilde{n}_H {\rm Im} \tilde{\phi}
!|      - {\rm Im} \tilde{n}_H {\rm Re} \tilde{\phi} ) k_y / B. \]
!| Using the quasi-linear assumption, the saturation level is
!| approximated by
!| \[ \frac{e \tilde{\phi}}{T_e}
!|  \approx \frac{1}{k_x \rho_{sH}} \frac{ \gamma }{ k_y c_{sH} }
!|  = \frac{2}{R k_x} \frac{\gamma}{\omega_{De}}. \]
!| Hence, the hydrogen particle flux is given by
!| \[ \frac{ F_H }{ n_H } \frac{ R k_x^2 }{ \omega_{De} }
!|   = 2 \hat{\gamma}^2
!|   ( {\rm Re} \hat{n}_H {\rm Im} \hat{\phi}
!|   - {\rm Im} \hat{n}_H {\rm Re} \hat{\phi} ) / | \hat{\phi} |^2.
!|    \]
!| Correspondingly, the hydrogen heat (here, $ n_H T_H $) flux is given by
!| \[ \frac{ F_{p_H} }{ n_H T_H } \frac{ R k_x^2 }{ \omega_{De} }
!|   = 2 \hat{\gamma}^2
!|   ( {\rm Re} ( \hat{n}_H + \hat{T}_H ) {\rm Im} \hat{\phi}
!|   - {\rm Im} ( \hat{n}_H + \hat{T}_H ) {\rm Re} \hat{\phi} )
!|     / | \hat{\phi} |^2.
!|    \]
!|
!| The case of the electron heat flux (here, the flux of $ n_e T_e $)
!| is a little more complicated.  Note that
!|
!
!..compute effective diffusivities directly from eigenvalues
!  assume eigenvectors are arranged in the order of
!  e\phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z
!
!
!... Calculate norm of solution in terms of (e\phi/T_e) **2
!
        ztemp1 =  REAL(zevec(1,j)) *  REAL(zevec(1,j))                 &  
     &         + AIMAG(zevec(1,j)) * AIMAG(zevec(1,j))
!
!... Check if current mode j is unstable
!
        if ( zgamma(j).gt.zepsqrt .and. ztemp1.gt.zepsqrt ) then
 
           nmodes = nmodes + 1
!
! Define fluxes : Thermal hydrogen flux
!
          zreal =  REAL(zevec(2,j)) +  REAL(zevec(3,j))
          zimag = AIMAG(zevec(2,j)) + AIMAG(zevec(3,j))
!
!...phase difference
!
          zphsph = - ( zimag * REAL(zevec(1,j))                         &  
     &             -   zreal * AIMAG(zevec(1,j)) ) / ztemp1
!
!...flux from j:th mode
!
          zflxph = 2.0_R8* zphsph * zgamma(j) * zgamma(j)
 
!
!...flux summed over all unstable modes
!
          zflxm(1) = zflxm(1) + zflxph
!
! Define hydrogen density flux - phase shift
!
          zphsnh = - ( AIMAG(zevec(3,j)) * REAL(zevec(1,j))             &  
     &             - REAL(zevec(3,j)) * AIMAG(zevec(1,j)) ) / ztemp1
!
!...Flux from j:th mode
!
          zflxnh = 2.0_R8* zphsnh * zgamma(j) * zgamma(j)
!
!...Flux summed over all unstable modes
!
          zflxm(2) = zflxm(2) + zflxnh
!
! Define thermal electron flux
!
          if ( ieq .gt. 3 ) then
            zreal =  REAL(zevec(4,j))                                   &  
     &            + ( zone - zfnz - zfns ) * REAL(zevec(2,j))           &  
     &            + zfnz * REAL(zevec(5,j))
            zimag = AIMAG(zevec(4,j))                                    &  
     &            + ( zone - zfnz - zfns ) * AIMAG(zevec(2,j))           &  
     &            + zfnz * AIMAG(zevec(5,j))
!
!  Note, the electron heat flux is reduced by the fraction of
!  trapped electrons - phase shift
!
            zphspe = - zft * ( zimag * REAL(zevec(1,j))                 &  
     &                       - zreal * AIMAG(zevec(1,j)) ) / ztemp1
 
!
!... Flux from the j:th mode
!
            zflxpe = 2.0_R8* zphspe * zgamma(j)* zgamma(j)
!
!... Flux summed over all unstable modes
!
            zflxm(3) = zflxm(3) + zflxpe
!
          endif
!
!... Impurity density flux
!
          if ( ieq .gt. 4 ) then
!
!... phase difference
!
            zphsnz = - (AIMAG(zevec(5,j))*REAL (zevec(1,j))             &  
     &               -  REAL (zevec(5,j))*AIMAG(zevec(1,j)) ) / ztemp1
 
!
!... Flux from the j:th mode
!
            zflxnz =  2.0_R8* zphsnz * zgamma(j) * zgamma(j)
!
!... Flux summed over all unstable modes
!
            zflxm(4) = zflxm(4) + zflxnz
!
          endif
!
! Define the impurity heat flux
!
          if ( ieq .gt. 5 ) then
!
            zreal =  REAL(zevec(5,j)) +  REAL(zevec(6,j))
            zimag = AIMAG(zevec(5,j)) + AIMAG(zevec(6,j))
!
!... phase difference
!
            zphspz = - ( zimag * REAL (zevec(1,j))                      &  
     &               -   zreal * AIMAG(zevec(1,j)) ) / ztemp1
!
!... Flux from the j:th mode
!
            zflxpz = 2.0_R8* zphspz * zgamma(j) * zgamma(j)
!
!... Flux summed over all unstable modes
!
            zflxm(5) = zflxm(5) + zflxpz
!
          endif
!
!..header for diagnostic printout of frequencies, fluxes, and phases
!
        if ( letain(29) .gt. 0 .and. lprint .gt. 0) then
          write (nout,134)
!
 134  format (//' Diagnostic printout of frequencies, fluxes and         &  
     & phases'                                                           &  
     &  /' Note: fluxph = flux of hydrogenic thermal energy'             &  
     &  ,' = 2.0 * gamma^2 * phaseph'                                    &  
     &  /9x,' (phaseph is related to the phases of the perturbations)'   &  
     &  /7x,'fluxnh = flux of hydrogenic ions = 2.0 * gamma^2 *          &  
     & phasenh'                                                          &  
     &  /7x,'fluxpe = flux of electron thermal energy'                   &  
     &  ,' = 2.0 * gamma^2 * phasepe'                                    &  
     &  /7x,'fluxnz = flux of impurity ions = 2.0 * gamma^2 * phasenz'   &  
     &  /7x,'fluxpz = flux of impurity thermal energy'                   &  
     &  ,' = 2.0 * gamma^2 * phasepz'                                    &  
     &  //1x,'radius',t10,'omega',t20,'gamma'                            &  
     &  ,t30,'fluxph',t40,'phaseph',t50,'fluxnh',t60,'phasenh'           &  
     &  ,t70,'fluxpe',t80,'phasepe',t90,'fluxnz',t100,'phasenz'          &  
     &  ,t110,'fluxpz',t120,'phasepz  #m')
!
!..diagnostic printout of frequencies and fluxes mode by mode
!
        write (nout,135) cetain(29), zomega(j), zgamma(j)                &  
     &    , zflxph, zphsph, zflxnh, zphsnh, zflxpe, zphspe               &  
     &    , zflxnz, zphsnz, zflxpz, zphspz
        endif
!
 135  format (f7.3,1p12e10.2,' #m')
!
        endif
!
!...End of flux and diffusivity definitions loop
!
      enddo
 
!
!... Error diagnostics
!
      if ( lprint .gt. 0 ) then
!        if (abs(zerrmax) .gt. 100*zepsqrt) nerr = max (nerr , 1)
        if (abs(zerrmax) .gt. zepsqrt) nerr = max (nerr , 1)
        write (nout,*)
        write (nout,*) zerrmax,' = zerrmax'
        write (nout,*) nerr,                                             &  
     &       ' = nerr, error in eigenvalue in weiland14flux'
      endif
 
 142  format (1p10e12.4)
!
!..compute effective total diffusivities
!
      zchim(1) = zflxm(1) / sign( max( abs( zgth ), zepsqrt ),  zgth )
      zchim(2) = zflxm(2) / sign( max( abs( zgnh ), zepsqrt ),  zgnh )
      zchim(3) = zflxm(3) / sign( max( abs( zgte ), zepsqrt ),  zgte )
      zchim(4) = zflxm(4) / sign( max( abs( zgnz ), zepsqrt ),  zgnz )
      zchim(5) = zflxm(5) / sign( max( abs( zgtz ), zepsqrt ),  zgtz )
!
!..save effective diffusivities and fluxes
!
      do j1=1,min(ieq-1,5)
        chieff(j1)  = zchim(j1)
        vfluxout(j1) = zflxm(j1)
      enddo
 
!
      if ( lprint .gt. 2 ) then
!
!..print eigenvalues and eigenfunctions
!
        write (nout,121)
        do j=1,ieq
          write (nout,122) zomega(j), zgamma(j)
        enddo
 
 121    format (/' Solution of the eigenvalue equations'                 &  
     &   /t4,'zomega',t18,'zgamma')
 122    format (1p2e14.5,i5)
!
        write (nout,*)
        write (nout,*) ' Effective diffusivities'                        &  
     &    ,' normalized by omega_{De} / k_y^2'
!
        write (nout,130)
        write (nout,132) (zchim(j1),j1=1,ieq-1)
!
      endif
!
      if ( lprint .gt. 99 ) then
!
        write (nout,*)
        write (nout,*) ' Eigenvectors zevec(j1,j2) j2->'
        do j1=1,ieq
          write (nout,124) (zevec(j1,j2),j2=1,ieq)
 124      format (2(1p12e11.3))
        enddo
!
        write (nout,*)
      endif
!
      if ( ifail .gt. 0 ) then
        nerr = 5
        return
      endif
!
!
!
!
      if ( imatrx .lt. 1 ) go to 90
!
!..diagnostic printout
!
      if ( lprint .gt. 6 ) then
!
        write (nout,*)
        write (nout,*) ' vector zflxm(j1)'
        do j1=1,ieq-1
          write (nout,132) zflxm(j1)
        enddo
!
      endif
!
  90  continue
!
      if ( lprint .gt. 0 ) then
!
        write (nout,*)
        write (nout,*) '--------------------------------------'          &  
     &    ,'  Final results from sbrtn weiland14flux'
!
        write (nout,*)
        write (nout,*) ' Effective diffusivities from eigenvectors'      &  
     &    ,' normalized by omega_{De} / k_y^2'
!
        write (nout,130)
        write (nout,132) (chieff(j1),j1=1,imatrx)
!
      endif
!
!
      return
!
 130    format (/t3,'chi_i',t16,'dif_h',t29,'chi_e'                      &  
     &    ,t42,'dif_Z',t55,'chi_Z')
 132    format (1p11e12.4)
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
