!     -----------------------------------------------------------------
      SUBROUTINE FeConst
!                                       constants needed in solution of fe
!     In the Valeo-Eder paper
!     \nu_0 = \beta_z lnLambda 4 \pi e^4 n /( m^2 v_t^3)
!
!     where \beta_z is of order 1/2 and is about (1+Z)/5
!     where lnLambda is the Coulomb logarithm
!     so
!     \nu_0 = \beta_z lnLambda \cdot
!             TWOPI/5  1.6^4  3.0^1 / 9.11^2   n_{13}/(v/c)^3
      USE CGSetc
      USE DqlBins
      USE FeBins
      USE MKSetc
      USE params
      USE PIetc
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      VthNorm = sqrt( TeUnits * (ERGperEV/1.E-12_R8) /                   &  
     &              (EMASS/1.E-28_R8) )                                  &  
     &             /  (CLIcgs/1.E10_R8)    * 1.E-2_R8
      fe0 = 1._R8/ sqrt(TWOPI)
!     nu0 = betaZ * lnLambda * 4. * PI * ((ECHARG/1.e-10) ** 4) /
!    ^     ((EMASS/1.e-28) * * 2) * 1.e+16 Enright took out this for this:
      nu0 = 4._R8* PI * ((ECHARG/1.E-10_R8) ** 4) /                      &  
     &     ((EMASS/1.E-28_R8) * (EMASS/1.E-28_R8)) * 1.E+16_R8
!cj  &     ((EMASS/1.E-28_R8) * * 2) * 1.E+16_R8
 
      PwrNorm = ELECMS * CLIGHT * CLIGHT * 1.E-31_R8* 1.E+16_R8*         &  
     &          1.E+6_R8
!     normalization for computation of QL power deposition
!     (see files [FW]power.F, subroutine RFheat[Ele])
!     The heating rate is:
!     3/2 n dT/dt = \int S_w \p \eps/\p v dv^3
!     where S_w is the wave induced flux in velocity space, and
!     where \eps is the energy per particle  == 1/2 mv^2
!     S_w = - D_{QL} \cdot \p f/\p v
!     so heating  = \int D_{QL} mv df/dv dv^1 (integrating over v-perps)
!                 = \int mc^2 (D_{QL} v d(cf)/dv dv (normalizing D and v to c)
!     The expression for PwrNorm gives mc^2 and then a conversion to m^-3
!     because f has density in it at cm-3 units.
!
!
      nuNorm =    nu0 / ((CLIGHT * 1.E10_R8) ** 3)
      DcollNorm = nuNorm
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
