!     -----------------------------------------------------------------
      SUBROUTINE DqlInit
!     initialize d quasilinear constants
      USE CGSetc
      USE DqlBins
      USE FeBins
      USE MKSetc
      USE params
      USE PIetc
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8                                                             &  
     &     EcgsbyEmks
      DATA EcgsbyEmks / 3.333333E-5_R8/
!     Valeo original:
!     DqlNorm = (PI / 2.) * ((ECHARG / EMASS) ** 2) * EcgsbyEmks *
!    ^     EcgsbyEmks / (omega * omega) * (3.e10)
!     see Valeo-Eder, Eq. (22)
!     but why did he have omega twice
!     Here is an MKS version based on
!     Dql = \Sigma (\pi/2) q^2/m^2 E_z^2 /({\Delta k_\parallel} v_\parallel)
!         = \Sigma (\pi/2) q^2/m^2 E_z^2 / omega * v/delta-v
!         = \Sigma DqlNorm E_z^2 /(v/c *delta-n)
      DqlNorm = (PI/2._R8) * (ECOULB/ELECMS)**2 * 1.0E24_R8/ omega
!     Now divide out c^2
      DqlNorm = DqlNorm / (CLIGHT * 1.0E08_R8)**2
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
