      MODULE numerics
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
!     numerics.inc
      INTEGER npoints_max
      REAL*8    eps, PI, dlxray

!cj   DATA npoints_max / MAXPOINTS /
!cj   DATA eps / ACCEPTABLE_ERROR /
!cj   DATA dlxray / 0.005_R8/


!     COMMON / numcom0 / npoints_max
!     COMMON / numcom1 / eps, PI, dlxray
!     numerics.inc
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE numerics
