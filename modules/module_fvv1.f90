      MODULE FVV1
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 :: plascur, dipdt, vvl, vvr, cvvpol, vvfin
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  bpolr, bpolz, fwirer, fwirez
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  vvdum
!     ------------------------------------------------------------------
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE FVV1
