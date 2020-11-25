      MODULE NONCOR
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 :: begcim
      INTEGER :: nby
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  nzspec
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  fracim
      REAL*8 :: fincim
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE NONCOR
