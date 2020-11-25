      MODULE SCR8
      USE PARAM
      USE SCRATCH
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  afac, bfac, cfac
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  sn, cn
!     REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  psi 
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  x 
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  z 
!     REAL*8, ALLOCATABLE, DIMENSION(:) ::  dert, derb 
!     REAL*8, ALLOCATABLE, DIMENSION(:) ::  derl, derr, derlt, derrt
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  derlt, derrt
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  abnd, bbnd
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  bv, cv
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR8
