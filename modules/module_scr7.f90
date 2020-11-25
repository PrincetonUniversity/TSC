      MODULE SCR7
      USE PARAM
      USE SCRATCH
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  aa1
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  dert, derb
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  derl, derr
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR7
