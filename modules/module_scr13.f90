      MODULE SCR13
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  a
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xnf
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  ynf
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR13
