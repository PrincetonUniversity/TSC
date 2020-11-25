      MODULE SCADVAN
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER :: icalcm
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::   iforce
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::   dxa,dxc,dzb,dzd,face
!============
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCADVAN
