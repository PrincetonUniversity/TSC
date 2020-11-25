      MODULE SCR21
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER :: mmax
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  ipath 
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xpc1,zpc1,xpc2,zpc2,orient,  &  
     &                             polsegc,polwir
!============
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR21
