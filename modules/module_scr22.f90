      MODULE SCR22
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER :: mplmax,kgrplmax
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  kgrouplw, jplexv
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 :: plsegmax
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xplw1,zplw1,xplw2,zplw2,     &  
     &                               polsegpl
!============
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR22
