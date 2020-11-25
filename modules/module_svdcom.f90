      MODULE SVDCOM
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
!     INTEGER nptsp,nordp
!============
! idecl:  explicitize implicit REAL declarations:
!============
!
!     commons and parameters for svd routine
!
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  am
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  uu
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  vv
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  datav
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  sigma, rv1, coef
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  rcocom, zcocom, rnormv,      &  
     &                                 znormv
      REAL*8 :: rs, rsi, zs
!
      REAL*8 :: eps1a,eps1c,eps10,eps2a,eps2c,eps20,eps3a,eps3c,eps30,   &  
     &          eps4a,eps4c,eps40
!
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SVDCOM
