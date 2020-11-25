      MODULE SCR3
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ::  npsit6
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  amat,bmat,cmat,emat
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  fvec, dvec, pvec
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  tmp1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  tmp4, tmp5
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  tmp3
!     ------------------------------------------------------------------
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xs6
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  as6, bs6, cs6, ds6
!============
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR3
