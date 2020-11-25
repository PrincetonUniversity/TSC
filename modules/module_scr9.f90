      MODULE SCR9
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  izer,jzer
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) ::  ans1v,sns1v
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) ::  ans2v,sns2v
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  r1,z1,r2,z2,g1,g2,g3,g4,g5,  &  
     &                                 g6,                               &  
     &                            plgf1,plgf2,signzer
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  sumi
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  sumj
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  gfun4

!============
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR9
