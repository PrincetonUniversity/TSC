      MODULE WALLCL
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 :: begcwa
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  xw,zw,bmagfc,sqg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  g33,g11,g22,g21,g12
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  xcentfc,apzone,aplost,     &  
     &                                 zcentfc,                          &  
     &                                 ripcon,zws,aplostr
      REAL*8 :: fincwa
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE WALLCL
