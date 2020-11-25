      MODULE SCR10
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8,ALLOCATABLE, DIMENSION(:) ::  prpest2,pppest2,ajpest2
      REAL*8,ALLOCATABLE, DIMENSION(:) ::  ppest,gpest
      REAL*8,ALLOCATABLE, DIMENSION(:,:) ::  opest
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR10
