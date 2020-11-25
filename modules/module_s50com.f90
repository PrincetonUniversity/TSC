      MODULE S50COM
      USE PARAM
      USE SCRATCH
      IMPLICIT NONE
!============
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  iv,itera
      REAL*8,  ALLOCATABLE, DIMENSION(:) ::  r1,z1,r2,z2,gr,gz,gg,grz,   &  
     &                                 gzz,grr,betaa,                    &  
     &                           wk,alfr,alfi,v
!============
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE S50COM
