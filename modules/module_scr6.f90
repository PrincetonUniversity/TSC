      MODULE SCR6
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER, ALLOCATABLE, DIMENSION(:) ::   jtemp
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  wrk1, wrk3
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  wrk2
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xtmp1,ztmp1,xtmp2,ztmp2,     &  
     &                             atrn1,atrn2,dxtmp1,dztmp1,            &  
     &                             dxtmp2,dztmp2,                        &  
     &                             g4,g5,g6,g1,g2,g3
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  atans
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  a1,a1i,a2,a3,a4
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  rr, vdt,veg0,vdiff 
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR6
