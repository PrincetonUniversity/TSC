      MODULE WALLP
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 :: begcwp
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  cpw, tc
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  tcond, csubp
      REAL*8 :: tcmlt,rhol,cpl,rhos,heatf,tmelt,tsps,dcps1,dcps2,fcps1,  &  
     &          fcps2,tshfg,achfg1,achfg2,bchfg1,bchfg2,cchfg1,cchfg2,   &  
     &          tevap
      REAL*8 :: fincwp
!     ------------------------------------------------------------------
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE WALLP
