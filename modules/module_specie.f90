      MODULE SPECIE
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER :: imppel,impbnd
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 :: begcsp
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  nq,nqo
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  dperi, dpari, ainz, rec
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  rad
      REAL*8 :: radpel, vxpel, vzpel, xpel, zpel, dndtpel ,              &  
     &          tepel , totradp , mwght
      REAL*8 :: fincsp
!     ------------------------------------------------------------------
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SPECIE
