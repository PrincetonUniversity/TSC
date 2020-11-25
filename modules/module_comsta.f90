      MODULE COMSTA
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER :: ipl
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 :: fstab, fdestab
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  tpl,fstabpl,fdestpl,fsumpl,  &  
     &                                 taupl,                            &  
     &                             cpasopl,cpasupl
!============
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE COMSTA
