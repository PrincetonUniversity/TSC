      MODULE FIXSHP
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 ::  begfix
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  amatsc, uuscfb, vvscfb
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  bvecsc,gsumg,grsumg,gzsumg,  &  
     &                              grrsumg,grho,grho2,sigma,xvecsc,     &  
     &                              xvecsco,bvecsco,rv1scfb
      REAL*8 :: finfix
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE FIXSHP
