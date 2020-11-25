      MODULE SCRATCH
      USE PARAM 
      IMPLICIT NONE
!============
      INTEGER :: nbig, pnxz1, pnxz2, pnxz3
      INTEGER :: pnxz4,pnxz5,pnxz6,pnxz7,pnxz8,pnxz9,                    &  
     &        pnxz10,pnxz11,pnxz12,pnxz13,pnxz14,pnxz15,pnxz16,          &  
     &        pnxz17,pnxz18,pnxz19,pnxz20,pnxz21
      INTEGER :: pnxx1

      REAL*8, ALLOCATABLE, DIMENSION(:) ::  bigcom
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  vecx, vecz, bo2
!

! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCRATCH

