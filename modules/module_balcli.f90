      MODULE BALCLI
      USE PARAM
      IMPLICIT NONE
!============
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alfa1,alfa2,alfav,beta1,   &  
     &                                 betav,                            &  
     &                                 gamav,gama1
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  di
      REAL*8                        :: baleh
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  idi, idr, idn, node1, idf
      INTEGER                       :: ibal1, ibal2, ifcount
!============
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE BALCLI
