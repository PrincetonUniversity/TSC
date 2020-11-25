      MODULE MKSetc
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
      INTEGER, PRIVATE, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     File: MKSetc.inc                          starts
      REAL*8                                                             &  
     &        CLIGHT, ECOULB, ELECMS, PROTMS

      DATA    CLIGHT, ECOULB, ELECMS, PROTMS /                           &  
     &        2.9979_R8, 1.6022_R8, 9.1095_R8, 1.6726_R8/

!     COMMON /MKSetc/
!    ^        CLIGHT, ECOULB, ELECMS, PROTMS
!     DATA    CLIGHT, ECOULB, ELECMS, PROTMS /
!    ^        2.9979, 1.6022, 9.1095, 1.6726 /
!            10^8 m/s 10^-19c 10^-28g 10^-24g
!     File: MKSetc.inc                          ends
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE MKSetc
