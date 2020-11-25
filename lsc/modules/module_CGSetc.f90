      MODULE CGSetc
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
      INTEGER, PRIVATE, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

!     File: CGSetc.inc                          starts
      REAL*8                                                             &  
     &        CLIcgs, ECHARG, EMASS,                                     &  
     &        dompesqdn,                                                 &  
     &        ERGperEV, keVperERG, RESTenergy

      DATA    CLIcgs, ECHARG, EMASS   /                                  &  
     &        3.0E10_R8,4.8E-10_R8, 9.1E-28_R8 /

!     COMMON /CGSetc/
!    ^        CLIcgs, ECHARG, EMASS,
!    ^        dompesqdn,
!    ^        ERGperEV, keVperERG, RESTenergy
 
!     File: CGSetc.inc                          ends
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE CGSetc
