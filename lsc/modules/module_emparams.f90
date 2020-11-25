      MODULE emparams
      IMPLICIT NONE
!
!     emparams.inc  -----------------------------------------------------------
!
      INTEGER MAXPOINTS, SPACEDIM, NPIXDIM, NRDIM, NZDIM, NRZMXDIM,      &  
     &        NCHORDIM
      REAL*8    ACCEPTABLE_ERROR
      PARAMETER(MAXPOINTS = 1000, ACCEPTABLE_ERROR = 1.D-3,              &  
     &        SPACEDIM = 3, NPIXDIM = 50, NCHORDIM = 6,                  &  
     &        NRDIM = 90, NZDIM = 90, NRZMXDIM=90)
!
!     emparams.inc  ----------------------------------------------------------
!
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE emparams
