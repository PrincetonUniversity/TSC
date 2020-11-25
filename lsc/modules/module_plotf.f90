      MODULE plotf
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
!     File: Plotf.inc   begins
      INTEGER NDIM, NLEVELS, LINEARSC, LOGSC, YSHIFT, NTDIM, GEOM,       &  
     &     LEFT, RIGHT, UPPER, LOWER, ULC, URC, LLC, LRC
      PARAMETER(NDIM = 100, NLEVELS = 50,                                &  
     &          LINEARSC = 1, LOGSC = 2, GEOM = 3,                       &  
     &          YSHIFT = 50, NTDIM = 1000,                               &  
     &          LEFT = 1, RIGHT = 2, UPPER = 3, LOWER = 4, ULC = 5,      &  
     &          URC = 6, LLC = 7, LRC = 8)
      INTEGER PlType
!     COMMON / PlotTCom / PlType
!     File: Plotf.inc   ends
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE plotf
