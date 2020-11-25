      MODULE PIetc
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
      INTEGER, PRIVATE, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

!     File: PIetc.inc                           starts
      REAL*8                                                             &  
     &        PI    , TWOPI , RTPI  , PIO4  , R3O2  , R1O2
      DATA    PI    , TWOPI , RTPI  , PIO4  , R3O2  , R1O2    /          &  
     &        3.1416_R8, 6.2832_R8, 1.77_R8,0.7854_R8,1.22475_R8,        &  
     &        0.707107_R8/

!     COMMON /PIetc/
!    ^        PI    , TWOPI , RTPI  , PIO4  , R3O2  , R1O2
!     File: PIetc.inc                           ends
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE PIetc
