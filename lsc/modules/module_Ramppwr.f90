      MODULE Ramppwr
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
!     Ramppwr.inc begins
!     NrampDIM  Max number of different power levels in the
!               ramp up process aimed at iterating to f_e(v_\parallel)
!
!     NrampUp   Num of ramp ups including the flat portion starting from
!                    a Maxwellian
!     Nflat     Num of iterations during the ramp which are flat
!     PwrLevel  Array of power levels. First increasing; then flat.
 
      INTEGER NRAMPDIM, NFLATDEF
      PARAMETER(NRAMPDIM = 200)
      PARAMETER(NFLATDEF =  10)
      INTEGER nRampUp, nFlat
      DATA     nrampup,    nflat /                                       &  
     &        NRAMPDIM, NFLATDEF /

      REAL*8 pwrlevel(NRAMPDIM), FeCvgAry(NRAMPDIM)
!     COMMON / rampc0 / nRampUp, nFlat
!     COMMON / rampc1 / pwrlevel, FeCvgAry
!     Ramppwr.inc ends
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE Ramppwr
