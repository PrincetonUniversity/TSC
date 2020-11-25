 
 
!
!     block data programs                   ---------------------------|
!                                                                      |
!                                                                      |
      BLOCK DATA AcDcBk
!
!     USE CGSetc
!     USE Doflags
!     USE DqlBins
!     USE Escan
!     USE FeBins
!     USE Jrf
!     USE MKSetc
!     USE params
!     USE PIetc
!     USE PlPr
!     USE Ramppwr
!     USE RayBins
!     USE RayWrk
!     USE xparams
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     DATA    PI    , TWOPI , RTPI  , PIO4  , R3O2  , R1O2    /
!    ^        3.1416_R8, 6.2832_R8, 1.77_R8,0.7854_R8,1.22475_R8,
!    &        0.707107_R8/
 
!     DATA    CLIGHT, ECOULB, ELECMS, PROTMS /
!    ^        2.9979_R8, 1.6022_R8, 9.1095_R8, 1.6726_R8/
!            10^8 m/s 10^-19c 10^-28g 10^-24g
 
!     .                                 iSMO = 2 --> Do smoothing as
!     .                                 EValeo did in original code
!     .                                 using smoothed Dql=Dql(iv,ip,2)
!     .                                 iSMO = 1 -->
!     .                                        Use unsmoothed Dql(.,.,1)
!     .                                 and smoothed dfdv(iv,ip) for
!                                              Pql & Jrf
!                smooth Dql is ALWAYS used in calculating fe(iv,ip,iITR)
!
!     DATA    CLIcgs, ECHARG, EMASS,   TeUnits /
!    ^        3.0E10_R8,4.8E-10_R8, 9.1E-28_R8, 1000._R8/
!     DATA    iSMO,      nv /
!    ^           1, NVELDIM /
 
!     DATA    iLastRy / -1 /
!     DATA    lfast, nfreq, npsi, nstep /
!    ^          0,   100, NPSIDIM,20000 /
!     DATA      nrays,   ntors,   npols /
!    ^        NRAYDIM, NTORDIM,       1 /
!     DATA    nzones  /
!    ^        NZONDIM /
!     DATA    PhaseDeg(1), PhaseDeg(2), PhaseDeg(3), nslices /
!    ^                90._R8,        180._R8,        135._R8,     301 /
!     DATA    couplers(1), couplers(2), couplers(3)          /
!    ^        'PBXMFAST' , 'PBXMSLOW' , 'TOKDEVAR'           /
!
!     DATA     nrampup,    nflat /
!    ^        NRAMPDIM, NFLATDEF /
!     DATA    vnmax, dEdcAmnt / 0.99_R8, 0.0001_R8/
!     DATA    DiffuJrf, PrfSpred / 0.000_R8, 0.000_R8/
!
!     DATA    enpol,  fghz, HstpLH /
!    ^           0._R8,   4.6_R8, .005_R8/
!     DATA    nparmax, nparmin, thet0, ScatKdeg, TurnNegs /
!    ^           5.5_R8,     2.5_R8,   0.0_R8,     0.00_R8,       0  /
!     DATA    npolmax, npolmin /
!    ^           +1.0_R8,    -1.0_R8/
!     DATA    nGrps  , centers(1), centers(2), centers(3) /
!    ^           3   ,      4.0_R8,      4.0_R8,      4.0_R8/
!     DATA             widths(1) , widths(2) , widths(3)  /
!    ^                      1.0_R8,      1.0_R8,      1.0_R8/
!     DATA             powers(1) , powers(2) , powers(3)  /
!    ^                      1.0_R8,      0.1_R8,      0.1_R8/
!     DATA
!    ^     RAYPL,    SPECPL,   RFDPL,    RFDPSPL,  DAMPL,    JRFPL,
!    ^     PITPRFPL, DQLFEPL
!    ^     /
!    ^     1,        2,        3,        4,        5,        6,
!    ^     7,        8
!    ^     /
!
!     DATA PlFlg      /
!    ^     nPlFlg * 0 /
!     DATA
!    ^     RAYWR,    FSTFRCWR, NPAPWRWR
!    ^     /
!    ^     1,        2,        3
!    ^     /
!     DATA PrFlg      /
!    ^     nPrFlg * 0 /
!
!     DATA idiag(1),idiag(2),idiag(3),idiag(4),idiag(5) /
!    ^     NRAMPDIM,      0 ,      0 ,      0 ,      0  /
!
!     DATA idiag(6),idiag(7),idiag(8),idiag(9),idiag(10)/
!    ^           0 ,      0 ,      0 ,      0 ,      0  /
!     DATA DoBram, DoTRAN /
!    ^          1,      0 /
!     DATA DoXcam, Do1Rpr, Do0Edc /
!    ^          0,      0,      0 /
!     DATA FeWind               /
!    ^     1, NVELDIM, 1, NPSIDIM /
!     DATA DqlWind              /
!    ^     1, NVELDIM, 1, NPSIDIM /
!     DATA   nsmoo,    nsmw /
!    ^     NSMODEF, NSMWDEF /
!     DATA Vmin, Vmax, WeghtItr /
!    ^      -1._R8,  1._R8,     0.20_R8/
!     DATA iEdc, EdcInp / 0, 0.001_R8/
 
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
