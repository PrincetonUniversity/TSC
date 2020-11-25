      MODULE Doflags
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
 
c     file:  Doflags.inc begins
c
c     Contains compute switches and global io switches
c
c     DoBram = 1 computes spectrum from JEStevens Brambilla code
c              0 makes a spectrum out of arbitrary Gaussians
c     DoTRAN = 1 LSC being called by TRANSP, which passes namelist variables
c                thru commons
c     DoXcam = 1 give pictures and plots like the 2d x ray camera
c     Do1Rpr = 1 retrace just 1 ray per call with iRayTrs = 1, after the
c                first call; Do1RayPerCall
c     Do0Edc = 1 zero out Edc as read from lhcdoua/jardin.d using flag
c                iEdc and value EdcInp in Escan.inc
 
 
      INTEGER        DoBram, DoTRAN,
     ^               DoXcam, Do1Rpr, Do0Edc

      DATA DoBram, DoTRAN /
     ^          1,      0 /
      DATA DoXcam, Do1Rpr, Do0Edc /
     ^          0,      0,      0 /
!     COMMON /DoCom0/DoBram, DoTRAN,
!    ^               DoXcam, Do1Rpr, Do0Edc
c
c     file:  Doflags.inc ends
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE Doflags
