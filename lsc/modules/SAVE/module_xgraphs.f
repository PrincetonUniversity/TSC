      MODULE xgraphs
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
c
c     xgraphs.inc --------------------------------------------------------
c     INTEGER*4 NWKDIx
      INTEGER   NWKDIx
      PARAMETER(NWKDIx = 1000)
c     INTEGER*4 iwkarx(NWKDIx)
      INTEGER   iwkarx(NWKDIx)
c     REAL*4 wkarx(NWKDIx, 4)
      REAL*8   wkarx(NWKDIx, 4)
!     COMMON / grecom1 / iwkarx
!     COMMON / grecom2 / wkarx
c     xgraphs.inc --------------------------------------------------------
c
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE xgraphs
