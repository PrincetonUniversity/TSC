      MODULE power_mod
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
!     power.inc --------------------------------------------------------
!     -                                                                |
      REAL*8    PRay(NVELDIM, NPSIDIM), Pql(NVELDIM, NPSIDIM),           &  
     &        PRaytot(NPSIDIM),       Pqltot(NPSIDIM),                   &  
     &        PrIntgrl(NPSIDIM),      PqIntgrl(NPSIDIM),                 &  
     &        PraySum,                PqlSum,             PPwrSum
!     COMMON / PwrQLCom /
!    ^        PRay,                   Pql,
!    ^        PRaytot, Pqltot, PrIntgrl, PqIntgrl,
!    ^        PraySum,                PqlSum,             PPwrSum
 
!     -                                                                |
!     power.inc ends ---------------------------------------------------
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE power_mod
