      subroutine fixup
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     INTEGER missionc,missionw
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 reginit,runaway,rawmeas,regler,vforce,restcon,vforcepl
!     REAL*8 fpplot,appvolto,spdpbx,spdtftr,spdpbxm,spdpbxn,fedtsc
!     REAL*8 spdasdex,growth,ripplot,colorc,ripple,d02bae,clock
!============
      entry reginit
!     entry walld
      entry runaway_sub
      entry rawmeas
      entry regler
      entry vforce
      entry restcon
      entry vforcepl
      entry fpplot
!     entry hyper
!     entry lingin1
!     entry lingin2
!     entry lingin3
      entry appvolto
      entry spdpbx
      entry spdtftr
!     entry spdd3d
      entry spdpbxm
      entry spdpbxn
      entry fedtsc
      entry spdasdex
      entry missionc
      entry missionw
      entry growth
      entry ripplot
!     entry lingplt
!     entry lingpl2
!     entry lingpl3
      entry colorc
!     entry balloon
      entry ripple
!     entry tridiag
      entry d02bae
!     for double precision nag use d02baf
      entry clock
      return
      end
! 25May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
