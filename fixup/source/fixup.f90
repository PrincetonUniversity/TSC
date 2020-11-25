      subroutine fixup
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     INTEGER missionc,missionw,lsc
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 reginit,rawmeas,regler,vforce,restcon,vforcepl,fpplot
!     REAL*8 appvolto,spdpbx,spdtftr,spdpbxm,spdpbxn,fedtsc
!     REAL*8 spdasdex,getlscrestart,putlscrestart,growth,ripplot
!     REAL*8 colorc,ripple,d02bae,clock
!============
      entry reginit
!     entry walld
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
      entry getlscrestart
      entry putlscrestart
      entry growth
      entry ripplot
!     entry lingplt
!     entry lingpl2
!     entry lingpl3
      entry colorc
!     entry balloon
      entry ripple
      entry lsc
!     entry tridiag
      entry d02bae
!     for double precision nag use d02baf
      entry clock
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
