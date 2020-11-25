      MODULE RUNAWAY
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
!--------
!       ajpre   - cell vertice runaway electron current density (advan23)
!       ajprecc - cell center runaway electron current density
!       ajphisf - tsc surface average current density
!       ajpresf - surface average runaway current density
!       anre    - surface average runaway density
!       adnre   - anre*vp(j)
!       recur   - runaway current
!       sresf   - surface average source term
!       sreav   - avalanche production
!       sumre   - total runaway production
!--------
!
      REAL*8 :: begrun
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  ajpre, ajprecc
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  ajpresf, anre, sresf,        &  
     &                                 etafac, adnre
      REAL*8 :: recur, sumre, sreav
      REAL*8 :: finrun 
 
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE RUNAWAY
