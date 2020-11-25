      subroutine ohmcalc( pohmp,pohmh)
!
!......calculate the ohmic power dissipated in plasma and halo
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pohmh,pohmp,sump,sumh,pohm
!============
      sump = 0._R8
      sumh = 0._R8
      do 551 i=2,nxp
      do 550 j=2,nzp
!
!.....ohmic power dissipated at that point
      pohm = etay(i,j)*ajphi(i,j)**2*udsr*udsi**2*deex*deez*tpi*xary(i)
      if(isym.eq.1 .and.j.ne.2) pohm = pohm*2._R8
!.....check if outside of vacuum vessel
      if(iexv(i,j).eq.1) go to 500
!.....check if on wrong side of separatrix
      if(iexs(i,j).eq.1) go to 532
      go to 498
  532 continue
      if(abs(psi(i,j)-psilim) .lt. abs(phalo-psilim)) go to 499
      go to 500
  498 if(psi(i,j).gt.phalo) go to 500
      if(psi(i,j).gt.psilim) go to 499
      if(igone.eq.1) go to 499
!
!.....point lies in plasma
      sump = sump + pohm
      go to 550
!
!.....point lies in halo
  499 continue
      sumh = sumh + pohm
      go to 550
  500 continue
!
!.....point lies in vacuum
  550 continue
  551 continue
!
      pohmp = sump
      pohmh = sumh
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
