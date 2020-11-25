!#include "f77_dcomplx.h"
      subroutine sprop
!====================================================================
!.....2.35 advan23
!
!
      USE CLINAM
      USE RADTAB
      USE SAPROP
      USE SCR3
      USE SPECIE
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!  * * this subroutine defines the surfaced averaged proerties of
!  * * ion and electron density and temperature in n/m**3 and ev
!  * * specie one is hydrogen
!
!
!      ----- note -----
!
!      Hydrogen here is treated as having
!      charge zgas and mass amgas(amu)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER jstrt,j,nimp,nchrgs,n, ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sumz,sumzb,sumz2,sumvol,sumq,sumq2,sumq2b
      REAL*8 sumq2s,sumni,sumni0,amuimp,argz,animpj,antot,tedge
      REAL*8 AREAL
!============
!     dimension zeffbs(ppsi)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zeffbs
!============      
      IF(.not.ALLOCATED(zeffbs)) ALLOCATE( zeffbs(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : sprop  ' 
!============      
      jstrt = 1
      if(iimp.eq.0) go to 41
      sumz = 0._R8
      sumzb = 0._R8
      sumz2 = 0._R8
      sumvol = 0._R8
      do 10 j = 2,npsi
      ane(j) = (adn(j)/vp(j))*udsd
         if(acoef(4994) .le. 0 .and. allocated(nbeami)) then
            do ii = 1, nspec_beam
               ane(j) = ane(j) + nbeami(j,ii) * q_snbi(ii)
            enddo
         endif
!
      sumq = 2._R8*anhe(j)
      sumq2 = 4._R8*anhe(j)
      sumq2b = 4._R8*anhe(j)
      sumq2s = 4._R8*anhe(j)/2._R8
      sumni = anhe(j)
!
      do 35 nimp = 1,pimp
      nchrgs = nchrgsr(nimp)
      sumni0 = sumni
      amuimp = (nchrgs-1)*2
      do 30 n = 2,nchrgs
      sumq = sumq+AREAL(n-1)*nq(n,nimp,j)/vp(j)
      sumq2 = sumq2+AREAL(n-1)**2*nq(n,nimp,j)/vp(j)
      sumq2s = sumq2s + AREAL(n-1)**2*nq(n,nimp,j)/vp(j)/amuimp
      sumni = sumni+nq(n,nimp,j)/vp(j)
   30 continue
!
      if(nimp.ge.7) go to 115
      go to (100,105,110,111,112,115),nimp
  100 anox(j) = sumni - sumni0
      go to 115
  105 anca(j) = sumni - sumni0
      go to 115
  110 anfe(j) = sumni - sumni0
      go to 115
  111 anbe(j) = sumni - sumni0
      go to 115
  112 anne(j) = sumni - sumni0
  115 continue
!
   35 continue
      aneimp(j) = sumq
      animp(j) = sumni
      if(j.gt.npsit) aneimp(j)=0.0_R8
      anhy(j) = ((adn(j)/vp(j))*udsd - sumq)/zgas
      sumnqa(j) = sumq
      if(anhy(j).lt.0) then
!     write(6,3111) j,anhy(j),ane(j),sumq
!3111 format(" j,anhy(j),ane(j),sumq",i3,1p3e12.4)
!     stop
      ineg=60
      anhy(j) = 1._R8
      endif
      if(sumni .ne. 0) avezimpa(j) = sumq/sumni
      amassimpa(j) = avezimpa(j)*2._R8
      amasshyda(j) = amgas
      sumni = sumni + anhy(j)
      sumnia(j) = sumni
      aimassa(j) = (amgas*anhy(j) + 4._R8*anhe(j) + (sumq*2._R8))        &  
     &             /sumni
      sumq2 = sumq2 + anhy(j)*zgas**2
         if(acoef(4994) .le. 0 .and. allocated(nbeami)) then
            do ii = 1, nspec_beam
               sumq2 = sumq2 + nbeami(j,ii) * q_snbi(ii) * q_snbi(ii)
            enddo
         endif
      zeffa(j) = sumq2/ane(j)
      sumq2b = sumq2b + anhy(j)*zgas**2
      zeffbs(j) = sumq2b/ane(j)
      sumq2s = sumq2s + anhy(j)*zgas**2/amgas
      zeffa2(j) = sumq2s/ane(j)
!
      ti(j) = udsh*(adp(j)-ade(j))/(vpg(j)*sumni)
      te(j) = udsh*ade(j)/(vpg(j)*ane(j))
      if(te(j).lt.tevv) te(j) = tevv
      if(ti(j).lt.tevv) ti(j) = tevv
      if(whalos.gt.0 .and. te(j).lt.thalos) te(j) = thalos
      if(whalos.gt.0 .and. ti(j).lt.thalos) ti(j) = thalos
      avez(j) = ane(j)/sumni
      if(j.gt.npsit) go to 10
      sumz = sumz + zeffa(j)*vp(j)
      sumzb = sumzb + zeffbs(j)*vp(j)
      sumz2 = sumz2 + zeffa2(j)*vp(j)
      sumvol = sumvol + vp(j)
!
   10 continue
      zeff = sumz/sumvol
      zeffb = sumzb/sumvol
!     dzeff2 = sumz2/sumvol
      te(1) = te(2)
      ane(1) = ane(2)
      avez(1) = avez(2)
      anox(1) = anox(2)
      anca(1) = anca(2)
      anfe(1) = anfe(2)
      anbe(1) = anbe(2)
      anne(1) = anne(2)
      aneimp(1) = aneimp(2)
      anhy(1) = anhy(2)
      zeffa(1) = zeffa(2)
      zeffa2(1) = zeffa2(2)
      ti(1) = ti(2)
!
!
      return
!...................................................................
!.....simplified treatment based on z-effective
!...................................................................
   41 continue
      do 51 j=jstrt,npsi
      ane(j) = (adn(j)/vp(j))*udsd
!
!.....density of ions other than impurities
      argz=(xsv(j)-psimin)/(psilim-psimin)
      if(argz .gt. 1.0_R8) argz=1.0_R8
      if(argz .lt. 0.0_R8) argz=0.0_R8
!
       zeffa(j) = zeff
!.....commented out 12/17/10  (scj)
!                    *(1._R8+acoef(3012)*(1._R8-argz)**acoef(3013)       &  
!     & *(argz)**acoef(3014))
!
      anhy(j) = (ane(j)*(zimp-zeffa(j)) + (4._R8-2._R8*zimp)*anhe(j))    &  
     &        / ((zimp-zgas)*zgas)
      animpj = (ane(j) - zgas*anhy(j) - 2._R8*anhe(j))/zimp
      antot = anhy(j) + anhe(j) + animpj
      te(j) = udsh*ade(j)/(vpg(j)*ane(j))
      tedge = tevv
      if(whalos.gt.0 .or. teflat_time.gt.0) tedge = max(tevv,thalos)
      if(acoef(880) .gt. 0) tedge=acoef(880)
!
      if(te(j) .lt. tedge) te(j)=tedge
      ti(j) = udsh*(adp(j)-ade(j))/(vpg(j)*antot)
      if(ti(j) .lt. tedge*(acoef(882)-1._R8)) then
        ti(j) = tedge*(acoef(882)- 1._R8)
      endif
      avez(j) = ane(j)/antot
      zeffa2(j) = (anhy(j)*zgas**2/amgas + anhe(j) + animpj*zimp/2._R8)  &  
     &          / ane(j)
!     zeffa(j) = zeff
      aneimp(j) = (2._R8*anhe(j) + zimp*animpj)*usdd
      sumnqa(j) = (2._R8*anhe(j) + zimp*animpj)
      animp(j) = (anhe(j) + animpj)
      sumnia(j) = animp(j) + anhy(j)
      avezimpa(j) = zimp
      amassimpa(j) = zimp*2._R8
      amasshyda(j) = amgas
      aimassa(j) = (amgas*anhy(j) + 4._R8*anhe(j) + (zimp*2._R8)*animpj)  &  
!    &                                                                   &  
     &           / antot
         if(acoef(4994) .le. 0 .and. allocated(nbeami)) then
            do ii = 1, nspec_beam
               ane(j) = ane(j) + nbeami(j,ii) * q_snbi(ii)
            enddo
         endif
   51 continue
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
