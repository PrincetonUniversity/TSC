      subroutine icalc
!......3.30 icalc
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ifcalc,ii,i,jmax,j,imin,imax
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 gisa,fac,psii,xc,gisc,xb,ajb,gjsb,psij
      REAL*8 xd,ajd,gjsd,sumx,sumz,sumxx,sumzz,sumxz,curf
      REAL*8 sum, sumcl
!============
!     dimension xl(5),pl(5),faccoil(pncoil)
      data ifcalc/1/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xl
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pl
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: faccoil
      REAL*8, DIMENSION(pncoil) :: faccoil
!============      
      IF(.not.ALLOCATED(xl)) ALLOCATE( xl(5), STAT=istat)
      IF(.not.ALLOCATED(pl)) ALLOCATE( pl(5), STAT=istat)
!     IF(.not.ALLOCATED(faccoil)) ALLOCATE( faccoil(pncoil), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : icalc  ' 
!============      
!
!     if(ifcalc .eq. 0) go to 801
!     ifcalc = 0
! 801 continue
      do 802 ii=1,nwire
      faccoil(ii) = 1._R8
      if(isym.eq.1 .and. zcoil(ncoil-nwire+ii) .ne. 0)                   &  
     &   faccoil(ii) = 2._R8
  802 continue
!cj oct-28-2009 moved to fix itrmod=26 bug:  801 continue
!...........................................................
!
!.....calculate total toroidal current inside grid using
!.....integral form of ampere's law.
!
!.....let ai4 = total toroidal grid current
!
!...........................................................
!
      ai3o = ai3
      aplo = apl
      ai3  = 0._R8
      ai4  = 0._R8
      apl  = 0._R8
!
!...right boundary :
!
      i      = nx
      gisa = deez/(deex*xarh(i+1))
      jmax   = nz
      fac = 1._R8
      if(isym.eq.0) go to 809
      fac = 2._R8
      ai4 = ai4 + (psi(i+1,2)-psi(i,2))*gisa
  809 continue
      do 810 j=3,jmax
      psii   = psi(i+1,j)-psi(i,j)
      ai4    = ai4 + fac*(psii*gisa )
  810 continue
!
!...left boundary :
!
      i      = 3
      xc = xarh(i)
      gisc = deez/(deex*xarh(i))
      jmax   = nz
      fac = 1
      if(isym.eq.0) go to 819
      fac = 2._R8
      ai4 = ai4 - (psi(i,2)-psi(i-1,2))*gisc
  819 continue
      do 820 j = 3,jmax
      psii   = psi(i,j)-psi(i-1,j)
      ai4    = ai4 - fac*(psii*gisc )
  820 continue
!
!...top boundary :
!
      j      = nz
      imin   = 3
      imax   = nx
      fac = 1._R8
      if(isym.eq.1) fac = 2._R8
      do 830 i = imin,imax
      xb     = .5_R8*(xary(i)+xary(i))
      ajb    = dxdz*xb
      gjsb   = (deex**2)/ajb
      psij   = psi(i,j+1) - psi(i,j)
      ai4    = ai4 + fac*( psij*gjsb)
  830 continue
      if(isym.eq.1) go to 841
!
!...bottom boundary :
!
      j      = 3
      imin   = 3
      imax   = nx
      do 840 i = imin,imax
      xd     = 0.5_R8*(xary(i)+xary(i))
      ajd    = dxdz*xd
      gjsd   = (deex**2)/ajd
      psij   = psi(i,j)-psi(i,j-1)
      ai4    = ai4-(psij*gjsd)
  840 continue
  841 continue
!...........................................................
!
!.....ai4 is total current in domain
!.....apl is plasma current in standard units
!.....ai6 is total current in coils
!
!...........................................................
      ai4 = ai4*udsi
      ai6 = 0._R8
      do 842 ii=1,nwire
  842 ai6 = ai6 + ccoil(ncoil-nwire+ii)*faccoil(ii)
      ai6 = ai6*udsi
      apl = ai4-ai6
!...........................................................
!
!......calculate moments of current to use in bounda routine
!
!...........................................................
      tcuodtp = tcurdtp
      do 97 i=1,2
      do 97 j=1,2
   97 cmomo(i,j) = cmom(i,j)
!
!.....relax current to zero for lrswtch.gt.0
!->   tcurdtp = (1.-acoef(70))*tcurdtp
!->   do 96 i=1,2
!->   do 96 j=1,2
!->96 cmom(i,j) = (1.-acoef(70))*cmom(i,j)
!
!->   if(lrswtch.gt.0) go to 98
      sum    = 0._R8
      sumx   = 0._R8
      sumz   = 0._R8
      sumxx  = 0._R8
      sumzz  = 0._R8
      sumxz  = 0._R8
      sumcl  = 0._R8
      do 99 i=iminn,imaxx
      do 99 j=jminn,jmaxx
      if(iexv(i,j).eq.1 ) go to 99
      if(whalos.eq.0._R8 .and. iexs(i,j).eq.1) go to 99
      fac = 1.0_R8
      if(isym.eq.1 .and. j.ne.2) fac = 2._R8
      if(psi(i,j).ge.phalo) go to 99
      curf   = fac*ajphi(i,j)*dxdz/(2.0_R8*pi)
      sum    = sum   + curf
      sumx   = sumx  + curf*xary(i)
      sumz   = sumz  + curf*zary(j)*(1-isym)
      sumxx  = sumxx + curf*xary(i)*xary(i)
      sumzz  = sumzz + curf*zary(j)*zary(j)
      sumxz  = sumxz + curf*xary(i)*zary(j)*(1-isym)
      if(psi(i,j).gt.psisep) go to 99
      sumcl = sumcl + curf
   99 continue
      tcurdtp = sum
      tcurdtpcl = sumcl
      if(sum.le.0._R8 .and. (lrswtch.gt.0 .or. igone.gt.0)) go to 107
      if(sum.eq.0._R8) go to 105
      xcurf  = sumx/sum
      zcurf  = sumz/sum
      if(lrswtch.eq.0 .and. igone.eq.0)    then
      if(xcurf.ge. alx) ineg=11
      if(xcurf.le.ccon) ineg=11
      if(zcurf.ge. alz) ineg=11
      if(zcurf.le.-alz) ineg=11
      else
      if(xcurf.ge. alx-2._R8*deex ) xcurf = alx - 2._R8*deex
      if(xcurf.le.xlim) xcurf = xlim
      if(zcurf.ge. (alz-2._R8*deez)) zcurf =  (alz-2._R8*deez)
      if(zcurf.le.-(alz-2._R8*deez)) zcurf = -(alz-2._R8*deez)
                          endif
      cmom(1,1) = sumxx - sum*xcurf*xcurf
      cmom(1,2) = sumxz - sum*xcurf*zcurf
      cmom(2,1) = sumxz - sum*xcurf*zcurf
      cmom(2,2) = sumzz - sum*zcurf*zcurf
   98 continue
      if(iflux.ge.2) go to 102
      cmom(1,1) = 0._R8
      cmom(1,2) = 0._R8
      cmom(2,1) = 0._R8
      cmom(2,2) = 0._R8
  102 continue
      if(iflux.ge.1) go to 103
      tcurdtp = 0.0_R8
      xcurf = xplas
      zcurf = 0._R8
  103 continue
      if(irst1.ne.2 .or. kcycle.gt.0) go to 107
      tcuodtp = tcurdtp
      do 106 i=1,2
      do 106 j=1,2
  106 cmomo(i,j) = cmom(i,j)
  107 continue
!
!.....average current moments over two time steps
      tcurdtp = .5_R8*(tcurdtp+tcuodtp)
      do 104 i=1,2
      do 104 j=1,2
  104 cmom(i,j) = .5_R8*(cmom(i,j)+cmomo(i,j))
!
!
      return
  105 continue
      write(*,*) "from icalc"
      ineg=10
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
