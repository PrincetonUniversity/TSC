      subroutine advan23
!=======================================================================
!.....2.35 advan23
!
!.....advance diffusive part of g and fast wave parts of u and g
!.....by substepping over ndiv intervals
      USE CLINAM
      USE POLCUR
      USE RUNAWAY
      USE SCADVAN
      USE SCR1
      USE SCR11
      USE SCR4
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ii,iw,jw,iig,ig,npol
      INTEGER i,j,n,ios92,i1,i2,j1,j2,jminnt,jmaxxt,mdif,mm,m,iii
      INTEGER iabs
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dtt,etag,delt,f1
      REAL*8 g1,g2,g3,g4,polsub,dth,denom,th2mmax,th2m,th2
      REAL*8 v28t,rdenom,fac4,x5,x6,x4,gzx5sq,cgroup1,etatemp,areaz
      REAL*8 rel
!============
!     dimension fcda(pnx,pnz),fcdb(pnx,pnz)
!     dimension volta(pnx,pnz),voltb(pnx,pnz),utemp(pnx,pnz)
!     dimension i1pol(pnwire),j1pol(pnwire),i2pol(pnwire),j2pol(pnwire)
!     dimension polcur_local(pnwire)
!     data ifrstpol /1/
      INTEGER  :: ifrstpol = 1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: etaa,etab,ptsav2,fac3u,     &  
     &                                       fac4u,fac3p,fac4p,fac41p,   &  
     &                                       fac42p,fac43p,fac44p,       &  
     &                                       fac3g,fac4g,psios,pso2s,    &  
     &                                       psave,gos,go2s,gsave,       &  
     &                                       trm1i,trm2i,trm3i
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: fcda
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: fcdb
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: volta
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: voltb
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: utemp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: i1pol
      INTEGER, ALLOCATABLE, DIMENSION(:) :: j1pol
      INTEGER, ALLOCATABLE, DIMENSION(:) :: i2pol
      INTEGER, ALLOCATABLE, DIMENSION(:) :: j2pol
      REAL*8, ALLOCATABLE, DIMENSION(:) :: polcur_local
!============      
      IF(.not.ALLOCATED(fcda)) ALLOCATE( fcda(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(fcdb)) ALLOCATE( fcdb(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(volta)) ALLOCATE( volta(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(voltb)) ALLOCATE( voltb(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(utemp)) ALLOCATE( utemp(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(i1pol)) ALLOCATE( i1pol(pnwire), STAT=istat)
      IF(.not.ALLOCATED(j1pol)) ALLOCATE( j1pol(pnwire), STAT=istat)
      IF(.not.ALLOCATED(i2pol)) ALLOCATE( i2pol(pnwire), STAT=istat)
      IF(.not.ALLOCATED(j2pol)) ALLOCATE( j2pol(pnwire), STAT=istat)
      IF(.not.ALLOCATED(polcur_local)) ALLOCATE( polcur_local(pnwire),   &  
     &                                 STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : advan23  ' 
!============      
!
      if (.NOT.ALLOCATED(etaa)) ALLOCATE (etaa(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(etab)) ALLOCATE (etab(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(ptsav2)) ALLOCATE (ptsav2(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac3u)) ALLOCATE (fac3u(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac4u)) ALLOCATE (fac4u(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac3p)) ALLOCATE (fac3p(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac4p)) ALLOCATE (fac4p(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac41p)) ALLOCATE (fac41p(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac42p)) ALLOCATE (fac42p(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac43p)) ALLOCATE (fac43p(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac44p)) ALLOCATE (fac44p(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac3g)) ALLOCATE (fac3g(pnx,pnz), STAT=istat)
      if (.NOT.ALLOCATED(fac4g)) ALLOCATE (fac4g(pnx,pnz), STAT=istat)
      if (istat .ne. 0) stop "Allocation Error : advan23"

!
      dtt=(dt+dtold)/ndiv
!
!
!.....define etay array for wire points from wire resistance
      do 51 ii=1,nwire
      iw = iwire(ii)
      jw = jwire(ii)
   51 etay(iw,jw) = rswire(ii)*dxdz/(tpi*xary(iw))
      if(kcycle.ge.npert) go to 56
!
!.....enhance coil resistance so plasma sees perturbation
      do 52 ii=1,nwire
      iw = iwire(ii)
      jw = jwire(ii)
      etay(iw,jw)=etav
   52 continue
   56 continue
!
!
!.... max resistivity for coil gap
      etag=etav
      if(dtt.ne.0._R8) etag=resgap*dxdz*xplas/(2.0_R8*dtt)
      if(ngroups.le.0) go to 611
      do 610 iig=1,ngroups
      ig = nogroups(iig)
      gapr(ig) = resgs(ig)
      if(resgs(ig).le.0.0_R8) gapr(ig)=etag*tpi/                         &  
     &                        (usdr*sgroup(ig)*(isym+1)                  &  
     &                                        *dxdz)
  610 continue
  611 continue
      delt=0.0_R8
      f1 = 1._R8
      if(kcycle.le.1000) f1 = (kcycle/1000._R8)
      if(acoef(79) .gt. 0.0_R8) delt = f1*dzerw*acoef(504)*usdv
      deltmks = delt*udsv
!
!
!     SPECIAL FOR HALO CURRENT FEEDBACK
      if(acoef(501).eq.0.0_R8) go to 214
      delt = acoef(501)*gcurfb(15)*udsi*.001_R8*usdv
      if(delt .gt. deex*acoef(504)*usdv) delt = deex*acoef(504)*usdv
      if(delt .lt.-deex*acoef(504)*usdv) delt =-deex*acoef(504)*usdv
!
      if(ifrstpol .eq. 0) go to 132
!     ifrstpol = 0
      npol = 0
      i = int((acoef(502) - ccon )/deex + 2.5_R8)
      j = int((acoef(503) - zzero)/deez + 2.5_R8)
      do 109 n=1,nwire
      if(icoil(i,j).eq.20 .and. icoil(i,j+1).eq.20) then
      npol = npol + 1
      i1pol(npol) = i+1
      j1pol(npol) = j+1
      i2pol(npol) = i
      j2pol(npol) = j+1
      j = j+1
       else
      go to 111
       endif
  109 continue
  111 continue
      do 120 n=1,nwire
      if(icoil(i,j).eq.20 .and. icoil(i+1,j).eq.20) then
      npol = npol + 1
      i2pol(npol) = i+1
      j2pol(npol) = j+1
      i1pol(npol) = i+1
      j1pol(npol) = j
      i = i+1
       else
      go to 121
       endif
  120 continue
  121 continue
      do 130 n=1,nwire
      if(icoil(i,j).eq.20 .and. icoil(i,j-1).eq.20) then
      npol = npol + 1
      i1pol(npol) = i
      j1pol(npol) = j
      i2pol(npol) = i+1
      j2pol(npol) = j
      j = j-1
       else
      go to 131
       endif
  130 continue
  131 continue
!
!--> Diagnostic Printout
      write(nout,3110)
 3110 format("      n        i1      j1      i2      j2")
      do 142 n=1,npol
      write(nout,3111) n, i1pol(n),j1pol(n),i2pol(n),j2pol(n)
 3111 format(5i8)
  142 continue
      if (numargs .lt. 1) then
         filename = 'halo.out'
      else
         filename = 'halo.out' // '.' // trim(suffix)
      end if
      open(92,file=trim(filename),status='unknown',iostat=ios92)
  132 continue
      if(iplt2.le.nskip2) go to 214
      do 133 n=1,npol
      i1 = i1pol(n)
      i2 = i2pol(n)
      j1 = j1pol(n)
      j2 = j2pol(n)
      polcur_local(n) = tpi*udsi*(xsqoj(i1)*g(i1,j1) -                   &  
     &                  xsqoj(i2)*g(i2,j2))
  133 continue
      deltmks = delt*udsv
      write(92,3112) kcycle,npol, deltmks,times,zmag
      write(92,3113) (polcur_local(n),n=1,npol)
 3112 format(2i5,1p3e15.7)
 3113 format(1p5e15.7)
  214 continue
!.....END OF SPECIAL SECTION FOR HALO FEEDBACK
!
!.... define etay arrays for poloidal current calculation
      ilower = 0
      jlower = 0
      iupper = 0
      jupper = 0
      do 50 j=2,nzp
      etay(1,j) = etay(2,j)
      do 50 i=2,nxp
      etaa(i,j) = .5_R8*(etay(i,j)+etay(i,j-1))*dxisq/xary(i)
      etab(i,j) = .5_R8*(etay(i,j)+etay(i-1,j))*dzisq
      ptsav2(i,j) = 0._R8
      fcda(i,j)=.5_R8*(rjcdg(i,j)+rjcdg(i,j-1))
      fcdb(i,j)=.5_R8*(rjcdg(i,j)+rjcdg(i-1,j))
!...for halo current feedback
      volta(i,j) = 0._R8
      voltb(i,j) = 0._R8
      if (acoef(504).eq.0.0_R8) go to 50
      if(icoil(i,j) .eq. 20 .and. icoil(i,j-1) .eq.20) then
      volta(i,j) = -delt
                     endif
      if(icoil(i,j) .eq. 20 .and. icoil(i-1,j) .eq. 20) then
      voltb(i,j) = delt
      if(zary(j).lt.0.0_R8) then
      ilower = i
      jlower = j+1
      else
      iupper = i
      jupper = j
      endif
                                                        endif
      if(icoil(i,j) .eq. 21 .and. icoil(i-1,j) .eq. 21) then
      voltb(i,j) = -delt
      iupper = i
      jupper = j
      iforce(iupper,jupper) = 1
                                                        endif
! kdm - global diagnostic for watching CHI boundary condition
      chivolt(i,j) = volta(i,j) - volta(i-1,j) - voltb(i,j) + voltb(i,j-1)
   50 continue
      if(ilower*jlower*iupper*jupper .ne. 0 ) then
!....diagnostic printout
!....inserted from version 10.9 on 8/9/2013
!
!.....place a point in the middle of the vessel
      g4 = xsqoj(nx/3)*g(nx/3,2*nz/3)
!
!....modified 11/04/10   (scj)
!.....points near the injector region
      g1 = xsqoj(ilower-1)*g(ilower-1,jlower+1)
      g2 = xsqoj(ilower  )*g(ilower  ,jlower+1)
      g3 = xsqoj(ilower+1)*g(ilower+1,jlower+1)
      polcurchi = tpi*udsi*(.5*(g1+g2)-g4)
      polcurchi2= tpi*udsi*(.5*(g3+g2)-g4)
!     g2 = xsqoj(iupper+1)*g(iupper+1,jupper-1)
!     g4 = xsqoj(ilower+1)*g(ilower+1,jlower+1)
!
!.....modified 9/24/07 (scj)
!     g1 = xsqoj(ilower-1)*g(ilower-1,jlower)
!     g3 = xsqoj(ilower-1)*g(ilower-1,jlower-1)
!     polcurchi = tpi*udsi*(g1-g3)
!     polcurchi2= tpi*udsi*(g4-g2)
      polsub = delt*dt
      if( mod(kcycle,100) .eq.0)                                         &  
     &write(nterm,6666) ilower,jlower,iupper,jupper,kcycle,polcurchi,    &  
     &      g1,g2,g3,g4,gzero,polsub
 6666 format(4i4,i7,1p7e12.4)
                                          endif
!
!
!...> diag output
!     do j=3,nzp
!     write(nterm,6667) j,(iforce(i,j),i=2,nxp)
!     enddo
!6667 format(i3,2x,80i1)
!     write(nterm,6668) phalo,psimin,psilim
!     psidiagmin = 1.e6
!     psidiagmax = -1.e6
!     do i=2,nxp
!     do j=2,nzp
!        if(iexv(i,j) .eq. 0) then
!        psidiagmin = amin1(psidiagmin,psi(i,j))
!        psidiagmax = amax1(psidiagmax,psi(i,j))
!                             endif
!     enddo
!     enddo
!     write(nterm,6669) psidiagmin,psidiagmax
!6668 format(" phalo, psimin, psilim =",1p3e12.4)
!6669 format(" psidiagmin, psidiagmax=",1p2e12.4)
      if(lrswtch.ne.0) go to 66
      if(acoef(40).ne.0.0_R8) go to 66
!
!.....zero out terms outside vacuum vessel
      do 49 i=1,nxp
      do 49 j=1,nzp
      ptsave(i,j) = (1-iexv(i,j))*ptsave(i,j)
   49 continue
   66 continue
!
!     insert two lines to call appvolt to use thyrister model
!
      if (acoef(290).eq.1._R8.and. idata .eq. 3) call appvolt
      if (acoef(290).eq.1._R8.and. idata .ne. 3) call appvolto
!
      if (acoef(290).eq.2._R8.or. acoef(290).eq.4._R8) then
        call appvolt2
        go to 54
                endif
!
!
!.....store voltage needed to shoot for cwire0
      if(acoef(290).ne.0.0_R8) go to 155
      do 153 ii=1,nwire
      iw = iwire(ii)
      jw = jwire(ii)
      ptsave(iw,jw)=-etay(iw,jw)*xary(iw)*cwire0(ii)/dxdz
  153 continue
  155 continue
      do 53 ii=1,nwire
      iw = iwire(ii)
      jw = jwire(ii)
      resave(ii) = -ptsave(iw,jw)
   53 continue
   54 continue
!
!
      do 420 j=2,nzp
      psiobl(j) = psi(2,j)
  420 psiobr(j) = psi(nxp,j)
      do 430 i=2,nxp
  430 psiobt(i) = psi(i,nzp)
      if(isym.eq.1) go to 429
      do 434 i=2,nxp
      psiobb(i) = psi(i,2)
  434 continue
  429 continue
!
 
!......obtain boundary values for poloidal flux psi
!            note: this defines the boundary flux from all coils
!                  except the oh coils : returns psibl,psibr,psibb,psibt
!
!
!.....define toroidal current array
!
      if(iflux.ne.4) go to 405
      do 301 i=3,nx
      do 301 j=3-isym,nz
  301 vd(i,j) = 0._R8
!
!.....symmetrize bounds
      jminnt = jminn
      jmaxxt = jmaxx
      if(isym.eq.1) go to 7020
      mdif = min(nzp-jmaxx,jminn-2)
      jmaxxt = nzp - mdif
      jminnt = 2*nh-jmaxxt
 7020 continue
      do 302 i=iminn,imaxx
      do 302 j=jminnt,jmaxxt
  302 vd(i,j) = ajphi(i,j)*xary(i)*dxdz
      do 303 ii=1,nwire
      i = iwire(ii)
      j = jwire(ii)
  303 vd(i,j) = ccoil(ncoil-nwire+ii)*xary(i)
!     write(6,*) " after 303 ", ineg
!     call flush(6)
      call elliptic(4,psizer)
  405 continue
!     write(6,*) " after 405 ", ineg
!     call flush(6)
      call bounda(1)
!     write(6,*) " after bounda ", ineg
!     call flush(6)
!
      dth = 0.5_R8*(dtold+dt)
      if(dth .le. 0) go to 440
      do 431 j=2,nzp
      psidbl(j) = .5_R8*(psidbl(j)+(psibl(j)+psib-psiobl(j))/dth)
  431 psidbr(j) = .5_R8*(psidbr(j)+(psibr(j)+psib-psiobr(j))/dth)
      do 432 i=2,nxp
  432 psidbt(i) = .5_R8*(psidbt(i)+(psibt(i)+psib-psiobt(i))/dth)
      if(isym.eq.1) go to 433
      do 435 i=2,nxp
      psidbb(i) = .5_R8*(psidbb(i)+(psibb(i)+psib-psiobb(i))/dth)
  435 continue
  433 continue
  440 continue
!
!
      dtt = (dt+dtold)/ndiv
      denom = (8._R8*dtt*amux)
      th2mmax = 0.5_R8
      if(denom.gt.0)                                                     &  
     &th2mmax = (deex**2+deez**2)/denom
      th2m = min(0.5_R8,th2mmax)
!
!....test 6.28.98
      th2m= 0._R8
      th2 = 1._R8- th2m
!
      do 185 i=3,nxp
!
      ama6(i) = amux*xary(i)/xarh(i+1)
      amc4(i) = amux*xary(i-1)/xarh(i-1)
      gze6a(i) = gzero*xarh(i+1)/xary(i)
      gze5a(i) = gzero*xarh(i)/xary(i)
      gze4c(i) = gzero*xarh(i-1)/xary(i-1)
      gze5c(i) = gzero*xarh(i)/xary(i-1)
!
      do 186 j=3,nzp
      denom = 1._R8+.5_R8*dtt*face(i,j)*(th2*amux*(                      &  
     & (xary(i)*dxa(i,j)+xary(i-1)*dxc(i,j))/xarh(i)                     &  
     & + dzb(i,j) + dzd(i,j)) + acoef(91) )
!
      fac3u(i,j) = (2._R8-denom)/denom
      fac4u(i,j) = dtt*face(i,j)/denom
!
!     u(i,j) = fac3*uo(i,j) + fac4*(utsave(i,j) +
!    1 dxisq*(ama6*uip(j  )-gze6a*g(i+1,j  )+gze5a*g(i  ,j  ))
!    2+dzisq*(th2*amux*uiz(j+1)-gzero*(g(i  ,j+1)-g(i  ,j  )))
!    3+dxisq*(amc4*uim(j  )+gze5c*g(i  ,j  )-gze4c*g(i-1,j  ))
!    4+dzisq*(th2*amux*uiz(j-1)+gzero*(g(i  ,j  )-g(i  ,j-1))) )
!
  186 continue
  185 continue
!
      do 1400 i=3,nx
      v28t = dzisq
      v4tv(i) = xary(i)*dxisq/xarh(i)
      v6tv(i) = xary(i)*dxisq/xarh(i+1)
      v5tv(i) = -(v4tv(i)+v6tv(i)+2._R8*v28t)
      do 1300 j=3-isym,nz
!
      denom = 1._R8-th*0.5_R8*dtt*v5tv(i)*etay(i,j)
      rdenom = 1._R8/denom
      fac3p(i,j) = (2.0_R8-denom)*rdenom
      fac4 = dtt*rdenom
      fac4p(i,j) = fac4
      fac41p(i,j) = fac4*etay(i,j)*v28t*th
      fac42p(i,j) = fac4*etay(i,j)*v4tv(i)*th
      fac43p(i,j) = fac4*etay(i,j)*v6tv(i)*th
      fac44p(i,j) = fac4*etay(i,j)*thm
!
!     psi(i,j) = (2.0-denom)*rdenom*pso2(i,j)
!    1         + dtt*rdenom*(ptsave(i,j) + etay(i,j)*(
!    2         v28t *(psio(i,j-1)+psio(i,j+1))
!    3         + v4t*psio(i-1,j) + v6t*psio(i+1,j)
!    4         + thm*ajp2(i,j) ) )
!
 1300 continue
 1400 continue
!
      do 1200 i=3,nxp
      x5 = xarh(i)
      do 1100 j=3,nzp
      denom = 1._R8+.5_R8*dtt*th2*((etaa(i,j)+etaa(i-1,j))*x5            &  
     &                  + etab(i,j)+etab(i,j-1)    )
      fac3g(i,j) = (2._R8-denom)/denom
      fac4g(i,j) = dtt/denom
!     g(i,j) = fac3*go2(i,j) + fac4*(gtsave(i,j)
!    1  -gzero*uo(i,j)/x5sq + etaa(i  ,j)*x6*go(i+1,j)
!    2                      + etaa(i-1,j)*x4*go(i-1,j)
!    3                      + etab(i,j  )*go(i,j+1)
!    4                      + etab(i,j-1)*go(i,j-1)  + term )
 1100 continue
 1200 continue
!
!.....start substepping loop on ndiv substeps per time cycle
!
!*******************************************************************
!
      do 500 n=1,ndiv
!
      if(lrswtch.ge.1) go to 201
!
!*******************************************************************
!....................................................................
!
!.....2.10 velocity divergence u
!
!....................................................................
      do 89 j=1,nzp+1
      uiz(j) = u(2,j)
      uip(j) = u(3,j)
   89 continue
      do 85 i=3,nxp
!
      do 96 j=1,nzp+1
      uim(j) = uiz(j)
      uiz(j) = uip(j)
      uip(j) = u(i+1,j)
   96 continue
!
      do 95 j=3,nzp
!
      u(i,j) = fac3u(i,j)*uo(i,j) + fac4u(i,j)*(utsave(i,j) +            &  
     &        dxa(i,j)*(ama6(i)*(th2*uip(j  )+th2m*(uo(i+1,j)-uo(i,j)))  &  
     &                  -gze6a(i)*g(i+1,j)+gze5a(i)*g(i  ,j))            &  
     &       +dzb(i,j)*(amux*(th2*uiz(j+1)+th2m*(uo(i,j+1)-uo(i,j)))     &  
     &                  -gzero*(g(i  ,j+1)-g(i  ,j  )))                  &  
     &       +dxc(i,j)*(amc4(i)*(th2*uim(j  )+th2m*(uo(i-1,j)-uo(i,j)))  &  
     &                  +gze5c(i)*g(i  ,j)-gze4c(i)*g(i-1,j))            &  
     &       +dzd(i,j)*(amux*(th2*uiz(j-1)+th2m*(uo(i,j-1)-uo(i,j)))     &  
     &                  +gzero*(g(i  ,j  )-g(i  ,j-1))))
!
      utemp(i,j) = uiz(j)
   95 continue
      if(isym.eq.1) u(i,2) = u(i,3)
      if(isym.eq.1) utemp(i,2) = utemp(i,3)
!
   85 continue
      do 86 i=3,nxp
      do 86 j=3,nzp
   86 uo(i,j) = utemp(i,j)
!
  201 continue
      do 210 i=3,nxp
      do 110 j=3-isym,nzp
      go2(i,j) = go(i,j)
      go(i,j) = g(i,j)
  110 continue
  210 continue
!....................................................................
!
!.....2.20 toroidal flux g
!
!....................................................................
      do 200 i=3,nxp
      x5 = xarh(i)
      x6 = xarh(i+1)
      x4 = xarh(i-1)
!
      if(acoef(870).ne.0._R8) go to 1101
      do 100 j=3,nzp
      gzx5sq = face(i,j)* gzero/xarh(i)**2
      g(i,j)=fac3g(i,j)*go2(i,j)+fac4g(i,j)*(gtsave(i,j)-gzx5sq *uo(i,j)  &  
!    &                                                                   &  
     & + etaa(i  ,j)*(th2*x6*go(i+1,j)+th2m*(x6*go2(i+1,j)-x5*go2(i,j))  &  
     &     + fcda(i,j)*dxdz                                              &  
     &    *.25_R8*(psi(i+1,j  )+psi(i+1,j-1)-psi(i-1,j  )-psi(i-1,j-1)))  &  
!    &                                                                   &  
     & + etaa(i-1,j)*(th2*x4*go(i-1,j)+th2m*(x4*go2(i-1,j)-x5*go2(i,j))  &  
     &     + fcda(i-1,j)*dxdz                                            &  
     &    *.25_R8*(psi(i-2,j  )+psi(i-2,j-1)-psi(i  ,j  )-psi(i  ,j-1)))  &  
!    &                                                                   &  
     & + etab(i,j  )*(th2*   go(i,j+1)+th2m*(   go2(i,j+1)-   go2(i,j))  &  
     &     + fcdb(i,j)*dxdz/x5                                           &  
     &    *.25_R8*(psi(i-1,j+1)+psi(i  ,j+1)-psi(i-1,j-1)-psi(i  ,j-1)))  &  
!    &                                                                   &  
     & + etab(i,j-1)*(th2*   go(i,j-1)+th2m*(   go2(i,j-1)-   go2(i,j))  &  
     &     + fcdb(i,j-1)*dxdz/x5                                         &  
     &    *.25_R8*(psi(i-1,j-2)+psi(i  ,j-2)-psi(i-1,j  )-psi(i  ,j  )))  &  
!    &                                                                   &  
     & + volta(i,j) - volta(i-1,j) - voltb(i,j) + voltb(i,j-1) )
  100 continue
      go to 1102
 1101 continue
      do 1105 j=3,nzp
      gzx5sq = face(i,j)* gzero/xarh(i)**2
      g(i,j) = fac3g(i,j)*go2(i,j) + fac4g(i,j)*(gtsave(i,j)             &  
     &  -gzx5sq*uo(i,j)                                                  &  
     & + etaa(i  ,j)*(x6*go(i+1,j) + acoef(870)*dxdz                     &  
     &    *.25_R8*(psi(i+1,j  )+psi(i+1,j-1)-psi(i-1,j  )-psi(i-1,j-1)))  &  
!    &                                                                   &  
     & + etaa(i-1,j)*(x4*go(i-1,j) + acoef(870)*dxdz                     &  
     &    *.25_R8*(psi(i-2,j  )+psi(i-2,j-1)-psi(i  ,j  )-psi(i  ,j-1)))  &  
!    &                                                                   &  
     & + etab(i,j  )*(   go(i,j+1) + acoef(870)*dxdz/x5                  &  
     &    *.25_R8*(psi(i-1,j+1)+psi(i  ,j+1)-psi(i-1,j-1)-psi(i  ,j-1)))  &  
!    &                                                                   &  
     & + etab(i,j-1)*(   go(i,j-1) + acoef(870)*dxdz/x5                  &  
     &    *.25_R8*(psi(i-1,j-2)+psi(i  ,j-2)-psi(i-1,j  )-psi(i  ,j  )))  &  
     &  )
 1105 continue
 1102 continue
!
!......impose boundary conditions at this cell if iexvc(i,j).eq.2
      do 101 j=3,nzp
      g(i,j) = (gzero/xsqoj(i))*iforce(i,j) + (1-iforce(i,j))*g(i,j)
  101 continue
      if(isym.eq.1) g(i,2) = g(i,3)
!..boundary condition for antisymmetric toroidal field
      if(jsym.eq.-1) g(i,2) = -g(i,3)
  200 continue
!
      do 410 i=2,nxp
      do 310 j=2-isym,nzp
      pso2(i,j) = psio(i,j)
      psio(i,j) = psi(i,j)
  310 continue
  410 continue
!
      do 316 i=3,nx
      do 315 j=3-isym,nz
      ajp2(i,j) = dzisq*(pso2(i,j-1)+pso2(i,j+1))                        &  
     &          + v4tv(i)*pso2(i-1,j) + v6tv(i)*pso2(i+1,j)              &  
     &          + v5tv(i)*pso2(i,j)
  315 continue
  316 continue
!
!.....accumulate current in wire groups for resistive gap calculation
!..rxw/06/10/87
      if(etag.eq.0._R8) go to 41
!...rxw/end
      do 10 mm=1,ngroup
      m = nogroup(mm)
      cgroup(m) = 0._R8
   10 continue
      if(idata .eq. 3) go to 21
      if(isym.eq.1) go to 25
      do 20 iig=1,ngroups
      ig = nogroups(iig)
      do 20 iii=1,icmaxs(iig)
      ii = icoils(iii,iig)
      i = iwire(ii)
      j = jwire(ii)
      cgroup(ig) = cgroup(ig)+(ajp2(i,j)/xary(i)-cwire0(ii)/dxdz)        &  
     &            *iseries(ig)
   20 continue
      go to 27
   25 continue
      do 26 iig=1,ngroups
      ig = nogroups(iig)
      do 26 iii=1,icmaxs(iig)
      ii = icoils(iii,iig)
      i = iwire(ii)
      j = jwire(ii)
      cgroup(ig) = cgroup(ig)+(ajp2(i,j)/xary(i)-cwire0(ii)/dxdz)        &  
     &            *iseries(ig)
   26 continue
   27 continue
!
!.....exclude single coil groups..
      if(ngroups.le.0) go to 29
      do 30 mm=1,ngroups
      m = nogroups(mm)
   30 cgroup(m) = cgroup(m)/sgroup(m)
      go to 29
   21 continue
!
!.....special for idata=3    (D-III-D comparisons)
      sgroup(1) = 0._R8
      do 22 ii=1,nwire
      i = iwire(ii)
      j = jwire(ii)
      ig = iabs(igroupw(ii))
      if(iseries(ig).eq.0) go to 22
      cgroup(1) = cgroup(1) + ajp2(i,j)/xary(i)
      sgroup(1) = sgroup(1) + 1._R8
   22 continue
      cgroup1 = cgroup(1)
      do 23 iig=1,ngroups
      ig=nogroups(iig)
      cgroup(ig) = cgroup1
      cgroup(ig) = cgroup(ig)/sgroup(1)
   23 continue
   29 continue
!
!
      if(iresgs.eq.1) go to 143
      do 40 ii=1,nwire
      i = iwire(ii)
      j = jwire(ii)
      ig = iabs(igroupw(ii))
!
      etatemp = etag
   40 ptsav2(i,j) = etatemp * cgroup(ig)
      go to 41
  143 continue
      do 140 iig=1,ngroups
      ig = nogroups(iig)
      do 140 iii=1,icmaxs(iig)
      ii = icoils(iii,iig)
      i = iwire(ii)
      j = jwire(ii)
!
      etatemp = gapr(ig)*usdr*sgroup(ig)*(isym+1)*dxdz/tpi
  140 ptsav2(i,j) = etatemp * cgroup(ig)
   41 continue
!
!....................................................................
!
!.....2.30 poloidal flux psi
!
!....................................................................
!
!......advance flux due to oh coils
!
      psi2b = psiob
      psiob = psib
      psib = psi2b + dtt*reboun
!
      if(acoef(870) .ne. 0.0_R8) go to 1401
      do 400 i=3,nx
      do 300 j=3-isym,nz
!
!
!......advance psi (leap frog method)
!
      psi(i,j) = fac3p(i,j)*pso2(i,j)                                    &  
     &       + fac4p(i,j)*(ptsave(i,j)+ptsav2(i,j)-etay(i,j)*rjcd(i,j)   &  
     &                            - etay(i,j)*ajpre(i,j)*xary(i)  )      &  
     &       + fac41p(i,j)*(psio(i,j-1)+psio(i,j+1))                     &  
     &       + fac42p(i,j)*psio(i-1,j) + fac43p(i,j)*psio(i+1,j)         &  
     &       + fac44p(i,j)*ajp2(i,j)
!
  300 continue
      if(isym.eq.1) psi(i,1) = psi(i,3)
  400 continue
      go to 1402
 1401 continue
      do 1405 i=3,nx
      do 1305 j=3-isym,nz
!
!
!......advance psi (leap frog method)
!
      psi(i,j) = fac3p(i,j)*pso2(i,j)                                    &  
     &       + fac4p(i,j)*(ptsave(i,j)+ptsav2(i,j)-etay(i,j)*(rjcd(i,j)-  &  
!    &                                                                   &  
     &         acoef(870)*.25_R8*(xsqoj(i  )*(g(i  ,j)+g(i  ,j+1))       &  
     &                        +xsqoj(i+1)*(g(i+1,j)+g(i+1,j+1))) )       &  
     &                            - etay(i,j)*ajpre(i,j)*xary(i)  )      &  
     &       + fac41p(i,j)*(psio(i,j-1)+psio(i,j+1))                     &  
     &       + fac42p(i,j)*psio(i-1,j) + fac43p(i,j)*psio(i+1,j)         &  
     &       + fac44p(i,j)*ajp2(i,j)
!
 1305 continue
      if(isym.eq.1) psi(i,1) = psi(i,3)
 1405 continue
 1402 continue
!
!
!.....compute relaxation factor for boundary flux
      f1 = acoef(806)
      if(irst1.eq.2 .and. kcycle.le.200) f1 = f1*(kcycle/200._R8)**4
!
!  ---> added 9/13/86
      areaz = max(deex,deez)**2
      rel = min(1._R8,0.1_R8*dtt*etav/areaz)*f1
!
!
!.....let boundary flux relax toward computed value
      do 320 j=2,nzp
      psi(2,j) =(psibl(j) + psib)*rel + psi(2,j)*(1._R8-rel)
  320 psi(nxp,j) =(psibr(j) + psib)*rel + psi(nxp,j)*(1._R8-rel)
      if(isym.eq.1) go to 331
      do 229 i=3,nx
      psi(i,2) =(psibb(i) + psib)*rel + psi(i,2)*(1._R8-rel)
  229 continue
  331 continue
      do 330 i=3,nx
  330 psi(i,nzp) =(psibt(i) + psib)*rel + psi(i,nzp)*(1._R8-rel)
  500 continue
!.....save gap values
      if(ngroups.le.0) go to 511
      do 510 iig=1,ngroups
      ig = nogroups(iig)
      if(sgroup(ig).eq.0) go to 510
      gapi(ig) = cgroup(ig)*dxdz*udsi*sgroup(ig)*(isym+1)*.001_R8
      gapv(ig) = gapr(ig)*gapi(ig)
  510 continue
  511 continue
!
!
!
!     DEALLOCATE(etaa,etab,ptsav2,fac3u,fac4u,fac3p,fac4p,fac41p,fac42p,
!    &           fac43p,fac44p,fac3g,fac4g
!    &           ,psios,pso2s,psave,gos,go2s,gsave,trm1i,trm2i,trm3i
!    &           )
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
