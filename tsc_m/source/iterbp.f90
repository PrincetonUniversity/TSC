      subroutine iterbp(iflagbp)
!
!
!           6.94 surfplot
 
!           6.94 surfplot
!.....6.906 coildr
!......6.70 outpts
!
!......6.80 outpl
!
!
!
 
!
!...routine for calculating the poloidal magnetic field at the PF
!...coils in ITER
!
!  iflagbp=1 contrib from ext. coils only
!         =2 plasma + ext coils
!         =3 plasma+ext coils +structure wires
!
      USE CLINAM
      USE ITER1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iflagbp,nrz,nzz,ipfsrt,i,j,kk,k,l,jj
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 bpcalc,brave,bzave
      REAL*8, DIMENSION(pncoil) :: dg
!============
!     common/iter1/ rz(pncoil,20),zz(pncoil,20),bpolav(pncoil),
!    +     bpolr(5,pncoil,pncoil),bpolz(5,pncoil,pncoil),
!    1     sumr(5),sumz(5),ansr(5),ansz(5),
!    1     rfz(pncoil,5), zfz(pncoil,5),bpolmax(pncoil)
!
!...perform these calculations only once
!
      nrz = 20
      nzz = 20
      if(ipfsrt .ge. 1) go to 901
      ipfsrt=1
!      open(83,file='bpout',status='unknown',iostat=ios)
      open(83,file='bpmax',status='unknown',position='append')
      open(84,file='bpave',status='unknown',position='append')




!
!...calculate coordinates of the Rayleigh quadrature points for each coil
!...and the filament locations inside each coil
!
      do 100 i=1,ncoil-nwire
      if(dxcoil(i).le.0 .or. dzcoil(i).le.0) return
      rfz(i,1) = xcoil(i)+dxcoil(i)/2.0_R8
      rfz(i,2) = xcoil(i)-dxcoil(i)/2.0_R8
      rfz(i,3) = xcoil(i)
      rfz(i,4) = xcoil(i)
      rfz(i,5) = xcoil(i)
      zfz(i,1) = zcoil(i)
      zfz(i,2) = zcoil(i)
      zfz(i,3) = zcoil(i)+dzcoil(i)/2.0_R8
      zfz(i,4) = zcoil(i)-dzcoil(i)/2.0_R8
      zfz(i,5) = zcoil(i)
      do 101 j=1,nrz
      rz(i,j)=xcoil(i)-dxcoil(i)/2.0_R8+(j-.5_R8)*dxcoil(i)/nrz
101   continue
      do 102 j=1,nzz
      zz(i,j)=zcoil(i)-dzcoil(i)/2.0_R8+(j-.5_R8)*dzcoil(i)/nzz
102   continue
100   continue
!
!...calculate the poloidal magnetic field per unit current at
!...each of the 5 Rayleigh quadrature points
!...due to all the filaments in every other coil (as well as the
!...coil itself)
!
      do 200 i=1,ncoil-nwire
      do 201 j=1,ncoil-nwire
      do 210 kk=1,5
      sumr(kk)=0.0_R8
      sumz(kk)=0.0_R8
      do 202 k=1,nrz
      do 203 l=1,nzz
      call gradgf(ineg,0,rfz(i,kk),zfz(i,kk),rz(j,k),zz(j,l),            &  
     &            ansr(kk),ansz(kk))
      sumr(kk)=sumr(kk)+ansz(kk)/rfz(i,kk)/ (nrz*nzz*tpi)
      sumz(kk)=sumz(kk)-ansr(kk)/rfz(i,kk)/ (nrz*nzz*tpi)
203   continue
202   continue
      bpolr(kk,i,j)=sumr(kk)
      bpolz(kk,i,j)=sumz(kk)
  210 continue
201   continue
200   continue
!
901   continue
!
!...using the coil currents, plasma currents, and structure currents,
!...calculate the poloidal magnetic field at the Rayleigh quadrature
!...points and combine to give the volume averaged poloidal field
!...at the PF coil
!
!...PF coils
!
      do 300 i=1,ncoil-nwire
      bpolmax(i)=0._R8
      do 310 kk=1,5
      sumr(kk)=0.0_R8
      sumz(kk)=0.0_R8
      do 301 j=1,ncoil-nwire
      sumr(kk)=sumr(kk)+bpolr(kk,i,j)*ccoil(j)
      sumz(kk)=sumz(kk)+bpolz(kk,i,j)*ccoil(j)
301   continue
      bpcalc=sqrt(sumr(kk)**2+sumz(kk)**2)
      bpolmax(i)=max(bpolmax(i),bpcalc)
  310 continue
 
!
!...plasma
!
      if(iflagbp.eq.1) goto 1234
      bpolmax(i)=0._R8
      do 302 j=iminn,imaxx
      do 303 k=jminn,jmaxx
      if(iexv(j,k) .eq. 1) go to 303
      if(iexs(j,k) .eq. 1) go to 303
      if(psi(j,k) .ge. phalo) go to 303
      do 313 kk=1,5
      call gradgf(ineg,0,rfz(i,kk),zfz(i,kk),xary(j),zary(k),            &  
     &            ansr(kk),ansz(kk))
      sumr(kk)=sumr(kk)+ansz(kk)/rfz(i,kk)*ajphi(j,k)*dxdz/tpi
      sumz(kk)=sumz(kk)-ansr(kk)/rfz(i,kk)*ajphi(j,k)*dxdz/tpi
      bpcalc=sqrt(sumr(kk)**2+sumz(kk)**2)
      bpolmax(i)=max(bpolmax(i),bpcalc)
  313 continue
303   continue
302   continue
 1234 continue
!
!...structure wires
!
      if(iflagbp.eq.3) then
      bpolmax(i)=0._R8
      do 304 j=1,nwire
      jj=ncoil-nwire+j
      k=iwire(j)
      l=jwire(j)
      do 314 kk=1,5
      call gradgf(ineg,0,rfz(i,kk),zfz(i,kk),xary(k),zary(l),            &  
     &            ansr(kk),ansz(kk))
      sumr(kk)=sumr(kk)+ansz(kk)/rfz(i,kk)*ccoil(jj)/tpi
      sumz(kk)=sumz(kk)-ansr(kk)/rfz(i,kk)*ccoil(jj)/tpi
      bpcalc=sqrt(sumr(kk)**2+sumz(kk)**2)
      bpolmax(i)=max(bpolmax(i),bpcalc)
  314 continue
304   continue
      endif
!
      brave = sumr(5)
      bzave = sumz(5)
      do 324 kk=1,5
      brave = brave + sumr(kk)
  324 bzave = bzave + sumz(kk)
      bpolavr(i) = brave/6._R8
      bpolavz(i) = bzave/6._R8
      bpolav(i) = sqrt((brave**2 + bzave**2))/6._R8
300   continue
!
      call gaps3(dg(1),dg(2),dg(3),dg(4),dg(5),dg(6))
      write(83,9111) times
      write(83,9112) (ccoil(i)*udsi*1.e-6,i=1,ncoil-nwire)
      write(83,9112) (bpolmax(i),i=1,ncoil-nwire)
      write(83,9112) (bpolav(i),i=1,ncoil-nwire)
      write(83,9112) (bpolavr(i),i=1,ncoil-nwire)
      write(83,9112) (bpolavz(i),i=1,ncoil-nwire)
      write(83,9113) (dg(i),i=1,6)
      write(83,9112) (vegplot(i)*udsv*1.e-3,i=1,ncoil-nwire)
      write(83,9114) apl, 5.*palpha, xcurf, zcurf, xmag, zmag
      write(83,9111) powerr+poweri
 9111 format(e13.5)
 9112 format(18e13.5)
 9113 format(6e13.5)
 9114 format(6e13.5)
!      write(83,904) name(1),kcycle,times
!  904 format(2x,a8,i10,1pe20.12)
!      write(83,905) (bpolav(i),i=1,ncoil-nwire)
!      write(nterm,905) (bpolav(i),i=1,ncoil-nwire)
!905   format(1p5e16.8)
 
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
