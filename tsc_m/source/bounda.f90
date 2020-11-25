      subroutine bounda(itype)
!***********************************************************************
!                                                                      *
!...supplies boundary conditions                                       *
!                                                                      *
!.....itype = 1 ... poloidal flux "psi"                                *
!           = 2 ... velocity stream function "abig"                    *
!           = 3 ... unused                                             *
!                                                                      *
!***********************************************************************
!
      USE CLINAM
      USE SCR9
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,ncoilte,ider,l,jb,i,n,ib,j,m,mmax,jmax,mm,imin
      INTEGER imax,mmaxmm,iz,jz,indx,ii,jj,jndx
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zaryjb,ans,xaryib,rc,zc,zs,ans2,ans1,xb,zb,gdot,efield
      REAL*8 sum1,sum2,sum3,xsq,gzx,fac,psir,psiz,psimid,psich
      REAL*8 psidot,bsq,vcom,acon
      REAL*8 sum
!============
      go to(1,2,3),itype
!
!
!**********************************************************
!
!....... poloidal flux "psi"
!
!**********************************************************
    1 continue
      ncoilte = ncoil
      if(iflux.eq.4) ncoilte = ncoil-nwire
!...........................................................
!
!...first time thru, calculate and store green's functions
!...needed for discrete coils contribution
!
!...........................................................
      if(ifrst(3).eq.0) go to 10
      ifrst(3) = 0._R8
!
!.....bottom-top boundaries
      do 11 ider=1,2
      do 11 l=1+isym,2
      if(l.eq.1) jb = 2
      if(l.eq.2) jb = nzp
      zaryjb = zary(jb)
      if(ider.eq.2) zaryjb = zary(jb) - (2*l-3)*deez
      do 12 i=2,nxp
      do 13 n=1,ncoilte
      call gf(ineg,nmult,xary(i),zaryjb,xcoil(n),zcoil(n),ans)
      ans1v(i,l,n,ider) = ans
      if(isym.eq.0) go to 13
      call gf(ineg,nmult,xary(i),zaryjb,xcoil(n),-zcoil(n),ans)
      if(zcoil(n).eq.0) ans = 0._R8
      sns1v(i,l,n,ider) = ans
   13 continue
   12 continue
   11 continue
!
!.....left-right boundaries
      do 14 ider=1,2
      do 14 l=1,2
      if(l.eq.1) ib = 2
      if(l.eq.2) ib = nxp
      xaryib = xary(ib)
      if(ider.eq.2) xaryib = xary(ib) - (2*l-3)*deex
      do 15 j=2,nzp
      do 16 n=1,ncoilte
      call gf(ineg,nmult,xaryib,zary(j),xcoil(n),zcoil(n),ans)
      ans2v(j,l,n,ider) = ans
      if(isym.eq.0) go to 16
      call gf(ineg,nmult,xaryib,zary(j),xcoil(n),-zcoil(n),ans)
      if(zcoil(n).eq.0) ans = 0._R8
      sns2v(j,l,n,ider) = ans
   16 continue
   15 continue
   14 continue
!
      if(iflux.ne.4) go to 10
      m = 0
      do 217 ider=1,2
      do 217 l=1+isym,2
      if(l.eq.1) jb = 2
      if(l.eq.2) jb = nzp
      do 218 i=2,nxp
      m = m+1
      r2(m) = xary(i)
      z2(m) = zary(jb)
      if(ider.eq.2) z2(m) = zary(jb) - (2*l-3)*deez
  218 continue
  217 continue
!
      do 219 ider=1,2
      do 219 l=1,2
      if(l.eq.1) ib = 2
      if(l.eq.2) ib = nxp
      do 220 j=2,nzp
      m = m+1
      r2(m) = xary(ib)
      if(ider.eq.2) r2(m) = xary(ib) - (2*l-3)*deex
      z2(m) = zary(j)
  220 continue
  219 continue
      mmax = m
!
!.....right boundary
!
      i = nx
      jmax = nz
      mm = 0
      do 810 j=3-isym,jmax
      mm = mm+1
      izer(mm) = i
      jzer(mm) = j
      signzer(mm) = -deez
      rc = xarh(i+1)
      zc = zary(j)
      zs = -zc
      do 811 m=1,mmax
      ans2 = 0
      call gf(ineg,nmult,r2(m),z2(m),rc,zc,ans1)
      if(isym.eq.0 .or. j.eq.2) go to 812
      call gf(ineg,nmult,r2(m),z2(m),rc,zs,ans2)
  812 gfun4(m,mm) = (ans1+ans2)/(tpi*rc*deex)
  811 continue
  810 continue
!
!.....left boundary
!
      i = 3
      rc = xarh(i)
      jmax = nz
      do 820 j=3-isym,jmax
      mm = mm+1
      izer(mm) = i
      jzer(mm) = j
      signzer(mm) = -deez
      zc = zary(j)
      zs = -zc
      do 821 m=1,mmax
      ans2 = 0._R8
      call gf(ineg,nmult,r2(m),z2(m),rc,zc,ans1)
      if(isym.eq.0 .or. j.eq.2) go to 822
      call gf(ineg,nmult,r2(m),z2(m),rc,zs,ans2)
  822 gfun4(m,mm) = (ans1+ans2)/(tpi*rc*deex)
  821 continue
  820 continue
!
!.....top boundary
!
      j = nz
      imin = 3
      imax = nx
      do 830 i=imin,imax
      rc = xary(i)
      mm = mm+1
      izer(mm) = i
      jzer(mm) = j
      signzer(mm) = -deex
      zc = .5_R8*(zary(nzp)+zary(nz))
      zs =-zc
      do 831 m=1,mmax
      ans2 = 0._R8
      call gf(ineg,nmult,r2(m),z2(m),rc,zc,ans1)
      if(isym.eq.0) go to 832
      call gf(ineg,nmult,r2(m),z2(m),rc,zs,ans2)
  832 gfun4(m,mm) = (ans1+ans2)/(tpi*rc*deez)
  831 continue
  830 continue
      if(isym.eq.1) go to 829
!
!.....bottom boundary
      j = 3
      imin = 3
      imax = nx
      do 840 i=imin,imax
      rc = xary(i)
      mm = mm+1
      izer(mm) = i
      jzer(mm) = j
      signzer(mm) = -deex
      zc = .5_R8*(zary(3)+zary(2))
      do 841 m=1,mmax
      ans2 = 0._R8
      call gf(ineg,nmult,r2(m),z2(m),rc,zc,ans1)
      gfun4(m,mm) = ans1/(tpi*rc*deez)
  841 continue
  840 continue
!
  829 continue
      mmaxmm = mm
   10 continue
!
      if(iflux .eq. 3) go to 31
      if(iflux .eq. 4) go to 32
!...........................................................
!
!...derivatives of green's function needed for plasma
!...contribution to psi from moments method
!
!...........................................................
      m = 0
!
!.....bottom-top boundaries
      do 17 ider=1,2
      do 17 l=1+isym,2
      if(l.eq.1) jb = 2
      if(l.eq.2) jb = nzp
      do 18 i=2,nxp
      m = m+1
      r1(m) = xcurf
      z1(m) = zcurf
      r2(m) = xary(i)
      z2(m) = zary(jb)
      if(ider.eq.2) z2(m) = zary(jb) - (2*l-3)*deez
   18 continue
   17 continue
!
!.....left-right boundaries
      do 19 ider=1,2
      do 19 l=1,2
      if(l.eq.1) ib = 2
      if(l.eq.2) ib = nxp
      do 20 j=2,nzp
      m = m+1
      r1(m) = xcurf
      z1(m) = zcurf
      r2(m) = xary(ib)
      if(ider.eq.2) r2(m) = xary(ib) - (2*l-3)*deex
      z2(m) = zary(j)
   20 continue
   19 continue
!
!.....derivatives
      call gvect(r1,z1,r2,z2,m,g1,g2,g3,g4,g5,g6,nmult,ineg)
!...........................................................
!
!...plasma contribution
!
!......note that tcurdtp is zero for iflux=0
!                cmom is zero for iflux=0 or 1
!
!...........................................................
      m = 0
!
!.....bottom-top boundaries
      do 21 ider=1,2
      do 21 l=1+isym,2
      if(l.eq.1) jb = 2
      if(l.eq.2) jb = nzp
      do 22 i=2,nxp
      m = m+1
      plgf1(m) = g1(m)*tcurdtp
      plgf2(m) = .5_R8*(g4(m)*(cmom(1,2)+cmom(2,1))                      &  
     &              +g5(m)*cmom(2,2) + g6(m)*cmom(1,1))
   22 continue
   21 continue
!
!.....left-right boundaries
      do 23 ider=1,2
      do 23 l=1,2
      if(l.eq.1) ib = 2
      if(l.eq.2) ib = nxp
      do 24 j=2,nzp
      m = m+1
      plgf1(m) = g1(m)*tcurdtp
      plgf2(m) = .5_R8*(g4(m)*(cmom(1,2)+cmom(2,1))                      &  
     &              +g5(m)*cmom(2,2) + g6(m)*cmom(1,1))
   24 continue
   23 continue
!
!...full volume integral method
!
      if(iflux .ne. 3) go to 29
   31 continue
      m = 0
!
!.....bottom-top boundaries
      do 25 ider=1,2
      do 25 l=1+isym,2
      if(l.eq.1) jb = 2
      if(l.eq.2) jb = nzp
      do 26 i=2,nxp
      m = m+1
      xb     = xary(i)
      zb     = zary(jb)
      if(ider.eq.2) zb = zary(jb) - (2*l-3)*deez
      call plfluxb(xb,zb,ans)
      plgf1(m) = ans
      plgf2(m) = 0._R8
   26 continue
   25 continue
!
!.....left-right boundaries
      do 27 ider=1,2
      do 27 l=1,2
      if(l.eq.1) ib = 2
      if(l.eq.2) ib = nxp
      do 28 j=2,nzp
      m = m+1
      xb     = xary(ib)
      zb     = zary(j)
      if(ider.eq.2) xb = xary(ib)-(2*l-3)*deex
      call plfluxb(xb,zb,ans)
      plgf1(m) = ans
      plgf2(m) = 0._R8
   28 continue
   27 continue
!
!
      go to 29
!
   32 continue
!
!
!
!------>>>  add coding here to implement Von-Hagenow method
!
      m = 0
!
!.....bottom-top boundaries
      do 117 ider=1,2
      do 117 l=1+isym,2
      if(l.eq.1) jb = 2
      if(l.eq.2) jb = nzp
      do 118 i=2,nxp
      m = m + 1
      sum = 0._R8
      do 121 mm=1,mmaxmm
      sum = sum + gfun4(m,mm)*psizer(izer(mm),jzer(mm))*signzer(mm)
  121 continue
      iz = i
      jz = jb - (2*l-3)
      plgf1(m) = sum + psizer(iz,jz)*(ider-1)
      plgf2(m) = 0._R8
  118 continue
  117 continue
!
!.....left-right boundaries
      do 119 ider=1,2
      do 119 l=1,2
      if(l.eq.1) ib = 2
      if(l.eq.2) ib = nxp
      do 120 j=2,nzp
      m = m + 1
      sum = 0._R8
      do 122 mm=1,mmaxmm
      sum = sum + gfun4(m,mm)*psizer(izer(mm),jzer(mm))*signzer(mm)
  122 continue
      iz = ib - (2*l-3)
      jz = j
      plgf1(m) = sum + psizer(iz,jz)*(ider-1)
      plgf2(m) = 0._R8
  120 continue
  119 continue
!
!
!
   29 continue
!...................................................................
!
!...discrete coils contribution
!
!....................................................................
      m = 0
!
!.....bottom-top boundaries
      do 53 ider=1,2
      do 53 l=1+isym,2
      if (l.eq.1) jb=2
      if (l.eq.2) jb=nzp
      do 51 i=2,nxp
   51 sumi(i) = 0._R8
      do 30 n=1,ncoilte
!cj dir$ ivdep
      do 30 i=2,nxp
      sumi(i) = sumi(i) + (ans1v(i,l,n,ider)+isym*(sns1v(i,l,n,ider)))   &  
     &                    *ccoil(n)/tpi
   30 continue
      do 50 i=2,nxp
      m = m+1
      if(ider.eq.2) go to 49
      if(l.eq.1) psibb(i) = sumi(i)+plgf1(m)+plgf2(m)
      if(l.eq.2) psibt(i) = sumi(i)+plgf1(m)+plgf2(m)
      go to 50
   49 continue
      if(l.eq.1) psitb(i) = sumi(i)+plgf1(m)+plgf2(m)
      if(l.eq.2) psitt(i) = sumi(i)+plgf1(m)+plgf2(m)
   50 continue
   53 continue
!
!.....left-right boundaries
      do 70 ider=1,2
      do 70 l=1,2
      if(l.eq.1) ib=2
      if(l.eq.2) ib=nxp
      do 71 j=2,nzp
   71 sumj(j) = 0._R8
      do 60 n=1,ncoilte
!cj dir$ ivdep
      do 60 j=2,nzp
      sumj(j) = sumj(j) + (ans2v(j,l,n,ider)+isym*sns2v(j,l,n,ider))     &  
     &                   *ccoil(n)/tpi
   60 continue
      do 70 j=2,nzp
      m = m+1
      if(ider.eq.2) go to 69
      if(l.eq.1) psibl(j) = sumj(j)+plgf1(m)+plgf2(m)
      if(l.eq.2) psibr(j) = sumj(j)+plgf1(m)+plgf2(m)
      go to 70
   69 continue
      if(l.eq.1) psitl(j) = sumj(j)+plgf1(m)+plgf2(m)
      if(l.eq.2) psitr(j) = sumj(j)+plgf1(m)+plgf2(m)
   70 continue
      do 73 j=2,nzp
      psisr(j) = psibr(j)
   73 psisl(j) = psibl(j)
      do 72 i=2,nxp
      psist(i) = psibt(i)
   72 psisb(i) = psibb(i)
!
      return
    2 continue
!***********************************************************
!
!....compute boundary condition for velocity stream function "abig"
!
!***********************************************************
!
!.....=> calculate boundary conditions on velocity to be compatible
!        with changing magnetic fields and ideal MHD
!
      gdot = 0._R8
      do 95 l=1,ntpts
   95 gdot = gdot + facd(l)*gzerov(l)
      efield = gdot*2._R8*alz*log(alx/ccon)/(2._R8*(alx-ccon)+4._R8*alz)     
!
!.....first right side
      sum1 = 0._R8
      sum2 = 0._R8
      sum3 = 0._R8
      indx = 0
      i = nxp
      xsq = xary(i)**2
      gzx = gzero*xary(i)
      do 100 j=3,nzp
      fac = 1._R8
      if(j.eq.nzp) fac = 0.5_R8
      if(j.eq.3 .and. isym.eq.0) fac = 0.5_R8
      indx = indx + 1
      psir = .5_R8*(psi(i,j)-psi(i-1,j)+psi(i,j-1)-psi(i-1,j-1))/deex
      psiz = .5_R8*(psi(i,j)-psi(i,j-1)+psi(i-1,j)-psi(i-1,j-1))/deez
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      psich = 0._R8
      if(dt.ne.0) psich = (psimid-psiold(indx))/dt
      psiold(Indx) = psimid
      psidot = acoef(18)*.5_R8*(psidbr(j)+psidbr(j-1))
      bsq = (psir**2+psiz**2+gzero**2)/xsq
      vn(indx)=(-psir*psidot/xsq-efield*(gzero**2+psiz**2)/gzx)/bsq
      sum3 = sum3 - tpi*deez*fac*psich*psir/xarh(i)
      sum1 = sum1 + xary(i)*vn(indx)*deez
  100 sum2 = sum2 + xary(i)*deez
!
!.....top
      j = nzp
      do 200 ii=3,nxp
      i = nxp+3-ii
      fac = 1._R8
      if(i.eq.nxp.or.i.eq.3) fac = 0.5_R8
      indx = indx + 1
      xsq = xarh(i)**2
      gzx = gzero*xarh(i)
      psir = .5_R8*(psi(i,j)+psi(i,j-1)-psi(i-1,j)-psi(i-1,j-1))/deex
      psiz = .5_R8*(psi(i,j)+psi(i-1,j)-psi(i,j-1)-psi(i-1,j-1))/deez
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      psich = 0._R8
      if(dt.ne.0) psich = (psimid-psiold(indx))/dt
      psiold(indx) = psimid
      psidot = acoef(18)*.5_R8*(psidbt(i)+psidbt(i-1))
      bsq = (psir**2+psiz**2+gzero**2)/xsq
      vn(indx)=(-psiz*psidot/xsq-efield*(gzero**2+psir**2)/gzx)/bsq
      sum3 = sum3 - tpi*deex*fac*psich*psiz/xarh(i)
      sum1 = sum1 + xarh(i)*vn(indx)*deex
  200 sum2 = sum2 + xarh(i)*deex
!
!.....left boundary
      i = 2
      xsq = xary(i)**2
      gzx = gzero*xary(i)
      do 300 jj=3,nzp
      j = nzp+3-jj
      fac = 1._R8
      if(j.eq.nzp) fac = 0.5_R8
      if(j.eq.3 .and. isym.eq.0) fac = 0.5_R8
      indx = indx+1
      psir = .5_R8*(psi(i+1,j)+psi(i+1,j-1)-psi(i,j)-psi(i,j-1))/deex
      psiz = .5_R8*(psi(i+1,j)+psi(i,j)-psi(i+1,j-1)-psi(i,j-1))/deez
      psimid = .25_R8*(psi(i,j)+psi(i+1,j)+psi(i,j-1)+psi(i+1,j-1))
      if(dt.ne.0) psich = (psimid-psiold(indx))/dt
      psiold(indx) = psimid
      psidot = acoef(18)*.5_R8*(psidbl(j)+psidbl(j-1))
      bsq = (psir**2+psiz**2+gzero**2)/xsq
      vn(indx)=( psir*psidot/xsq-efield*(gzero**2+psiz**2)/gzx)/bsq
      sum3 = sum3 + tpi*deez*fac*psich*psir/xarh(3)
      sum1 = sum1 + xary(i)*vn(indx)*deez
  300 sum2 = sum2 + xary(i)*deez
      if(isym.eq.1) go to 500
!
!.....bottom
      j = 2
      do 400 i=3,nxp
      fac = 1._R8
      if(i.eq.nxp.or.i.eq.3) fac = 0.5_R8
      indx = indx + 1
      xsq = xarh(i)**2
      gzx = gzero*xarh(i)
      psir = .5_R8*(psi(i,j)+psi(i,j+1)-psi(i-1,j)-psi(i-1,j+1))/deex
      psiz = .5_R8*(psi(i-1,j+1)+psi(i,j+1)-psi(i-1,j)-psi(i,j))/deez
      psimid=.25_R8*(psi(i,j)+psi(i,j+1)+psi(i-1,j)+psi(i-1,j+1))
      psich = 0._R8
      if(dt.ne.0) psich = (psimid-psiold(indx))/dt
      psiold(indx) = psimid
      psidot = acoef(18)*.5_R8*(psidbb(i)+psidbb(i-1))
      bsq = (psir**2+psiz**2+gzero**2)/xsq
      vn(indx)=( psiz*psidot/xsq-efield*(gzero**2+psir**2)/gzx)/bsq
      sum3 = sum3 + tpi*deex*fac*psich*psiz/xarh(i)
      sum1 = sum1 + xarh(i)*vn(indx)*deex
  400 sum2 = sum2 + xarh(i)*deex
  500 continue
      enerpb = (1+isym)*sum3*udsp/udst
!
      vcom = sum1/sum2
!
!.....decompose into boundary conditions for velocity
!     stream function and potential
      acon = 0._R8
      if(isym.ne.0) go to 99
      jndx = 0
      i = nxp
      do 98 j=3,nh
      jndx = jndx + 1
      acon = acon + xary(i)*(vn(jndx)-vcom)*deez
   98 continue
   99 continue
!
      avbndr(2) = -acon
      jndx = 0
      i = nxp
      do 101 j=3,nzp
      jndx = jndx + 1
      omderr(j) = vcom*deex
  101 avbndr(j) = avbndr(j-1) + xary(i)*(vn(jndx)-vcom)*deez
      j = nzp
      avbndt(nxp) = avbndr(nzp)
      do 201 ii=3,nxp
      i = nxp+3-ii
      jndx = jndx+1
      omdert(i) = vcom*deez
  201 avbndt(i-1) = avbndt(i) + xarh(i)*(vn(jndx)-vcom)*deex
      i = 2
      avbndl(nzp) = avbndt(2)
      do 301 jj=3,nzp
      j = nzp+3-jj
      jndx = jndx+1
      omderl(j) = -vcom*deex
  301 avbndl(j-1) = avbndl(j) + xary(i)*(vn(jndx)-vcom)*deez
      if(isym.eq.1) go to 501
      avbndb(2) = avbndl(2)
      j = 2
      do 401 i=3,nxp
      jndx = jndx + 1
      omderb(i) = -vcom*deez
  401 avbndb(i) = avbndb(i-1) + xarh(i)*(vn(jndx)-vcom)*deex
      if(abs(avbndb(nxp)-avbndr(2)).lt.1.E-6_R8) go to 504
      write(nout,1504) avbndb(nxp)
 1504 format(" error...avbndb(nxp)=",1pe12.4)
      ineg=23
      go to 504
  501 continue
      if(abs(avbndl(2)).lt.1.E-6_R8) go to 504
      write(nout,1501) itype,avbndl(2)
 1501 format(" error...itype,avbndl(2)=",i5,1pe12.4)
      write(nout,1503) sum1,sum2,vcom,avbndl(nzp),efield,gdot,ccon
 1503 format("sum1,sum2,vcom,avbndl(nzp),efield,gdot,ccon",1p7e12.4)
      ineg=23
  504 continue
      if(jndx.eq.indx) go to 502
      write(nout,1502) jndx,indx
 1502 format(" error...jndx and indx do not match ",2i5)
      ineg=23
  502 continue
!
      do 161 i=2,nxp
      psibb(i) = .5_R8*(avbndb(i)+psibbo(i))
      psibt(i) = .5_R8*(avbndt(i)+psibto(i))
      psibbo(i) = psibb(i)
  161 psibto(i) = psibt(i)
      do 162 j=2,nzp
      psibl(j) = .5_R8*(avbndl(j)+psiblo(j))
      psibr(j) = .5_R8*(avbndr(j)+psibro(j))
      psiblo(j) = psibl(j)
  162 psibro(j) = psibr(j)
!
!
      return
!***********************************************************
!
!...unused
!
!***********************************************************
    3 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
