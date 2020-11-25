!#include "f77_dcomplx.h"
       subroutine setup
!
!.....generate grid and plot grid lines
!
      USE CLINAM
      USE SCADVAN
      USE SCR15
      USE SCR16
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,jflip,i,n,ii,nindx,icheck,iw,jw,il,nl,igr,iabs,l,n1
      INTEGER n2,i2,ig,iig,mm,m,ll,isum,iplas,jplas,jj,nhof1,nhof2
      INTEGER ilct
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 fac4d,xl,xe,xx,zz,x1,x2,z1,z2,aturnsv,atnvws,atnvw2s
      REAL*8 atnvw3s,atnvw4s,facl,de,dis,dmid,aone,btmax
      REAL*8 bpmax,btave,dt1s,dt2s,dt3s,t1,t2,t3,t3d,t4,dplateh
      REAL*8 aindcs,ailr,timelr,densl,xmin,zmin,xmax,zmax,xave
      REAL*8 xdiff,zdiff,zfac,x3,z3,x4,z4,z1neg,zlneg,xhof,zhof,dx
      REAL*8 xp,ratio,alams,aindx
      REAL*8 AREAL, sum
!============
      dimension fac4d(4)
!============      
!
!
!.....compute number of zones in each region
!
      nzp = nz+1
      if(isym.eq.0) nzp=2*(int( AREAL(nz+1)/2._R8+ .1_R8))
      nz=nzp-1
      nzm=nz-1
      nxm=nx-1
      nxp=nx+1
      nh=(nzp+2)/2
      if(isym.eq.1) nh=2
!
!.......................................................................
!.....part 1:  generate coordinates of vertices, region by region
!.......................................................................
      deez = (2._R8-isym)*alz/nzm
      dzisq = 1._R8/deez**2
      deex =    (alx - ccon)/nxm
      dxisq = 1._R8/deex**2
      dxdz = deex*deez
      xl = ccon - deex
      do 105 j = nh-isym,nz+2
      zary(j) =      (j-nh)*deez
  105 continue
      if(isym.eq.1) go to 111
      do 110 j = 1,nh-1
      jflip = nz + 3 - j
      zary(j) = -zary(jflip)
  110 continue
  111 continue
!
      do 120 i = 1,nx+2
      xary(i) = xl + (i-1)*deex
            if(i.gt.1) xe = .5_R8*(xary(i)+xary(i-1))
      if(i.eq.1) xe=xary(1)-.5_R8*deex
      ajey(i) = xe*dxdz
      xsqoj(i) = xe/dxdz
      xarh(i) = xe
  120 continue
!
!
!.......................................................................
!.....part 2:  calculate indices for wires, separatrix and limiter
!.......................................................................
      iminsep = int((x1sep-ccon)/deex+2.5_R8)
      jminsep = int((z1sep-zzero)/deez+2.5_R8)
      imaxsep = int((x2sep-ccon)/deex+2.5_R8)
      jmaxsep = int((z2sep-zzero)/deez+2.5_R8)
      iminn = int((xlim-ccon)/deex+2.5_R8)
      jmaxx = int((zlim-zzero)/deez+2.5_R8)
      jminn = 2*nh-jmaxx
      if(isym.eq.1) jminn = 2
      imaxx = int((xlim2-ccon)/deex+2.5_R8)
      if(iminsep.lt.2) iminsep = 2
      if(jminsep.lt.2) jminsep = 2
      if(imaxsep.gt.nxp) imaxsep = nxp
      if(jmaxsep.gt.nzp) jmaxsep = nzp
      if(iminn.lt.4) iminn = 4
      if(jmaxx.gt.nz-1) jmaxx = nz-1
      if(imaxx.gt.nx-1) imaxx = nx-1
!
!.....calculate indices for observation points,wires,and limiters
!
      do 990 n=1,nobs
      xx = xobs(n)
      zz = zobs(n)
!     if(xx.lt.ccon.or.xx.gt.alx.or.abs(zz).gt.alz) ineg=7
      if(xx.lt.ccon.or.xx.gt.alx.or.abs(zz).gt.alz) then
      write(*,1990) n,nobs,xx,zz,ccon,alx,alz
 1990 format(" nn,nobs,xx,zz,ccon,alx,alz =",                           &
     &       /,2i3,1p5e12.4 )
      ineg=7
      endif

      iobs(n) = (xx-ccon)/deex+2.0_R8
      jobs(n) = (zz-zzero)/deez+2.0_R8
  990 continue
!
      do 991 ii=1,nwire
      n=ncoil-nwire+ii
      xx = xcoil(n)
      zz = zcoil(n)
      if(xx.lt.ccon.or.xx.gt.alx.or.abs(zz).gt.alz) ineg=8
      iwire(ii) = int((xx-ccon)/deex + 2.5_R8)
      jwire(ii) = int((zz-zzero)/deez + 2.5_R8)
      if(iwire(ii).le.3 .or. iwire(ii).ge.nx) ineg=8
      if(jwire(ii).lt.2 .or. jwire(ii).ge.nz) ineg=8
      if(jwire(ii).le.3 .and. isym.eq.0) ineg=8
  991 continue
!
!.....check for coil inside grid
      if(ncoil.le.nwire) go to 892
      do 891 n=1,ncoil-nwire
      xx = xcoil(n)
      zz = zcoil(n)
      if(xx.lt.ccon .or. xx.gt.alx .or. abs(zz).gt.alz) go to 891
      write(nout,1891) n,xx,zz
 1891 format(" coil",i4," with coordinates",2f6.2," is inside grid")
      ineg=26
  891 continue
  892 continue
!
!
!......check for neg group numbers
!......signifying distributed coil currents
      nindx = nwire
      do 997 ii=1,nwire
      n = ncoil-nwire+ii
      if(igroupw(ii).ge.0) go to 997
!
      icheck = (ncoil-nwire) + nindx + 3
      if(icheck.lt.pncoil) go to 982
      ineg=6
      write(nout,1982)
 1982 format("  new coils from neg group numbers exceed param size")
      return
  982 continue
!
      xx = xcoil(n)
      zz = zcoil(n)
      iw = (xx-ccon)/deex + 2
      jw = (zz+alz*(1-isym))/deez + 2
      x1 = xary(iw)
      x2 = xary(iw+1)
      z1 = zary(jw)
      z2 = zary(jw+1)
      aturnsv = aturnsw(ii)
      iwire(ii) = iw
      jwire(ii) = jw
      iwire(nindx+1) = iw
      jwire(nindx+1) = jw+1
      iwire(nindx+2) = iw+1
      jwire(nindx+2) = jw+1
      iwire(nindx+3) = iw+1
      jwire(nindx+3) = jw
      fac4d(1) = (x2-xx)*(z2-zz)/dxdz
      fac4d(2) = (x2-xx)*(zz-z1)/dxdz
      fac4d(3) = (xx-x1)*(zz-z1)/dxdz
      fac4d(4) = (xx-x1)*(z2-zz)/dxdz
      aturnsw(ii) = aturnsv*fac4d(1)
      do 986 i=1,3
      aturnsw(nindx+i) = aturnsv*fac4d(i+1)
  986 continue
      do 998 i=1,4
      il = nindx + (i-1)
      if(i.eq.1) il = ii
      nl = ncoil - nwire + il
      xwire(il) = xary(iwire(il))
      zwire(il) = zary(jwire(il))
      igroupw(il) = igroupw(ii)
      rswires(il) = rswires(ii)
      rswire(il) = rswire(ii)
      aindw(il) = aindw(ii)
      igr = iabs(igroupw(il))
      cwire0(il) = gcur(istart,igr)*aturnsw(il)*usdi
      ilwire(il) = ilwire(ii)
      ccoil(nl) = cwire0(il)
      ccoils(nl) = ccoil(nl)*udsi
      xcoil(nl) = xary(iwire(il))
      zcoil(nl) = zary(jwire(il))
      aturnsc(nl) = aturnsw(il)
      igroupc(nl) = igroupw(il)
      do 988 l=1,ntpts
      if(i.gt.1) go to 987
      atnvws = atnvw1(l,ii)
      atnvw2s = atnvw2(l,ii)
      atnvw3s = atnvw3(l,ii)
      atnvw4s = atnvw4(l,ii)
  987 continue
      atnvw1(l,il) = atnvws*fac4d(i)
      atnvw2(l,il) = atnvw2s*fac4d(i)
      atnvw3(l,il) = atnvw3s*fac4d(i)
      atnvw4(l,il) = atnvw4s*fac4d(i)
  988 continue
      ngrvw1(il) = ngrvw1(ii)
      ngrvw2(il) = ngrvw2(ii)
      ngrvw3(il) = ngrvw3(ii)
      ngrvw4(il) = ngrvw4(ii)
  998 continue
      nindx = nindx + 3
  997 continue
      ncoil = ncoil + (nindx-nwire)
      nwire = nindx
      if(nwire.gt.pnwire) ineg=6
      if(ncoil.gt.pncoil) ineg=6
      if(ineg.eq.0) go to 999
      n1 = pnwire
      n2 = pncoil
      write(nout,1999) nwire,ncoil,n1,n2
 1999 format("  nwire,ncoil,pnwire,pncoil = ",4i5)
  999 continue
!
      do 995 ii=1,nwire-1
      iw = iwire(ii)
      jw = jwire(ii)
      do 996 i2 = ii+1,nwire
      if(iw.ne.iwire(i2) .or. jw.ne.jwire(i2)) go to 996
      ineg   = 3
      write(nout,1992) ii,i2,iw,jw,xwire(i2),zwire(i2)
 1992 format("  * * * error * * *   wires no ",2i4," have i,j = ",       &  
     &  2i3," x,z =",1p2e12.4)
      write(nout,1993) xwire(ii),zwire(ii),xwire(i2),zwire(i2)
 1993 format(10x," x-z coordinates  ",1p2e12.4,5x,1p2e12.4,/ )
  996 continue
  995 continue
!
      do 992 ii=1,nlim
      xx = xlima(ii)
      zz = zlima(ii)
      if(xx.lt.ccon.or.xx.gt.alx.or.abs(zz).gt.alz) ineg=9
      ilima(ii) = int((xx-ccon)/deex +2.5_R8)
      jlima(ii) = int((zz-zzero)/deez +2.5_R8)
  992 continue
!
      if(iplate.le.0) go to 451
      do 456 n=1,nplate
      if(nseg(n).le.0) nseg(n) = pnseg
      do 452 l=1,nseg(n)+1
      if(xsega(n,l).le.0) go to 450
  452 continue
      go to 453
  450 continue
!
      do 454 l=1,nseg(n)+1
      facl = (l-1._R8)/nseg(n)
      xsega(n,l) = xlplate(n)+(xrplate(n)-xlplate(n))*facl
      zsega(n,l) = zlplate(n)+(zrplate(n)-zlplate(n))*facl
  454 continue
  453 continue
      de = 0._R8
      do 455 l=1,nseg(n)
      dis =                                                              &  
     & sqrt((xsega(n,l+1)-xsega(n,l))**2+(zsega(n,l+1)-zsega(n,l))**2)
      dmid = de + 0.5_R8*dis
      de = de + dis
      dplate(n,l) = dmid
  455 continue
  456 continue
  451 continue
!
!.....define resgs,ngroup,nogroup,ngroups,nogroups for convenience
!
!      ngroup is the total number of wire groups used
!      ngroups is the number of wire groups with series connections
!      ngroupt is the total number of coil (including wire) groups used
!
!      nogroup(n);n=1,ngroup are the wire group numbers
!      nogroups(n),n=1,ngroups are the series wire group numbers
!      nogroupt(n),n=1,ngroupt are the coil group numbers
!
!      sgroup(ig),ig=1,pngroup are the number of coils in each group
!
      do 870 ig=1,pngroup
      if(resgs(ig).gt.0) go to 872
  870 continue
      iresgs=0
      go to 873
  872 continue
      iresgs = 1
  873 continue
!
      ngroup = 0
      ngroups= 0
      if(nwire.eq.0) go to 876
      do 874 ii=1,nwire
      ig = iabs(igroupw(ii))
      if(ngroup.eq.0) go to 877
      do 875 n=1,ngroup
      if(ig.eq.nogroup(n)) go to 874
  875 continue
  877 continue
      ngroup = ngroup + 1
      nogroup(ngroup) = ig
      if(iseries(ig).eq.0) go to 874
      ngroups = ngroups+1
      nogroups(ngroups) = ig
  874 continue
      if(ngroups .gt. pngrps) ineg=6
      if(ineg .eq. 6)                                                    &  
     & print *, " in setup: ngroups = ", ngroups, " pngrps = ", pngrps
      if(ineg.ne.0) return
  876 continue
  973 continue
!
      do 864 iig=1,ngroups
      ig = nogroups(iig)
      icmaxs(iig) = 0
      do 865 ii=1,nwire
      if(ig .ne. iabs(igroupw(ii))) go to 865
      icmaxs(iig) = icmaxs(iig)+1
      icoils(icmaxs(iig),iig) = ii
  865 continue
  864 continue
!
      ngroupt=0
      if(ncoil.eq.0) go to 976
      do 974 ii=1,ncoil
      ig = iabs(igroupc(ii))
      if(ngroupt.eq.0) go to 977
      do 975 n=1,ngroupt
      if(ig.eq.nogroupt(n)) go to 974
  975 continue
  977 continue
      ngroupt = ngroupt + 1
      nogroupt(ngroupt) = ig
  974 continue
  976 continue
      if(numfb.eq.0) go to 1976
      do 1974 ii=1,numfb
      ig = nrfb(ii)
      if(ig.eq.0) go to 1974
      if(ngroupt.eq.0) go to 1977
      do 1975 n=1,ngroupt
      if(ig.eq.nogroupt(n)) go to 1974
 1975 continue
 1977 continue
      ngroupt = ngroupt + 1
      nogroupt(ngroupt) = ig
!     write(nterm,2222) ngroupt,ig
!2222 format("ngroup increased to",2i3)
 1974 continue
 1976 continue
!
      do 10 mm=1,ngroup
      m = nogroup(mm)
   10 sgroup(m) = 0._R8
      if(isym.eq.1) go to 25
      do 20 ii=1,nwire
      ig = iabs(igroupw(ii))
   20 sgroup(ig) = sgroup(ig)+ 1._R8
      go to 27
   25 continue
      do 26 ii=1,nwire
      aone = 1._R8
      if(zwire(ii).eq.0 ) aone=0.5_R8
      ig = iabs(igroupw(ii))
   26 sgroup(ig) = sgroup(ig)+ aone
   27 continue
!
!
!.......................................................................
!.....part 3:  printout setup info
!.......................................................................
!
      entry prset
      btmax = gzero/ccon
      aminor = axic
      if(axic.le.0) aminor = .25_R8*(alx-ccon)
      bpmax = 1._R8
      if(tcuro.gt.0) bpmax = tcuro/(tpi*aminor)
      btave = gzero/(.5_R8*(xary(3)+xary(nx)))
      if(gzero.le.0) btmax = tcuro/(tpi*xplas)
      if(btmax.le.0) btmax = 1._R8
      if(gzero.le.0) btave = btmax
      dt1s = dtfac*ndiv*min(deex,deez)**2/(etav*8._R8*(1._R8-th))*udst
      dt2s = dtfac*min(deex,deez)/(sqrt(2._R8)*bpmax)*udst
      dt3s = dtfac*ndiv*min(deex,deez)/(2.8_R8*btmax)*udst
      t1 = (alx-ccon)*2._R8*alz/etav*udst
      t2 = sqrt((alx-ccon)*2._R8*alz)/btave*udst
      amux = amu*.003_R8*sqrt((alx-ccon)*alz)
      t3 = .5_R8*(alx-ccon)*alz/amux*udst
      t3d = 0._R8
      if(acoef(90) .gt. 0._R8)  t3d = udst / acoef(90)
      t4 = t3/acoef(10)
      write(nout,1400)
!
      write(nout,1010) udsd,udst,udsv,udsi,udsr,udse
      write(nout,1120)
 1120 format(//," xary...(2) to (nx+1)")
      write(nout,1121) (xary(i),i=2,nx+1)
 1121 format(1x,1p10e12.4)
      write(nout,1122)
 1122 format(//," zary...(2) to (nz+1)")
      write(nout,1121) (zary(i),i=2,nz+1)
      write(nout,1000) nx,nz,nh
 1000 format(//,1x,"     geometry parameters and switches",//,1x         &  
     &       ," nx=",i3," nz=",i3," nh=",i3)
!
      write(nout,1420) irst1,irst2,ipest,isym,idata,lrswtch,             &  
     &idens,ipres ,ifunc,iflux,icube,isurf,itrmod,idiv,iwall,            &  
     &ialpha,iplate,ibalsw,icirc,isvd,isaw,irfp,itevv,irippl,            &  
     &iimp,ilte,ibootst,iffac
 1420 format(/,"    irst1   ",i2,"    irst2   ",i2,"    ipest   "        &  
     & ,i2,  "    isym    ",i2,"    idata   ",i2,"    lrswtch "          &  
     & ,i2,/,"    idens   ",i2,"    ipres   ",i2,"    ifunc   "          &  
     & ,i2,  "    iflux   ",i2,"    icube   ",i2,"    isurf   "          &  
     & ,i2,/,"    itrmod  ",i2,"    idiv    ",i2,"    iwall   "          &  
     & ,i2,  "    ialpha  ",i2,"    iplate  ",i2,"    ibalsw  "          &  
     & ,i2,/,"    icirc   ",i2,"    isvd    ",i2,"    isaw   "           &  
     & ,i2,  "    irfp    ",i2,"    itevv   ",i2,"    irippl  "          &  
     & ,i2,/,"    iimp    ",i2,"    ilte    ",i2,"    ibootst "          &  
     & ,i2,  "    iffac   ",i2)
!
      write(nout,1400)
 1400 format(/)
      write(nout,1414) 1000._R8*dtmins
      write(nout,1415) 1000._R8*dtmaxs
      write(nout,1416) dtfac
      write(nout,1417) ndiv
      write(nout,1418) ffac
      write(nout,1419) tevv
      write(nout,1400)
      write(nout,1401) 1000._R8*dt1s
      write(nout,1402) 1000._R8*dt2s
      write(nout,1412) 1000._R8*dt3s
      write(nout,1403) 1000._R8*t1
      write(nout,1404) 1000._R8*t2
      write(nout,1405) 1000._R8*t3
      write(nout,1413) 1000._R8*t3d
      write(nout,1406) 1000._R8*t4
 1401 format(" time step based on vacuum resistivity ",e13.5,"(msec)")
 1402 format(" time step based on alfven wave        ",e13.5,"(msec)")
 1414 format(" min time step (input).................",e13.5,"(msec)")
 1415 format(" max time step (input).................",e13.5,"(msec)")
 1416 format(" dtfac--time step safety factor               ",e13.5)
 1417 format(" ndiv--(sub-cycles for diff and fast wave  )  ",i10 )
 1418 format(" ffac--(factor alfven and fast waves slowed ) ",f13.0)
 1419 format(" tevv--(vacuum temp for resistivity calc      ",e13.5)
 1412 format(" time step based on fast wave          ",e13.5,"(msec)")
 1413 format(" diffusion transit time for drag       ",e13.5,"(msec)")
 1403 format(" diffusion transit time                ",e13.5,"(msec)")
 1404 format(" wave transit time                     ",e13.5,"(msec)")
!
 1405 format(" compressible viscous damping time     ",e13.5,"(msec)")
 1406 format(" incompressible viscous damping time   ",e13.5,"(msec)")
!
      write(nout,1805) iminsep,jminsep,imaxsep,jmaxsep,                  &  
     &          iminn,imaxx,jmaxx
 1805 format(/,"  iminsep",i3,5x,"jminsep",i3,/,                         &  
     &         "  imaxsep",i3,5x,"jmaxsep",i3,/,                         &  
     &         "  iminn",i3,5x,"imaxx",i3,5x,"jmaxx",i3)
!
      if(iplate.le.0) go to 1501
      write(nout,1502)
 1502 format(///," endpoints of divertor plates ",/,                     &  
     &      "  n    xl(n)   zl(n)   xr(n)   zr(n)   fraction")
      do 1503 n=1,nplate
 1503 write(nout,1504) n,xlplate(n),zlplate(n),xrplate(n),zrplate(n),    &  
     &                 fplate(n,1),fplate(n,2)
 1504 format(i5,6f8.3)
!
      do 1551 n=1,nplate
      write(nout,1561) n
 1561 format("1   * * * divertor plate ",i3," * * *",//,                 &  
     &"    i  xsega(n,i) zsega(n,i) dplate(n,i)")
      do 1571 l=1,nseg(n)+1
      dplateh = 0._R8
      if(l.gt.1) dplateh = dplate(n,l-1)
      write(nout,2571) l,xsega(n,l),zsega(n,l),dplateh
 2571 format(1x,i4,1p3e12.4)
 1571 continue
!
 1551 continue
 1501 continue
      if(nobs.le.0) go to 2807
!
!......print observation point info
      write(nout,1806)
 1806 format(//,"  observation points ",/,                               &  
     & "        i1   j1   i2   j2  nplot  ",                             &  
     &          "   x1          z1          x2          z2")
      do 807 n=1,nobs,2
      nn = (n+1)/2
  807 write(nout,1807)nn,iobs(n),jobs(n),iobs(n+1),jobs(n+1),npltobs(n)  &  
     &  ,xobs(n),zobs(n),xobs(n+1),zobs(n+1)
 1807 format(6i5,1p4e12.4)
      if(ncoil-nwire .le. 0) go to 3808
!
!.....print external coil info
      write(nout,4808)
 4808 format(//" external coil information:",//,                         &  
     &"   n          x          z group      turns",                     &  
     &"    rscoils     aindcs         dx         dz",                    &  
     &" fcu fss      tempc    l/rtime")
      do 3809 n=1,ncoil-nwire
      aindcs = aindc(n)/udsi
      if(rscoils(n).ne.0) ailr = aindcs/rscoils(n)
      write(nout,4809) n,xcoil(n),zcoil(n),igroupc(n),aturnsc(n),        &  
     &           rscoils(n),aindcs,dxcoil(n),dzcoil(n),fcu(n),fss(n),    &  
     &                 tempc(n),ailr
 4809 format(i4,1p2e12.4,i5,1p5e11.3,0p2f5.2,1p2e11.3)
 3809 continue
 3808 continue
!
 2807 if(nwire.le.0) go to 2808
!
!......print wire info
      write(nout,1808)
 1808 format(//" internal coil information:",                            &  
     &       /,"        iwire   jwire igroupw aturnsw   xwire       ",   &  
     &    "zwire       rswires     aindw       cwics      lr-time")
      do 809 n=1,nwire
      timelr = aindw(n)/rswires(n)
      write(nout,1809) n,iwire(n),jwire(n),igroupw(n),aturnsw(n),        &  
     &      xwire(n),zwire(n),rswires(n),aindw(n),cwics(n),timelr
  809 continue
 1809 format(i5,3i8,f8.1,1p6e12.4)
 2808 if(nlim.le.0) go to 2809
!
!......print limiter info
      write(nout,1810)
 1810 format(//"        ilima   jlima   xlima        zlima")
      do 810 n=1,nlim
      write(nout,1811) n,ilima(n),jlima(n),xlima(n),zlima(n)
  810 continue
 1811 format(i5,2i8,1p2e12.4)
 2809 if(numfb.le.0) go to 831
!
!......print feedback info
      write(nout,1820)
      do 822 l=1,numfb
  822 write(nout,1822) l,nrfb(l),ipext(l),nfeedo(l),nfeed2(l),fbfac(l),  &  
     &              fbcon(l),fbfac1(l),fbfacd(l),tfbons(l),              &  
     &              tfbofs(l),fbfaci(l),idelay(l)
 1822 format(5i5,1p7e12.4,i5)
!
      write(nout,1942)
 1942 format(//,"  variable feedback observation points",                &  
     &       //,"  inumfb  nfeedv(ll,inumfb),ll=1,ntpts",/)
      do 843 l=1,numfb
      write(nout,1843) l,(nfeedv(ll,l),ll=1,ntpts)
  843 continue
 1843 format(2x,i5,15i5,10(/,7x,15i5))
!
!......calculate and print normalized feedback gains
      if(ineg.ne.0) return
      call feednorm
!
      isum = 0
      do 922 ii=1,nwire
  922 isum = isum + ngrvw1(ii) + ngrvw2(ii) + ngrvw3(ii) + ngrvw4(ii)
      if(isum.le.0) go to 821
      write(nout,1922)
 1922 format(//,"  variable turn feedback(internal coils):",//,          &  
     &"  ii ngrvw   atnvw(l),l=1,ntpts")
      do 823 ii=1,nwire
      if(ngrvw1(ii).le.0) go to 723
      write(nout,1823) ii,ngrvw1(ii),(atnvw1(l,ii),l=1,ntpts)
  723 continue
      if(ngrvw2(ii).le.0) go to 724
      write(nout,1823) ii,ngrvw2(ii),(atnvw2(l,ii),l=1,ntpts)
  724 continue
      if(ngrvw3(ii).le.0) go to 725
      write(nout,1823) ii,ngrvw3(ii),(atnvw3(l,ii),l=1,ntpts)
  725 continue
      if(ngrvw4(ii).le.0) go to 823
      write(nout,1823) ii,ngrvw4(ii),(atnvw4(l,ii),l=1,ntpts)
  823 continue
 1823 format(2i5,1p8e10.2,10(/,10x,1p8e10.2) )
  821 continue
!
      if(ncoil.eq.nwire) go to 831
      isum = 0
      do 932 ii=1,ncoil-nwire
  932 isum = isum + ngrvc1(ii) + ngrvc2(ii) + ngrvc3(ii) + ngrvc4(ii)
      if(isum.le.0) go to 831
      write(nout,1932)
 1932 format(//,"  variable turn feedback(external coils):",//,          &  
     &"  ii ngrvc   atnvc(l),l=1,ntpts")
      do 833 ii=1,ncoil-nwire
      if(ngrvc1(ii).le.0) go to 733
      write(nout,1833) ii,ngrvc1(ii),(atnvc1(l,ii),l=1,ntpts)
  733 continue
      if(ngrvc2(ii).le.0) go to 734
      write(nout,1833) ii,ngrvc2(ii),(atnvc2(l,ii),l=1,ntpts)
  734 continue
      if(ngrvc3(ii).le.0) go to 735
      write(nout,1833) ii,ngrvc3(ii),(atnvc3(l,ii),l=1,ntpts)
  735 continue
      if(ngrvc4(ii).le.0) go to 833
      write(nout,1833) ii,ngrvc4(ii),(atnvc4(l,ii),l=1,ntpts)
  833 continue
 1833 format(2i5,1p8e10.2,10(/,10x,1p8e10.2) )
  831 continue
 1820 format(//"    l nrfb ipext  nfo  nf2     fbfac       fbcon",       &  
     &       "      fbfac1      fbfacd     tfbon       tfbof  ",         &  
     &       "   fbfaci      idelay")
!
!......time programming
      write(nout,1830)
      do 830 l=1,ntpts
      densl = rnorm(l)*udsd
      write(nout,1862) l,tpros(l),pcur(l),ppres(l),densl,beamp(l),       &
     &    gzerov(l),vloopv(l),zeffv(l),tevv0(l),ffac0(l)
!cj      write(nout,1832)
 1832 format(/)
  830 continue
 1830 format(//,"    l  time        pl cur"                         &  
     &,"      pl press    pl dens     beam power  gzerov"                 &  
     &,"      vloopv      zeffv       tevv0       ffac0" )
 1831 format(1x,i5,1p10e12.4,10(/,90x,1p3e12.4) )
!
!......time programming
      write(nout,1834)
      do 834 l=1,ntpts
      write(nout,1862) l,frcparv(l),plhamp(l),alpharv(l),betarv(l),      &
     &  rzerv(l),azerv(l),ezerv(l),dzerv(l),xmagz(l),zmagz(l)
!cj      write(nout,1832)
  834 continue
 1834 format(//,"    l  fracpar     plhamp"                         &  
     &,"      alpharv     betarv      rzerv       azerv"                 &  
     &,"       ezerv       dzerv       xmagz       zmagz" )
!
!......time programming
      write(nout,1836)
      do 836 l=1,ntpts
      write(nout,1862) l,alhd(l),dlhd(l),a1lhd(l),a2lhd(l),      &  
     &  aclhd(l),dclhd(l),a1clhd(l),a2clhd(l),picrh(l),thalov(l)
  836 continue
 1836 format(//,"    l  alhd        dlhd"                         &  
     &,"        a1lhd       a2lhd       aclhd       dclhd"                 &  
     &,"       a1clhd      a2clhd      picrh       thalov")
!
!......time programming
      write(nout,1838)
      do 838 l=1,ntpts
      write(nout,1862) l,whalov(l),fwcd(l),afwd(l),dfwd(l),      &
     &  a1fwd(l),a2fwd(l),acfwd(l),dcfwd(l),a1cfwd(l),a2cfwd(l)
  838 continue
 1838 format(//,"    l  whalov      fwcd"                         &
     &,"       afwd        dfwd        a1fwd       a2fwd"                 &  
     &,"       acfwd       dcfwd       a1cfwd      a2cfwd")
!
!......time programming
      write(nout,1840)
      do 840 l=1,ntpts
      write(nout,1862) l,heactv(l),rnorm(l),qaddv(l), &   !shmodi(l),      &  
     &  pwidthcv(l),chipedv(l),fhmodeiv(l),tpedv(l),sawtime(l)
  840 continue
 1840 format(//,"    l  heactv      rnorm"                         &  
     &,"       qaddv       pwidthcv    chipedv    fhmodeiv"                 &
     &,"     tpedv       sawtime" )
!
!
!......time programming
      write(nout,1839)
      do 839 l=1,ntpts
      write(nout,1862) l,nflagv(l),expn1v(l),expn2v(l), &   !shmodi(l),      &
     &  firitbv(l),secitbv(l),fracn0v(l),newdenv (l)
  839 continue
 1839 format(//,"    l  nflagv      expn1v"                         &
     &,"       expn2v       firitbv     secitbv    fracn0v    newdenv ")
!
!......time programming
      write(nout,1842)
      do 842 l=1,ntpts
      write(nout,1862) l,fraciv(1,l),fraciv(2,l),fraciv(3,l),fraciv(4,l),      &  
     &  fraciv(5,l),fraciv(6,l),fraciv(7,l),fraciv(8,l)
  842 continue
 1842 format(//,"    l  fraciv(1)   fraciv(2)"                         &  
     &,"   fraciv(3)   fraciv(4)   fraciv(5)"                 &  
     &,"   fraciv(6)   fraciv(7)   fraciv(8)" )
!
!......time programming
      write(nout,1845)
      do l=1,ntpts
         write(nout,1862) l,pecrh(l),eccd(l),aecd(l),decd(l),             &
     &     a1ecd(l),a2ecd(l)
      enddo
 1845 format(//,"    l  pecrh      eccd"                         &
     &,"       aecd        decd        a1ecd       a2ecd")
!
!2/22/10 tpros, pcur, ppres, rnorm, beamp, gzerov, vloopv, zeffv, tevv0, ffac0,
!2/22/10 frcparv, plhamp, alpharv, betarv, rzerv, azerv, ezerv, dzerv, xmagz,
!2/22/10 zmagz, alhd, dlhd, a1lhd, a2lhd, aclhd, dclhd, a1clhd, a2clhd, picrh,
!2/22/10 thalov, whalov, fwcd, afwd, dfwd, a1fwd, a2fwd, acfwd, dcfwd, a1cfwd,
!2/22/10 a2cfwd, heactv, sawtime, qadd, shmodi, pwidthc, chiped, tped, fraci(:8) 
!
!.....print out acoef array
      write(nout,1860)
 1860 format(/,"    l        acoef array",/)
      do 861 l=1,4999-9+pncoil,10
      sum = 0._R8
      do 862 i=1,10
  862 sum = sum + abs(acoef(l+i-1))
      if(sum.eq.0) go to 861
      write(nout,1862) l,(acoef(l+i-1),i=1,10)
 1862 format(i5,1p11e12.4)
  861 continue
!
      if(itemp.gt.0 .and. irst1.ne.1) call temprise
!
!.....initialize arrays to keep current from flowing outside v-vessel
      do 776 i=1,nxp
      jminny(i)=jminn
      jmaxxy(i)=jmaxx
      do 776 j=1,nzp
      iexvc(i,j) = 0
  776 iexv(i,j) = 0
      if(acoef(1).eq.4._R8) call itervv (ivvlo,ivvhi)
      if(acoef(1).eq.3._R8) call d3dvv(jvvlo,jvvhi)
      if(acoef(1).eq.2._R8.and. idata.ne.3) call vvexc(jvvlo,jvvhi)
      if(acoef(1).eq.2._R8.and. idata.eq.3) call vvexc1(jvvlo,jvvhi)
        if (acoef(1).eq.1._R8)   then
                           call  pbxvv (ivvlo, ivvhi)
                                           endif
      if(acoef(1).eq.0 .or. acoef(1).gt.4._R8) go to 1777
      iplas = (xplas-ccon)/deex + 2
      jplas = (zplas-zzero)/deez + 2
      if(iexv(iplas,jplas).eq.0) go to 1777
      ineg=44
      write(nout,1778)
 1778 format("error...masking arrays cover plasma.  This is caused",     &  
     &     /,"by acoef(1) > 0, and by not having enough type 10 cards",  &  
     &     /,"to completely define a vacuum vessel")
 1777 continue
      do 777 i=1,nxp+1
      do 777 j=1,nzp+1
      iforce(i,j) = 0
      if(i.le.iminn .or. i.ge.imaxx .or. j.ge.jmaxx) iexv(i,j)=1
      if(isym.eq.0 .and. j.le.jminn) iexv(i,j)=1
  777 continue
      do 778 j=3,nzp+1
      do 778 i=3,nxp+1
      if(iexv(i  ,j  ) .eq. 1 .and. iexv(i-1,j  ) .eq. 1 .and.           &  
     &   iexv(i  ,j-1) .eq. 1 .and. iexv(i-1,j-1) .eq. 1) iexvc(i,j)=1
  778 continue
      do 774 i=2,nxp+1
      if(isym.eq.1) iexvc(i,2) = iexvc(i,3)
      if(isym.eq.1) iexv (i,1) = iexv (i,3)
      if(isym.eq.0) iexvc(i,3) = 1
      iexvc(i,nzp) = 1
  774 iexvc(i,nzp+1) = 1
      do 779 j=2,nzp+1
      iexvc(3,j) = 1
      iexvc(nxp,j) = 1
  779 iexvc(nxp+1,j) = 1
      do 300 i=2,nxp+1
      do 300 j=2,nzp+1
  300 icoil(i,j) = 0
      do 400 n=1,nwire
  400 icoil(iwire(n),jwire(n))=igroupw(n)
!
      do 781 i=iminn,imaxx
      do 782 jj=jminn,jmaxx
      j = jminn+jmaxx-jj
      if(iexv(i,j).eq.0 .or. iexvc(i,j).eq.0) go to 783
  782 continue
      jmaxxy(i) = jminn-1
      go to 784
  783 jmaxxy(i) =j+1
  784 continue
      if(isym.eq.1) go to 781
      do 785 jj=jminn,jmaxx
      j = jj
      if(iexv(i,j).eq.0 .or. iexvc(i,j).eq.0) go to 786
  785 continue
      jminny(i) = jmaxx+1
      go to 781
  786 jminny(i)=j-1
  781 continue
!
!.....changed 9/4/00 for chi modeling
!....removed setting of iforce on 5/26/10
!....changed back 11/24/11
!
      do 795 j=3,nzp
!     do 793 ii=3,nxp
!     i = ii
!     if(icoil(i-1,j).gt.0 .and. icoil(i-1,j-1).gt.0) go to 791
!     if(iexvc(i,j).eq.0) go to 793
!     iforce(i,j)=1
! 793 continue
  791 continue
!     do 792 ii=3,nxp
!     i = 3+nxp-ii
!     if(icoil(i,j).gt.0 .and. icoil(i,j-1).gt.0) go to 795
!     if(iexvc(i,j).eq.0) go to 795
!     iforce(i,j)=1
! 792 continue
  795 continue
!
!
!......define boundary condition arrays dxa,dzb,dxc,dzd,face
      do 55 j=2,nzp
      do 55 i=2,nxp
      if(iexvc(i,j).eq.0) go to 48
      dxa(i,j) = 0
      dzb(i,j) = 0
      dxc(i,j) = 0
      dzd(i,j) = 0
      face(i,j) = 0._R8
      go to 55
   48 continue
      dxa(i,j) = dxisq
      dzb(i,j) = dzisq
      dxc(i,j) = dxisq
      dzd(i,j) = dzisq
      face(i,j) = 1._R8
   55 continue
      do 197 i=2,nxp
      do 197 j=2,nzp
      if(iexv(i,j).eq.0 .or. iexv(i,j-1).eq.0) go to 194
      dxa(i,j)=0
  194 continue
      if(iexv(i-1,j).eq.0 .or. iexv(i,j).eq.0) go to 195
      dzb(i,j)=0
  195 continue
      if(iexv(i-1,j-1).eq.0 .or. iexv(i-1,j).eq.0) go to 196
      dxc(i,j)=0
  196 continue
      if(iexv(i,j-1).eq.0 .or. iexv(i-1,j-1).eq.0) go to 197
      dzd(i,j)=0
  197 continue
!
      do 551 i=2,nxp
      dzb(i,nzp) = 0
  551 continue
      if(isym.ne.0) go to 553
      do 554 i=2,nxp
      dzd(i,3) = 0
  554 continue
  553 continue
      do 552 j=2,nzp
      dxc(3,j) = 0._R8
  552 dxa(nxp,j) = 0._R8
  555 continue
      call prarray
!
      if((imovie.gt.0 .and. imovie.lt.10)) return
!.......................................................................
!.....part 4:  now plot the grid
!.......................................................................
!
!..rxw/02/03/88
      if(noplot(1).gt.0) goto 1231
!
      xmin = ccon
      zmin = -alz
      xmax = alx
      zmax = alz
!
!.....external coils
      if(itemp.le.0) go to 761
      do 760 n=1,ncoil-nwire
      xmax = max(xmax,xcoil(n)+0.5_R8*dxcoil(n))
      xmin = min(xmin,xcoil(n)-0.5_R8*dxcoil(n))
      zmax = max(zmax,zcoil(n)+0.5_R8*dzcoil(n))
      zmin = min(zmin,-zmax)
  760 continue
  761 continue
      xave = .5_R8*(xmax+xmin)
      xdiff = xmax-xmin
      zdiff = zmax-zmin
      if((zdiff).gt.xdiff) xmax = xave+.5_R8*(zdiff)
      if((zdiff).gt.xdiff) xmin = xave-.5_R8*(zdiff)
      if((xdiff).gt.(zdiff)) zmax = .5_R8*(xdiff)
      if((xdiff).gt.(zdiff)) zmin = -.5_R8*(xdiff)
      call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.285_R8,1._R8)
      do 730 ll=1,isym+1
      zfac = 1._R8
      if(ll.eq.2) zfac = -1._R8
!
      if(itemp.eq.0) go to 790
!
!.....draw external coils
      do 780 n=1,ncoil-nwire
      if(xcoil(n).ge.100._R8) go to 780
      x1 = xcoil(n)-.5_R8*dxcoil(n)
      z1 = zcoil(n)-.5_R8*dzcoil(n)
      x2 = x1+dxcoil(n)
      z2 = z1
      x3 = x2
      z3 = z2+dzcoil(n)
      x4 = x1
      z4 = z3
      call setcrt(x1,z1*zfac)
      call vector(x2,z2*zfac)
      call vector(x3,z3*zfac)
      call vector(x4,z4*zfac)
      call vector(x1,z1*zfac)
  780 continue
  790 continue
!
!
!.....draw vertical lines
!
      do 610 i=1,nx
      xx = xary(i+1)
      zz = zfac*zary(2)
      call setcrt(xx,zz)
      do 620 j=2,nz
      xx = xary(i+1)
      zz = zfac*zary(j+1)
      call vector(xx,zz)
  620 continue
  610 continue
!
!.....draw horizontal lines
!
      do 710 j=1,nz
      xx = xary(2)
      zz = zfac*zary(j+1)
      call setcrt(xx,zz)
      do 720 i=2,nx
      xx = xary(i+1)
      zz = zfac*zary(j+1)
      call vector(xx,zz)
  720 continue
  710 continue
  730 continue
      do 81 ii=1,nwire
   81 call coildr(iwire(ii),jwire(ii),ilwire(ii))
      do 82 n=1,nlim
      x1 = xlima(n) + .5_R8*deex
      z1 = zlima(n) + .5_R8*deez
      x2 = xlima(n) - .5_R8*deex
      z2 = zlima(n) - .5_R8*deez
      x3 = xlima(n) + .5_R8*deex
      z3 = zlima(n) - .5_R8*deez
      x4 = xlima(n) - .5_R8*deex
      z4 = zlima(n) + .5_R8*deez
      call setcrt(x1,z1)
      call vector(x2,z2)
      call setcrt(x3,z3)
      call vector(x4,z4)
   82 continue
!
      do 83 n=1,nobs
      x1 = xobs(n)
      z1 = zobs(n) + .5_R8*deez
      x2 = xobs(n) - .5_R8*deex
      z2 = zobs(n)
      x3 = xobs(n)
      z3 = zobs(n) - .5_R8*deez
      x4 = xobs(n) + .5_R8*deex
      z4 = zobs(n)
      call setcrt(x1,z1)
      call vector(x2,z2)
      call vector(x3,z3)
      call vector(x4,z4)
      call vector(x1,z1)
   83 continue
      if(idiv.eq.0) go to 84
      call arrow(x2sep+2._R8*deex,z2sep+2._R8*deez,-2._R8*deex,-2._R8*   &  
     & deez)
      call arrow(x2sep+2._R8*deex,z1sep-2._R8*deez,-2._R8*deex,+2._R8*   &  
     & deez)
      call arrow(x1sep-2._R8*deex,z1sep-2._R8*deez,+2._R8*deex,+2._R8*   &  
     & deez)
      call arrow(x1sep-2._R8*deex,z2sep+2._R8*deez,+2._R8*deex,-2._R8*   &  
     & deez)
   84 continue
!
      if(iplate.eq.0 .or. nplate.eq.0) go to 8100
      do 8110 n=1,nplate
      call setcrt(xsega(n,1),zsega(n,1))
      do 8111 l=2,nseg(n)+1
 8111 call vector(xsega(n,l),zsega(n,l))
      if(isym.eq.0) go to 8110
      z1neg = -zsega(n,1)
      call setcrt(xsega(n,1),z1neg)
      do 8112 l=2,nseg(n)+1
      zlneg = -zsega(n,l)
 8112 call vector(xsega(n,l),zlneg)
 8110 continue
 8100 continue
      call frscj(3)
      if(idata.ne.7) goto 7989
! plot coils, wires, and preprogrammed shapes for Hofmann control
!
      xmin = ccon
      zmin = -alz
      xmax = alx
      zmax = alz
!
!.....external coils
      if(itemp.le.0) go to 7611
      do 7601 n=1,ncoil-nwire
      xmax = max(xmax,xcoil(n)+0.5_R8*dxcoil(n))
      xmin = min(xmin,xcoil(n)-0.5_R8*dxcoil(n))
      zmax = max(zmax,zcoil(n)+0.5_R8*dzcoil(n))
      zmin = min(zmin,-zmax)
 7601 continue
 7611 continue
      xave = .5_R8*(xmax+xmin)
      xdiff = xmax-xmin
      zdiff = zmax-zmin
      if((zdiff).gt.xdiff) xmax = xave+.5_R8*(zdiff)
      if((zdiff).gt.xdiff) xmin = xave-.5_R8*(zdiff)
      if((xdiff).gt.(zdiff)) zmax = .5_R8*(xdiff)
      if((xdiff).gt.(zdiff)) zmin = -.5_R8*(xdiff)
      call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.285_R8,1._R8)
      do 7301 ll=1,isym+1
      zfac = 1._R8
      if(ll.eq.2) zfac = -1._R8
!
      if(itemp.eq.0) go to 7901
!
!.....draw external coils
      do 7801 n=1,ncoil-nwire
      if(xcoil(n).ge.100._R8) go to 7801
      x1 = xcoil(n)-.5_R8*dxcoil(n)
      z1 = zcoil(n)-.5_R8*dzcoil(n)
      x2 = x1+dxcoil(n)
      z2 = z1
      x3 = x2
      z3 = z2+dzcoil(n)
      x4 = x1
      z4 = z3
      call setcrt(x1,z1*zfac)
      call vector(x2,z2*zfac)
      call vector(x3,z3*zfac)
      call vector(x4,z4*zfac)
      call vector(x1,z1*zfac)
 7801 continue
 7901 continue
 7301 continue
      do 811 ii=1,nwire
  811 call coildr(iwire(ii),jwire(ii),ilwire(ii))
! no. of preprogrammed bounndary points
      nhof1=acoef(380)
! no. of preprogrammed plasma shapes
      nhof2=acoef(340)
      do 8794 i=1,nhof2
      do 8793 j=1,nhof1
      xhof=acoef(2000+(i-1)*100+j)
      zhof=acoef(2030+(i-1)*100+j)
      if(j.eq.1) call setcrt(xhof,zhof)
      call vector(xhof,zhof)
8793  continue
      call vector(acoef(2001+(i-1)*100),acoef(2031+(i-1)*100))
      do 8792 j=1,nhof1
      xhof=acoef(2000+(i-1)*100+j)
      zhof=acoef(2030+(i-1)*100+j)
      x1=xhof+.5_R8*deex
      z1=zhof+.5_R8*deez
      x2=xhof-.5_R8*deex
      z2=zhof-.5_R8*deez
      x3=x1
      z3=z2
      x4=x2
      z4=z1
      if(j.eq.i) goto (51,52,53,54),i
      goto 8791
51    call setcrt(xhof+deex,z1)
      call vector(xhof+deex,z3)
      goto 8791
52    call setcrt(.5_R8*(xhof+x4)+deex,z1)
      call vector(.5_R8*(xhof+x1)+deex,z1)
      call vector(.5_R8*(xhof+x1)+deex,zhof)
      call vector(.5_R8*(xhof+x4)+deex,zhof)
      call vector(.5_R8*(xhof+x4)+deex,z3)
      call vector(.5_R8*(xhof+x1)+deex,z3)
      goto 8791
53    call setcrt(.5_R8*(xhof+x4)+deex,z1)
      call vector(.5_R8*(xhof+x1)+deex,z1)
      call vector(xhof+deex,zhof)
      call vector(.5_R8*(xhof+x1)+deex,z3)
      call vector(.5_R8*(xhof+x4)+deex,z3)
      goto 8791
54    call setcrt(xhof+deex,z3)
      call vector(xhof+deex,z1)
      call vector(x4+deex,zhof)
      call vector(x1+deex,zhof)
      goto 8791
8791  continue
      call setcrt(x1,z1)
      call vector(x2,z2)
      call setcrt(x3,z3)
      call vector(x4,z4)
 8792 continue
8794  continue
 7989 continue
      call frscj(3)
!
!.....plot zones where force balance is calculated
      call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.285_R8,1._R8)
      do 731 i=3,nxp
      do 731 j=3,nzp
      if(iexvc(i,j).gt.0) go to 731
      call setcrt(xary(i  ), zary(j  ))
      call vector(xary(i-1), zary(j  ))
      call vector(xary(i-1), zary(j-1))
      call vector(xary(i  ), zary(j-1))
      call vector(xary(i  ), zary(j  ))
      if(isym.eq.0) go to 731
      call setcrt(xary(i  ),-zary(j  ))
      call vector(xary(i-1),-zary(j  ))
      call vector(xary(i-1),-zary(j-1))
      call vector(xary(i  ),-zary(j-1))
      call vector(xary(i  ),-zary(j  ))
  731 continue
      call setld(1._R8,30._R8,1,0,1,1)
      write(s100,1731)
      call gtext(s100,80,0)
 1731 format(" interior zones where plasma can exist....")
      write(nsc1,1732)
 1732 format(" plasma zones")
      call frscj(6)
!
!.....coil resistivity for poloidal current calculation
      if(nwire.le.0) go to 1231

      call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.285_R8,1._R8)
      do 210 i=1,nxp
      do 210 j=1,nzp
  210 itap1(i,j) = 0
      do 211 n=1,nwire
      iw = iwire(n)
      jw = jwire(n)
  211 itap1(iw,jw) = 1
      do 212 n=1,nwire
      iw = iwire(n)
      jw = jwire(n)
      do 213 m=1,4
      call setcrt(xary(iw),zary(jw))
      go to(214,215,216,217),m
  214 if(itap1(iw+1,jw).eq.1) call vector(xary(iw+1),zary(jw))
!cj      write(76,7991) xary(iw+1),zary(jw)
!cj 7991 format('7991 ',1p2(2x,e15.7))     ! debug 
      go to 213
  215 if(itap1(iw-1,jw).eq.1) call vector(xary(iw-1),zary(jw))
!cj      write(76,7992) xary(iw+1),zary(jw)
!cj 7992 format('7992 ',1p2(2x,e15.7))     ! debug 
      go to 213
  216 if(itap1(iw,jw+1).eq.1) call vector(xary(iw),zary(jw+1))
!cj      write(76,7993) xary(iw+1),zary(jw)
!cj 7993 format('7993 ',1p2(2x,e15.7))     ! debug 
      go to 213
  217 if(itap1(iw,jw-1).eq.1) call vector(xary(iw),zary(jw-1))
!cj      write(76,7994) xary(iw+1),zary(jw)
!cj 7994 format('7994 ',1p2(2x,e15.7))     ! debug 
  213 continue
      if(isym.eq.0) go to 212
      do 223 m=1,4
      call setcrt(xary(iw),-zary(jw))
      go to(224,225,226,227),m
  224 if(itap1(iw+1,jw).eq.1) call vector(xary(iw+1),-zary(jw))
!cj      write(76,7995) xary(iw+1),zary(jw)
!cj 7995 format('7995 ',1p2(2x,e15.7))     ! debug 
      go to 223
  225 if(itap1(iw-1,jw).eq.1) call vector(xary(iw-1),-zary(jw))
!cj      write(76,7996) xary(iw+1),zary(jw)
!cj 7996 format('7996 ',1p2(2x,e15.7))     ! debug 
      go to 223
  226 if(itap1(iw,jw+1).eq.1) call vector(xary(iw),-zary(jw+1))
!cj      write(76,7997) xary(iw+1),zary(jw)
!cj 7997 format('7997 ',1p2(2x,e15.7))     ! debug 
      go to 223
  227 if(itap1(iw,jw-1).eq.1) call vector(xary(iw),-zary(jw-1))
!cj      write(76,7976) xary(iw+1),zary(jw)
!cj 7976 format('7976 ',1p2(2x,e15.7))     ! debug 
  223 continue
  212 continue
      call setld(1._R8,30._R8,1,0,1,1)
      write(s100,1216)
      call gtext(s100,80,0)
 1216 format(" current paths used in poloidal current calculation")
      write(nsc1,1217)
 1217 format(" poloidal current paths")
      call frscj(6)

!
!
!..rxw/02/03/88
 1231 continue
!
!.....switch and time step info:
      if(noplot(2).gt.0) goto 1232
!
      call map(.285_R8,1._R8,.285_R8,1._R8,0._R8,1._R8,0._R8,1._R8)
      call setld(1._R8,40._R8,1,0,1,0)
      write(s100,1001)(name(i),i=1,8)
      call gtext(s100,80,0)
 1001 format(1x,9a8,a7)
      write(s100,1010) udsd,udst,udsv,udsi,udsr,udse
      call gtextm(s100,80,0,1,6)
 1010 format(                                                            &  
     &" udsd =",1pe12.4,/,                                               &  
     &" udst =",1pe12.4,/,                                               &  
     &" udsv =",1pe12.4,/,                                               &  
     &" udsi =",1pe12.4,/,                                               &  
     &" udsr =",1pe12.4,/,                                               &  
     &" udse =",1pe12.4 )
 1421 format(/,"    irst1  ",i2,"    irst2  ",i2,"    ipest  "           &  
     & ,i2, "    isym   ",i2,"    idata  ",i2,"    lrswtch",i2)
 1422 format(/,                                                          &  
     &"    idens  ",i2,"    ipres  ",i2,"    ifunc  "                    &  
     & ,i2, "    iflux  ",i2,"    icube  ",i2,"    isurf  ",i2)
 1423 format(/,                                                          &  
     &"    itrmod ",i2,"    idiv   ",i2,"    iwall  "                    &  
     & ,i2, "    ialpha ",i2,"    iplate ",i2,"    ibalsw ",i2)
 1424 format(/,                                                          &  
     &"    icirc  ",i2,"    isvd   ",i2,"    isaw  "                     &  
     & ,i2, "    irfp   ",i2,"    itevv  ",i2,"    irippl ",i2)
 1425 format(/,                                                          &  
     &"    iimp   ",i2,"    ilte   ",i2,"    ibootst"                    &  
     & ,i2, "    iffac  ",i2)
!
      write(s100,1421) irst1,irst2,ipest,isym,idata,lrswtch
      call gtextm(s100,80,0,1,2)
      write(s100,1422) idens,ipres,ifunc,iflux,icube,isurf
      call gtextm(s100,80,0,2,2)
      write(s100,1423) itrmod,idiv,iwall,ialpha,iplate,ibalsw
      call gtextm(s100,80,0,2,2)
      write(s100,1424) icirc,isvd,isaw,irfp,itevv,irippl
      call gtextm(s100,80,0,2,2)
      write(s100,1425) iimp,ilte,ibootst,iffac
      call gtextm(s100,80,0,2,2)
!
!
      write(s100,1400)
      call gtext(s100,80,0)
      write(s100,1414) 1000._R8*dtmins
      call gtext(s100,80,0)
      write(s100,1415) 1000._R8*dtmaxs
      call gtext(s100,80,0)
      write(s100,1416) dtfac
      call gtext(s100,80,0)
      write(s100,1417) ndiv
      call gtext(s100,80,0)
      write(s100,1418) ffac
      call gtext(s100,80,0)
      write(s100,1419) tevv
      call gtext(s100,80,0)
      write(s100,1400)
      call gtext(s100,80,0)
      write(s100,1401) 1000._R8*dt1s
      call gtext(s100,80,0)
      write(s100,1402) 1000._R8*dt2s
      call gtext(s100,80,0)
      write(s100,1412) 1000._R8*dt3s
      call gtext(s100,80,0)
      write(s100,1403) 1000._R8*t1
      call gtext(s100,80,0)
      write(s100,1404) 1000._R8*t2
      call gtext(s100,80,0)
      write(s100,1405) 1000._R8*t3
      call gtext(s100,80,0)
      write(s100,1413) 1000._R8*t3d
      call gtext(s100,80,0)
      write(s100,1406) 1000._R8*t4
      call gtext(s100,80,0)
      call frscj(11)
!
!..rxw/02/03/88
 1232 continue
!
      if(ineg.eq.3) go to 912
!     if(isym.eq.1) go to 750
      if(lrswtch.ge.1) go to 741
      xmin = abs(xdist) + .0001_R8
      if(xmin.lt. .5_R8*(ccon+xplas)) xmin = .5_R8*(ccon+xplas)
      xmax = 2._R8*xplas-xmin
      if(xmax.gt. .5_R8*(alx+xplas))  xmax = .5_R8*(alx+xplas)
      dx = (xmax-xmin)/12._R8
!
!..rxw/02/02/88
!.....filiment growth model:
      if(noplot(3).gt.0) goto 741
!
!.......bypass growth calculation for nwire .gt. 50
      if(ncoil .ge. pnx .or. acoef(28) .gt. 0.0_R8) go to 741
      call map(.285_R8,1._R8,.285_R8,1._R8,0._R8,1._R8,0._R8,1._R8)
      call setld(1._R8,40._R8,1,0,1,0)
      write(nout,1925)
      write(s100,1925)
      call gtextm(s100,80,0,1,3)
 1925 format(//,"    position    growth rate  index")
!....................................................................
!......part 5:  calculate growth rates for filiment model located at origin
!.......................................................................
!
      do 740 n=1,12
      xp = xmin + (n-1)*dx
      call growth(xp,ratio,alams,aindx)
      write(s100,1926) xp,alams,aindx
      call gtext(s100,80,0)
      write(nout,1926) xp,alams,aindx
  740 continue
      call frscj(4)
  741 continue
      if(ineg.ne.0) return
!
!.....calculate field gradients and multipole moments of current groups
      call multipol
!
!
 1926 format(1x,1p5e12.4)
  750 continue
!
!..rxw/02/03/88
!.....initial coil and wire info:
      if(noplot(4).gt.0) return
!
      if(ncoil.le.0) return
      call map(.285_R8,1._R8,.285_R8,1._R8,0._R8,1._R8,0._R8,1._R8)
      call setld(1._R8,40._R8,1,0,1,0)
      ilct = 0._R8
      if(nwire.le.0) go to 910
      write(s100,1908)
      call gtext(s100,80,0)
 1908 format("   n   xwire      zwire     iw  jw igrw ntrn",             &  
     &         "    rswires    aindw     lr-time")
      do 909 n=1,nwire
      ilct = ilct + 1
      if(ilct.eq.29) call frscj(5)
      if(ilct.eq.29) call map(.285_R8,1._R8,.285_R8,1._R8,0._R8,1._R8,   &  
     & 0._R8,1._R8)
      if(ilct.eq.29) call setld(1._R8,40._R8,1,0,1,0)
      if(ilct.eq.29) ilct = 0
      if(rswires(n).ne.0) go to 908
      ineg=5
      go to 909
  908 continue
      timelr = aindw(n)/rswires(n)
      write(s100,1909) n,xwire(n),zwire(n),iwire(n),jwire(n),igroupw(n),  &  
!    &                                                                   &  
     &  aturnsw(n),rswires(n),aindw(n),timelr
      call gtext(s100,80,0)
  909 continue
  910 if(ncoil-nwire.le.0) go to 913
      write(s100,1910)
      call gtext(s100,80,0)
 1910 format("   n   xcoil      zcoil           igrc ntrn")
      do 911 n=1,ncoil-nwire
      ilct = ilct+1
      if(ilct.eq.29) call frscj(5)
      if(ilct.eq.29) call map(.285_R8,1._R8,.285_R8,1._R8,0._R8,1._R8,   &  
     & 0._R8,1._R8)
      if(ilct.eq.29) call setld(1._R8,40._R8,1,0,1,0)
      if(ilct.eq.29) ilct=1
      write(s100,1911) n,xcoil(n),zcoil(n),igroupc(n),aturnsc(n)
      call gtext(s100,80,0)
  911 continue
 1911 format(i4,1p2e11.3,10x,i2,0pf6.1,1p2e11.3)
  912 continue
 1909 format(i4,1p2e11.3,3i4,0pf6.1,1p4e11.3)
  913 if(ilct.ne.0) call frscj(5)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
