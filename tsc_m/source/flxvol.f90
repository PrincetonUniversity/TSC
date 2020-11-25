!#include "f77_dcomplx.h"
      subroutine flxvol(itype)
!......8.20 flxvol
!
!.....this subroutine does line integrals around constant flux surfaces
!.....to compute geometrical quantities
!
!
!.....for itype=1  equal increments in poloidal flux
!.......  itype=2  equal increments in toroidal flux
!........ itype=3  equal poloidal flux, initial equilibrium with ifunc=5 only
!
      USE CLINAM
      USE BALCLI
      USE RUNAWAY
      USE SAPROP
      USE SCR1 
      USE SCR3
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,itroub,isw,iinit,ifast,nplot,npartsm,j,i,ip,l
      INTEGER iqtrub,nparts2,nparts,k,ihalf,m,kk,kmaxs,iis,ii,ifs
      INTEGER ii66,ios66,ifail,lval,ll,llsave
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 czsv,cxsv2,xshape,fshape
      REAL*8 psi100,angl1,angl2,angl3
      REAL*8 angl4,pcon,vsave,gradsq,dpsidx,dpsidz,gsval,psval
      REAL*8 psixz,psixx,psizz,gsval1,psixx1,psizz1,gsval2,psixx2
      REAL*8 psizz2,akap,term,th1,qprof2o,qprof2n,qdiff,pinc
      REAL*8 qproflp1,pinco,qave,pcono,roffp,cz1,cr1,cz2,cr2,cz3
      REAL*8 cr3,cz4,cr4,av1,av2,av3,av4,av5,av11,av6,av7,av8,av9
      REAL*8 av10,av12,av13,av14,alam,alamf,area,dtheta,theta,bmax,cosine
      REAL*8 sine,y,x3,z3,psi3,x2,z2,psi2,dum1,dum2,dum3,dum4,dum5
      REAL*8 dum6,dum7,y2,y3,xpass,zpass,grs,gs1,pt,ysv1,grssv1
      REAL*8 pr1,pr2,pr1sv1,prave
      REAL*8 gs1sv1,pp,dy,ysave,xcont,zcont,btry,rtheta,rthetl
      REAL*8 xright,xleft,r3,z1,r1,z4,r4,am1,am2,dsdr,dxsq,grso,gs2
      REAL*8 gsave,grsa,rave,zave,r2,fac,bsq,bmag,bsqi,arg,brac
      REAL*8 psix,psiz,ajp3x,ajp3z,theval,sarym2,d1,term1,d2,term2
      REAL*8 d3,term3,t1,t2,t3,t4,rmaxfs,rminfs,zmaxfs,rzmaxfs
      REAL*8 zminfs,rzminfs,perim,pval,ppval,eval,epval,rval,rpval
      REAL*8 sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,gval
      REAL*8 gpval,gppval,dtfluxv,apld,sum10,tfcon,pvpgtempo,paxis
      REAL*8 ppaxis,prel,pboundary,ps,vptemp,vpgtemp,preturn
      REAL*8 pvpgtemp,dpvg
      REAL*8 err, AREAL
!============
      dimension czsv(3),cxsv2(3),                                        &  
     &          xshape(4),fshape(4)
!              tftemp(pnx+ppsi),
!    2          sary(pnparts),
!    3 yold(ppsi),x3old(ppsi),z3old(ppsi),psi3old(ppsi),term1aj(ppsi)
!     dimension grssv(pnparts),gs1sv(pnparts)
      data itroub/1/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tftemp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: yold
      REAL*8, ALLOCATABLE, DIMENSION(:) :: x3old
      REAL*8, ALLOCATABLE, DIMENSION(:) :: z3old
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psi3old
      REAL*8, ALLOCATABLE, DIMENSION(:) :: term1aj
      REAL*8, ALLOCATABLE, DIMENSION(:) :: avepres
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grssv
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gs1sv
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pr1sv
!============      
      IF(.not.ALLOCATED(tftemp)) ALLOCATE( tftemp(pnx+ppsi), STAT=istat)
      IF(.not.ALLOCATED(sary)) ALLOCATE( sary(pnparts), STAT=istat)
      IF(.not.ALLOCATED(yold)) ALLOCATE( yold(ppsi), STAT=istat)
      IF(.not.ALLOCATED(x3old)) ALLOCATE( x3old(ppsi), STAT=istat)
      IF(.not.ALLOCATED(z3old)) ALLOCATE( z3old(ppsi), STAT=istat)
      IF(.not.ALLOCATED(psi3old)) ALLOCATE( psi3old(ppsi), STAT=istat)
      IF(.not.ALLOCATED(term1aj)) ALLOCATE( term1aj(ppsi), STAT=istat)
      IF(.not.ALLOCATED(avepres)) ALLOCATE( avepres(ppsi), STAT=istat)
      IF(.not.ALLOCATED(grssv)) ALLOCATE( grssv(pnparts), STAT=istat)
      IF(.not.ALLOCATED(gs1sv)) ALLOCATE( gs1sv(pnparts), STAT=istat)
      IF(.not.ALLOCATED(pr1sv)) ALLOCATE( pr1sv(pnparts), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : flxvol  ' 
!============      
!
!     open(66,file='diagout',status='unknown',iostat=ios)
!     do i=2,nxp
!     write(66,6666) xary(i),psi(i,2)
!     enddo
!6666 format(1x,1p2e12.4)
!     stop
      if(irst1.ne.0) ifrst(8) = 0
      if(itype.eq.3) npsit = npsi/tfmult
      if(itype.eq.1) npsit = npsi
!     write(nout,9901) kcycle
!     write(nout,9902)
!     write(nout,9903) (xsv2(j),j=1,npsit)
!     write(nout,9905)
!     write(nout,9903) (qprof2(j),j=1,npsit)
!     write(nout,7701) kcycle, npsi,npsit,itype,ppsi,tfmult
!7701 format(" kcycle, npsi, npsit, itype, ppsi, tfmult =", 5i5, 1pe12.4)

!9901 format(" cycle=",i7)
!9903 format(1p10e12.4)
!9902 format(" xsv2 array")
!9905 format(" qprof2 array")
!
!.....locate index of 95% flux surface,   i95
      if(itype.ne.2) go to 97
!
      psi100 = psilim
      if(numsep.gt.0.and.psisep.ge.psilim) psi100=psisep
!
      psi95 = .05_R8*psimin + .95_R8*psi100
      if(psi95.gt.psilim) psi95 = psilim
!
      psi90 = .10_R8*psimin + .90_R8*psi100
      if(psi90.gt.psilim) psi90 = psilim
   97 continue
!
      isw=1
      iinit = 0
      ifast = 0
      nplot = 1
      npsim = npsi-1
      npartsm = ((nx-1)+(isym+1)*(nz-1))
      err = 1.E-12_R8
      angl1 = atan2(deex,deez)
      angl2 = pi - angl1
      angl3 = pi + angl1
      angl4 = 2._R8*pi - angl1
      pcon = psimin
      vsave = 0._R8
      ajpress = 0._R8
      xsv2(1) = psimin
      call grap(0,zmag,xmag,gradsq,dpsidx,dpsidz,gsval,psval,            &  
     &                psixz,psixx,psizz,1)
      if(isym.eq.0) go to 777
      j=2
      i = (xmag-ccon)/deex+2
      ip = i+1
      call grap(0,zmag,xary(i),gradsq,dpsidx,dpsidz,gsval1,psval,        &  
     &                psixz,psixx1,psizz1,1)
      call grap(0,zmag,xary(ip),gradsq,dpsidx,dpsidz,gsval2,psval,       &  
     &                psixz,psixx2,psizz2,1)
      psixx = psixx1*(xary(ip) - xmag)/deex                              &  
     &       +psixx2*(xmag - xary(i) )/deex
      psizz = psizz1*(xary(ip) - xmag)/deex                              &  
     &       +psizz2*(xmag - xary(i) )/deex
  777 continue
      akap = 1._R8
      if(psixx.gt.0 .and. psizz.gt.0)                                    &  
     &akap = sqrt(psixx/psizz)
      term = (akap + 1._R8/akap) / (psixx + psizz)
      qprof2p = gsval * term / xmag
!ccccccccccccc  2/22/93 SCJ
!     qprof2(1) = qprof2(2)
!..........5/12/02 SCJ
      qprof2(1) = 2._R8*qprof2(2) - qprof2(3)
!
!ccccccccccccc
!
!.....start of loop over psi values.
!
!     note:   npsi is max no of toroidal flux intervals available
!             npsit is the current no in the plasma
!             must have npsit .lt. npsi always
!
      i95 = 0
!     write(nterm,3333)
 3333 format(" i95 set to zero")
      i90 = 0
      iqtrubmax = 0
!
!.....start loop over flux surfaces
      do 600 l=1,npsi-1
      iqtrub = 0
      th1 = 0.5_R8
      go to 96
   95 continue
      iqtrub = iqtrub + 1
!
!....note fix 7/15/98 (SCJ)
!    th1 starts out at 0.5 and decreases as iqtrub increases
      th1 = th1*.96_R8
      if(iqtrub.ge.190)                                                  &  
     & write(nterm,1097) kcycle,iqtrub,l,qprof2o,qprof2n,qdiff,th1
      if(iqtrub.ge.190)                                                  &  
     & write(nout,1097) kcycle,iqtrub,l,qprof2o,qprof2n,qdiff,th1
 1097 format("  --> qtrub diagnositc",3i6,1p4e12.4)
      if(iqtrub.le.199) go to 96
      write(nout,1096) l
      write(nterm,1096) l
 1096 format(" ERROR --> more than 99 iterations in flxvol on q, l=",i4)     
      ineg=47
      return
   96 continue
!
      nparts2 =2 *sqrt(l*(npartsm-8)**2/(4._R8*(npsit-1))) + 8
      nparts = 2*nparts2
      if(nparts.gt.pnparts-6) nparts = pnparts-6
      go to(71,72,74),itype
   71 pinc = (psilim-psimin)/(npsi)
      go to 73
   72 pcon = xsv2(l)
      if(iqtrub.ge.1) then
        qproflp1 = th1*qprof2(l+1)+(1._R8-th1)*qprof2o
        pinco = pinc
        if (l.eq.1) qprof2(1) = qproflp1
        qave = 0.5_R8*(qprof2(l)+qproflp1)
        pinc = dpsi/qave/tpi
!
      else
        qave = .5_R8*(qprof2(l)+qprof2(l+1))
        pinc = dpsi/qave/tpi
      endif
!
      go to 73
   74 pinc = tfmult*(psilim-psimin)/npsi
   73 continue
      pcono = pcon
      pcon = pcon+pinc
!
      if(itype.ne.2) go to 75
      if(xsv2(l).lt.psi95 .and. pcon .ge. psi95) i95 = l+1
      if(xsv2(l).lt.psi90 .and. pcon .ge. psi90) i90 = l+1
      if(pcono.le.psilim .and. pcon.gt.psilim)                           &  
     &  fraclst = (psilim-pcono)/(pcon-pcono)
   75 continue
      roffp = abs(pcon*1.E-12_R8)
      if(pcon.gt.psilim+roffp) go to 598
      cz1  = 0._R8
      cr1  = 0._R8
      cz2  = 0._R8
      cr2  = 0._R8
      cz3  = 0._R8
      cr3  = 0._R8
      cz4  = 0._R8
      cr4  = 0._R8
      av1  = 0._R8
      av2  = 0._R8
      av3  = 0._R8
      av4  = 0._R8
      av5  = 0._R8
      av11 = 0._R8
      av6  = 0._R8
      av7 = 0._R8
      av8 = 0._R8
      av9 = 0._R8
      av10 = 0._R8
      av12 = 0._R8
      av13 = 0._R8
      av14 = 0._R8
!
!
!.....define hyper term here
      if(hypermult.le.0._R8) go to 6464
      alam = alamf(pcon)
 6464 continue
      area = 0._R8
      dtheta = 2._R8*pi/nparts
      theta = -(pi/2+2._R8*dtheta)
      bmax = 0._R8
!
!.....start of loop over angles.
!.....angles are measured clockwise starting from the vertical
      kmaxo = kmax
      kmax = nparts
      if(isym.eq.1) kmax = nparts/2+2
      do 500 k=1,kmax
      theta = theta + dtheta
      if(k.eq.1 .and. isym.eq.1) go to 500
      if(k.eq.2 .and. isym.eq.1) iinit=0
      ihalf = 1
      cosine = cos(theta)
      sine = sin(theta)
      ifast = 0
      if(iinit.ne.0) go to 89
   13 iinit = 1
      ifast = 1
      y = 0
      x3 = xmag
      z3 = zmag
      psi3 = psimin
      if(l.eq.1) go to 10
      y = yold(l-1)
      x3 = x3old(l-1)
      z3 = z3old(l-1)
      psi3 = psi3old(l-1)
   10 x2 = x3
      z2 = z3
      psi2 = psi3
      y = y+nx*deex/(ihalf*npsi*4)
      x3 = xmag+y*sine
      z3 = zmag+y*cosine
      if(x3.gt.alx        ) go to 41
      if(x3.lt.ccon+deex  ) go to 41
      if(z3.gt.alz        ) go to 41
      if(z3.lt.-alz        ) go to 41
      call grap(1,z3,x3,dum1,dum2,dum3,dum4,psi3,                        &  
     &          dum5,dum6,dum7,isw)
      isw=0
!     if(itype.eq.2) write(nout,2995) kcycle,l,iqtrub,k,ihalf,iinit,ifast
!2995 format(" kcycle,l,iqtrub,k,ihalf,iinit,ifast",7i5)
!     if(itype.eq.2) write(nout,2005) x3,z3,pcon,psi2,psi3
!2005 format(' gdiag,x3,z3,pcon,psi2,psi3',1p5e12.4)
      if(psi3.ge.pcon-roffp) go to 40
      go to 10
   41 continue
!     if(ihalf.gt.5)
!    1write(nout,3011) l,k,ihalf,x2,x3,z2,z3,psi2,psi3,pcon
!3011 format(' 41diag',3i3,1p7e12.4)
      ihalf = ihalf + 1
      if(ihalf.gt.8) go to 87
      go to 13
   40 continue
!.....y is distance along ray from origin.
      y2 = cosine*(z2 - zmag) + sine*(x2 - xmag)
      y3 = cosine*(z3 - zmag) + sine*(x3 - xmag)
!.....solve for intersection point using newtons method
      y =.5_R8*(y2+y3)
   89 continue
      do 85 m=1,8
      xpass = xmag+y*sine
      zpass = zmag+y*cosine
      if(xpass.le.ccon.or.xpass.ge.alx) go to 85
      if(zpass.le.-alz.or.zpass.ge.alz) go to 85
      call grap(2,zpass,xpass,grs ,dpsidx,dpsidz,gs1 ,pt,                &  
     &          dum3,dum4,dum5,isw)
      call evalpr(zpass,xpass,pr1)
      isw = 0
      if( m.gt.1 ) go to 84
      ysv1   = y
      grssv1 = grs
      gs1sv1 = gs1
      pr1sv1 = pr1
   84 continue
      pp = sine*dpsidx + cosine*dpsidz
!.....write(nout,2008) xpass,zpass,pt,pp
!2008 format(' hdiag,xpass,zpass,pt,pp',1p4e12.4)
      dy = (pcon-pt)/pp
!.....added 12/14/85
      if(abs(dy) .gt. deex) go to 185
      y = y+dy
      if(y.lt.0._R8) go to 185
      if( m.ge.2 .and. abs(dy)      .lt. 1.0E-8_R8*deex) go to 86
      if( m.ge.2 .and. abs(pcon-pt) .lt. 1.0E-6_R8*pinc) go to 86
   85 continue
  185 continue
!
      if(ifast.eq.0) go to 13
      go to 41
   87 continue
      if( pcon.ge.psi2 .and. pcon.le.psi3 ) go to 88
      write(nout,703) l,k,pcon,x2,z3,psi2,x3,z3,psi3
  703 format(' flxvoltrouble',2i4,1p7e12.4)
      itroub = itroub + 1
      if(itroub.le.0) go to 704
      ineg=14
      npsit = l
      return
  704 continue
      y = ysave
      go to 86
   88 continue
      y      = ysv1
      grs    = grssv1
      gs1    = gs1sv1
      pr1    = pr1sv1
!.....zcont and xcont are the coordinates of
!.....the point where ray number k intersects
!.....the psi surface number l.
   86 xcont = sine*y + xmag
      zcont = cosine*y + zmag
      ysave = y
      xplot(2,k) = xplot(1,k)
      zplot(2,k) = zplot(1,k)
      xplot(1,k) = xcont
      zplot(1,k) = zcont
      grssv(k)=grs
      gs1sv(k)=gs1
      pr1sv(k)=pr1
!
!.....max b on flux surface
      btry = sqrt((grs+gs1**2)/xcont**2)
      bmax = max(bmax,btry)
!.....evaluation of minor radius
      rtheta = abs(theta-(pi/2))
      rthetl = abs(theta+(pi/2))
      if(rtheta.lt.1.E-6_R8) xright = abs(sine*y)
      if(rthetl.lt.1.E-6_R8) xleft = abs(sine*y)
      if(isym.eq.1 .and. k.ne.2) go to 500
      if(isym.eq.0 .and. k.ne.1) go to 500
      yold(l) = y
      x3old(l) = xpass
      z3old(l) = zpass
      psi3old(l) = pt
  500 continue
      qamin(l+1)=(xright+xleft)/2
      if(isym.eq.0) go to 801
      zplot(2,kmax+1) = zplot(1,kmax+1)
      xplot(2,kmax+1) = xplot(1,kmax+1)
      zplot(1,kmax+1) = -zplot(1,kmax-1)
      xplot(1,kmax+1) = xplot(1,kmax-1)
      grssv(kmax+1) = grssv(kmax-1)
      gs1sv(kmax+1) = gs1sv(kmax-1)
      pr1sv(kmax+1) = pr1sv(kmax-1)
      zplot(2,1) = zplot(1,1)
      xplot(2,1) = xplot(1,1)
      zplot(1,1) = -zplot(1,3)
      xplot(1,1) = xplot(1,3)
      grssv(1) = grssv(3)
      gs1sv(1) = gs1sv(3)
      pr1sv(1) = pr1sv(3)
      go to 803
  801 continue
      do 802 kk=1,3
      zplot(2,kmax+kk) = zplot(1,kmax+kk)
      xplot(2,kmax+kk) = xplot(1,kmax+kk)
      zplot(1,kmax+kk) = zplot(1,kk)
      xplot(1,kmax+kk) = xplot(1,kk)
      grssv(kmax+kk) = grssv(kk)
      gs1sv(kmax+kk) = gs1sv(kk)
      pr1sv(kmax+kk) = pr1sv(kk)
  802 continue
  803 continue
!
      if(itype.ne.2 .or. l+1.ne.i95) go to 804
      do 805 k=1,kmax+3
      xplot(4,k) = xplot(1,k)
      zplot(4,k) = zplot(1,k)
  805 continue
      do 806 k=1,kmaxo+3
      xplot(3,k) = xplot(2,k)
      zplot(3,k) = zplot(2,k)
  806 continue
      kkm95 = kmax
      kkm95m = kmaxo
  804 continue
!
      if(itype.ne.2 .or. l+1.ne.i90) go to 904
      do 905 k=1,kmax+3
      xplot(6,k) = xplot(1,k)
      zplot(6,k) = zplot(1,k)
  905 continue
      do 906 k=1,kmaxo+3
      xplot(5,k) = xplot(2,k)
      zplot(5,k) = zplot(2,k)
  906 continue
      kkm90 = kmax
      kkm90m = kmaxo
  904 continue
!
!.....evaluate contour integrals
      cz2 = zplot(1,1)
      cr2 = xplot(1,1)
      cz3 = zplot(1,2)
      cr3 = xplot(1,2)
      cz4 = zplot(1,3)
      cr4 = xplot(1,3)
      grs = grssv(2)
      gs1 = gs1sv(2)
      pr1 = pr1sv(2)
      kmaxs = kmax+1
      if(isym.eq.0) kmaxs = kmax+3
      do 501 k=4,kmaxs
      xcont = xplot(1,k)
      zcont = zplot(1,k)
      cz1 = cz2
      cr1 = cr2
      cz2 = cz3
      cr2 = cr3
      cz3 = cz4
      cr3 = cr4
      cz4 = zcont
      cr4 = xcont
      sine = cz2-cz3
      cosine = -(cr2-cr3)
      r3 = -sine*(cz3-cz2) + cosine*(cr3-cr2)
      z1 =  cosine*(cz1-cz2) + sine*(cr1-cr2)
      r1 = -sine*(cz1-cz2) + cosine*(cr1-cr2)
      z4 =  cosine*(cz4-cz2) + sine*(cr4-cr2)
      r4 = -sine*(cz4-cz2) + cosine*(cr4-cr2)
!.....slope at points 2 and 3.
      am1 = .5_R8*z1/r1
      am2 = .5_R8*z4/(r4-r3)
!
!.....arc length multiplier
      dsdr = 1._R8+(2._R8*am1**2 + 2._R8*am2**2 - am1*am2)/30._R8
      dxsq = sine**2 + cosine**2
      grso = grs
      gs2 = gs1
      pr2 = pr1
      grs = grssv(k-1)
      gs1 = gs1sv(k-1)
      pr1 = pr1sv(k-1)
      gsave = .5_R8*(gs1+gs2)
      prave = .5_R8*(pr1+pr2)
      grsa = .5_R8*(grs+grso)
      rave = .5_R8*(cr2+cr3)
      zave = .5_R8*(cz2+cz3)
      if(rave.le.0._R8.or.grsa.le.0._R8) go to 501
      r2 = rave**2
      fac = sqrt(dxsq/grsa)*dsdr*tpi*rave
!
!.....averages needed for trapped particle fraction
      bsq = (grsa + gsave**2)/rave**2
      bmag = sqrt(bsq)
      bsqi = 1._R8/bsq
      arg = max(1._R8-bmag/bmax , 1.E-6_R8)
      brac = arg**0.5_R8- 0.33_R8*arg**1.5_R8
      av3 = av3 + fac*bsq
      av4 = av4 + fac*bsqi
      av7 = av7 + fac*bsqi*brac
!
      av1 = av1 + fac
      av2 = av2 + fac*gsave/r2
      av5 = av5 + fac/r2
      av11 = av11 + fac*gsave**2/r2
      av6 = av6 + sqrt(dxsq*grsa)*dsdr*tpi/rave
      av12 = av12 + sqrt(dxsq*grsa)*dsdr*tpi*rave
!     if(l+1.eq.50.and.kcycle.gt.5787) then
!     write(nterm,2222) kcycle,lp,itroub,rave,zave,pt,psilim
!2222 format("kcycle,lp,itroub,rave,zave,pt,psilim",3i5,1p4e12.4)
!     endif
      av8 = av8 + sqrt(dxsq)*dsdr
      av10 = av10 + sqrt(dxsq)*dsdr
      av13 = av13 + fac*prave
!
!
!.....compute average needed for hyperres heating term
      if(hypermult.le.0) go to 6465
      i = (xcont-ccon)/deex+2
      j = (zcont-zzero)/deez+2
      psix = .5_R8*(psi(i+1,j+1)+psi(i+1,j)                              &  
     &          -psi(i  ,j+1)-psi(i  ,j) )/deex
      psiz = .5_R8*(psi(i+1,j+1)+psi(i,j+1)                              &  
     &          -psi(i+1,j  )-psi(i,j  ) )/deez
      ajp3x = .5_R8*(ajp3(i+1,j+1)+ajp3(i+1,j)                           &  
     &          -ajp3(i  ,j+1)-ajp3(i  ,j) )/deex
      ajp3z = .5_R8*(ajp3(i+1,j+1)+ajp3(i,j+1)                           &  
     &          -ajp3(i+1,j  )-ajp3(i,j  ) )/deez
      term = fac*(psix*ajp3x + psiz*ajp3z)
      av9 = av9 + term
 6465 continue
      if((iwall.eq.0.and.ibalsw.eq.0.and.irippl.eq.0).or.itype.ne.2)     &  
     &  go to 501
!     sary(k-2) = av5
!.....equal arc length jacobian
      sary(k-2) = av8
  501 continue
!
!.....skip over unless wall calculations are included
      if((iwall.eq.0.and.ibalsw.eq.0.and.irippl.eq.0).or.itype.ne.2)     &  
     &  go to 503
      sary(1) = 0._R8
      fac = (2._R8-AREAL(isym))*pi/sary(kmaxs-2)
      do 505 k=1,kmaxs-2
  505 sary(k) = sary(k)*fac
      if(isym.eq.1) sary(kmaxs-1) = 2._R8*pi-sary(kmaxs-3)
      if(isym.eq.0) sary(kmaxs-1) = 2._R8*pi+sary(2)
      i = 2
      iis = 3
      do 506 ii = iis,nthe+1
      theval = AREAL(ii-iis+1)*(2._R8-AREAL(isym))*pi/nthe
  540 continue
      if(i.ge.kmaxs-2) go to 550
      if(theval.lt.sary(i)) go to 550
      i = i+1
      if(i.ge.kmaxs-2) go to 550
      go to 540
  550 continue
      sarym2 = sary(i-2)
      if(i.eq.2.and.isym.eq.0) sarym2 = sary(kmaxs-3)-2*pi
      if(i.eq.2.and.isym.eq.1) sarym2 =-sary(2)
!
!.....theval lies between sary(i-1) and sary(i)
!     cubic angle interpolation
      d1 = (sary(i)-sary(i-1))**3
      term1 = (theval-sary(i-1))**2*(3._R8*sary(i)-sary(i-1)-2._R8*      &  
     & theval)/d1
      d2 = (sary(i)-sary(i-1))**2*(sary(i)-sarym2)
      term2 = (theval-sary(i))**2*(theval-sary(i-1))/d2
      d3 = (sary(i)-sary(i-1))**2*(sary(i+1)-sary(i-1))
      term3 = (theval-sary(i-1))**2*(theval-sary(i))/d3
      t1 = -term2
      t2 = 1._R8-term1 - term3
      t3  = term1 + term2
      t4 = term3
      xw(ii,l+1)=t1*xplot(1,i-1)+t2*xplot(1,i)+t3*xplot(1,i+1)           &  
     &          +t4*xplot(1,i+2)
      zw(ii,l+1)=t1*zplot(1,i-1)+t2*zplot(1,i)+t3*zplot(1,i+1)           &  
     &          +t4*zplot(1,i+2)
  506 continue
      if(isym.gt.0) go to 1506
      xw(2,l+1) = xplot(1,2)
      zw(2,l+1) = zplot(1,2)
      xw(1,l+1) = xw(nthe+1,l+1)
      zw(1,l+1) = zw(nthe+1,l+1)
      xw(nthe+2,l+1) = xw(2,l+1)
      zw(nthe+2,l+1) = zw(2,l+1)
      xw(nthe+3,l+1) = xw(3,l+1)
      zw(nthe+3,l+1) = zw(3,l+1)
!
      go to 503
 1506 continue
      xw(2,l+1) = xplot(1,2)
      zw(2,l+1) = zplot(1,2)
      xw(1,l+1) = xw(3,l+1)
      zw(1,l+1) =-zw(3,l+1)
      xw(nthe+2,l+1) = xplot(1,kmax)
      zw(nthe+2,l+1) = zplot(1,kmax)
      xw(nthe+3,l+1) = xw(nthe+1,l+1)
      zw(nthe+3,l+1) =-zw(nthe+1,l+1)
!
  503 continue
      if(isym.gt.0) go to 2506
      rmaxfs = xplot(1,2)
      rminfs = xplot(1,2)
      zmaxfs = zplot(1,2)
      rzmaxfs= xplot(1,2)
      zminfs = zplot(1,2)
      rzminfs= xplot(1,2)
      do ifs=3,kmax
        rmaxfs = max(rmaxfs,xplot(1,ifs))
        rminfs = min(rminfs,xplot(1,ifs))
        if(zplot(1,ifs).gt.zmaxfs) then
          zmaxfs = zplot(1,ifs)
          rzmaxfs= xplot(1,ifs)
        endif
!
        if(zplot(1,ifs).lt.zminfs) then
          zminfs = zplot(1,ifs)
          rzminfs= xplot(1,ifs)
        endif
!
      enddo
      go to 2503
 2506 continue
      rmaxfs = xplot(1,kmax)
      rminfs = xplot(1,2)
      zmaxfs = zplot(1,2)
      rzmaxfs= xplot(1,2)
      do ifs=3,kmax
        if(zplot(1,ifs).gt.zmaxfs) then
          zmaxfs = zplot(1,ifs)
          rzmaxfs= xplot(1,ifs)
        endif
      enddo
      zminfs = -zmaxfs
      rzminfs= rzmaxfs
 2503 continue
!
!.....end of loop over angles
      qprof2o = qprof2(l+1)
      qprof2(l+1) =av2*(1+isym)/tpi**2
!
!.....added 5/17/2013 (SCJ)
      if(qprof2(l+1) .gt. 100.) qprof2(l+1) = 100.
!
      qprof2n = qprof2(l+1)
      if(itype.eq.2) then
      if(kcycle.ge.0) qprof2n = th1*qprof2n+(1._R8-th1)*qprof2o
!
!.....added 7/19/98
      if(l.eq.1 .and. iqtrub.gt.100) then
      qprof2n = qprof2(3)
      endif
!
      qprof2(l+1)=qprof2n
      qdiff = abs(qprof2n-qprof2o)
      if(abs(qprof2n-qprof2o) .gt. .0005_R8*qprof2o)                     &  
     &go to 95
      endif
!
      vp2(l+1) = av1*(1+isym)/(tpi*qprof2(l+1))
      xmja2(l+1) = av5*(1+isym)
      x2ave(l+1) = av1/av5
      gxmja2(l+1) = av6*(1+isym)
      gja2(l+1) = av12*(1+isym)
      avhyp(l+1) = alam*av9*(1+isym)*acoef(56)
      xsv2(l+1) = pcon
      vsave = vsave + av1*(1+isym)*(xsv2(l+1)-xsv2(l))
      vary(l+1) = vsave
      term1aj(l+1) =  av1*(1+isym)*(xsv2(l+1)-xsv2(l))                   &  
     &                  *(av11/av3 - 1._R8)*udsi/tpi
      perim = av10*(1+isym)
      bpolar(l+1) = gxmja2(l+1)*1.256637E-2_R8/(perim*usdi*tpi)
!.....added 1/28/02
      rminora(l+1) = (rmaxfs-rminfs)/2._R8
      rmajora(l+1) = (rmaxfs+rminfs)/2._R8
      elonga(l+1) =  (zmaxfs-zminfs)/(rmaxfs-rminfs)
      arg =( rmajora(l+1)-0.5_R8*(rzmaxfs+rzminfs))/rminora(l+1)
      deltaa(l+1) = arg
      if(abs(arg).lt.1._R8) deltaa(l+1) = asin(arg)
!.....added 6/09/03
      ajave(l+1)  = av1
      ajbtsq(l+1) = av11
      ajbsq(l+1)  = av3
      avepres(l+1) = av13/av1
!
      npsit = l+1
      if(itype.eq.2)npsitmx = max( npsit,npsitmx)
!
!.....particle trapping fraction
!                               ref: Hirshman and Jardin
!                                    Phys Fluids 22 (p 731) 1979
!
      bsq = av3/av1
      bsqi = av4/av1
      bsqar(l+1) = bsq*1.E8_R8
      bsqiar(l+1) = bsqi*1.E-8_R8
      brac = av7/av1
      ftrap(l+1) = 1._R8+bsq*( -bsqi + 1.5_R8*brac)
      if(acoef(109).gt.0) ftrap(l+1) = 0._R8
      if(vary(l+1).le.0) ineg=14
!
      iqtrubmax = max(iqtrub, iqtrubmax)
      go to 600
  598 continue
      xsv2(l+1) = pcon
      qprof2(l+1) = qprof2(l)
      vp2(l+1) = vp2(l)
      xmja2(l+1) = xmja2(l)
      x2ave(l+1) = x2ave(l)
      gxmja2(l+1) = gxmja2(l)
      gja2(l+1) = gja2(l)
      vary(l+1) = vary(l)
      if(vary(l+1).le.0) ineg=14
      avhyp(l+1) = 0._R8
      term1aj(l+1) = term1aj(l)
      bpolar(l+1) = bpolar(l)
      rminora(l+1) = rminora(l)
      rmajora(l+1) = rmajora(l)
      elonga(l+1) = elonga(l)
      deltaa(l+1) = deltaa(l)
      ajave(l+1) = ajave(l)
      ajbtsq(l+1) = ajbtsq(l)
      ajbsq(l+1) = ajbsq(l)
      avepres(l+1) = avepres(l)
  600 continue
      if(iqtrubmax.le.2) then
      nskipsf = min(nskipsf+1,nskipsfi)
             endif
      if(iqtrubmax.ge.4 ) then
      nskipsf = max(nskipsf-1,1)
             endif
!     write(nterm,1077) kcycle,  iqtrubmax,nskipsf
 1077 format("modify: ",3i6)
      if(npsit.le.3) ineg=14
      if(ineg.ne.0) return
!
      if(npsit .ge. i95) go to 618
      do 617 k=1,kmax+3
      xplot(4,k) = xplot(1,k)
      zplot(4,k) = zplot(1,k)
  617 continue
      do 616 k=1,kmaxo+3
      xplot(3,k) = xplot(2,k)
      zplot(3,k) = zplot(2,k)
  616 continue
      kkm95 = kmax
      kkm95m = kmaxo
  618 continue
      if(npsit .ge. i90) go to 718
      do 716 k=1,kmax+3
      xplot(6,k) = xplot(1,k)
      zplot(6,k) = zplot(1,k)
  716 continue
      do 717 k=1,kmaxo+3
      xplot(5,k) = xplot(2,k)
      zplot(5,k) = zplot(2,k)
  717 continue
      kkm90 = kmax
      kkm90m = kmaxo
  718 continue
!
      vp2(1) = 2._R8*vp2(2) - vp2(3)
      xmja2(1) = 2._R8*xmja2(2) - xmja2(3)
      gxmja2(1) = 0.0_R8
      gja2(1) = 0._R8
      x2ave(1) = xmag**2
      avhyp(1) = 0._R8
      vary(1) = 1.E-6_R8
      ftrap(1) = 0.0_R8
      ajave(1) = 2._R8*ajave(2) - ajave(3)
      ajbtsq(1) = 2._R8*ajbtsq(2) - ajbtsq(3)
      ajbsq(1) = 2._R8*ajbsq(2) - ajbsq(3)
!
!
!
!
      if(vp2(1) .lt. 0) vp2(1) = 0._R8
      if(xmja2(1).lt..5_R8*xmja2(2)) xmja2(1) = .5_R8*xmja2(2)
      xsv2(npsi+1) = 2._R8*xsv2(npsi) - xsv2(npsi-1)
      vp2(npsi+1) = 2._R8*vp2(npsi) - vp2(npsi-1)
      gxmja2(npsi+1) = 2._R8*gxmja2(npsi) - gxmja2(npsi-1)
      gja2(npsi+1) = 2._R8*gja2(npsi) - gja2(npsi-1)
      xmja2(npsi+1) = 2._R8*xmja2(npsi) - xmja2(npsi-1)
!
      if(itype.eq.1) go to 603
      ajpress = 0._R8
      do 620 l=2,npsi+1
      xsv(l) = .5_R8*(xsv2(l)+xsv2(l-1))
      vp(l) = .5_R8*(vp2(l)+vp2(l-1))
      vpg(l) = vp(l)**(5._R8/3._R8)
      gxmja(l) = .5_R8*(gxmja2(l)+gxmja2(l-1))
      xmja(l) = .5_R8*(xmja2(l)+xmja2(l-1))
!
!.....added 2/14/03 .... toroidal current due to diamagnetic and Pfirsch-Schluter effects
      call peval(xsv(l),2,pval,ppval,imag,jmag)
      ajpress = ajpress + term1aj(l)*ppval
      ajpary(l)= term1aj(l)*ppval
  620 continue
      xsv(1) = xsv2(1)
      vp(1) = vp(2)
      vpg(1) = vp(1)**(5._R8/3._R8)
      gxmja(1) = gxmja(2)
      xmja(1) = xmja(2)
      if(iwall.eq.0.and.ibalsw.eq.0.and.irippl.eq.0) go to 605
      do 604 i=1,nthe+3
      xw(i,1) = xmag
      zw(i,1) = zmag
      xw(i,npsit+1) = 2._R8*xw(i,npsit)-xw(i,npsit-1)
  604 zw(i,npsit+1) = 2._R8*zw(i,npsit)-zw(i,npsit-1)
      call metricw
  605 continue
      psis1 = xsv2(npsit)
      psis2 = xsv2(npsit-1)
      psis3 = psilim
!
!.....define initial adiabatic arrays for pressure,electron pressure
!.....density and q
      if(idefnpe .eq. 1) then
         do l=2,npsi
            adp(l) = avepres(l)*vpg(l)
            ade(l) = acoef(2)*avepres(l)*vpg(l)
            call reval(xsv(l),idens,0,rval,rpval,imag,jmag)
            adn(l) = rval*vp(l)
            adi(l) = 1._R8/(.5_R8*(qprof2(l)+qprof2(l-1)))
         enddo
         adp(1) = adp(2)
         ade(1) = ade(2)
         adn(1) = adn(2)
         adi(1) = adi(2)
         idefnpe = 0
      endif
!
      if(kcycle.gt.0 .or. irst1.eq.2) go to 607
      do 602 l=2,npsi
      call peval(xsv(l),1,pval,ppval,imag,jmag)
      adp(l) = pval*vpg(l)
      call eeval(xsv(l),1,eval,epval,imag,jmag)
      ade(l) = eval*vpg(l)
      call reval(xsv(l),idens,0,rval,rpval,imag,jmag)
      adn(l) = rval*vp(l)
  602 adi(l) = 1._R8/(.5_R8*(qprof2(l)+qprof2(l-1)))
      adp(1) = adp(2)
      ade(1) = ade(2)
      adn(1) = adn(2)
      adi(1) = adi(2)
!
!.....boundary conditions for transport
  607 continue
      do 1607 l=2,npsit
 1607 adi(l) = 2._R8/(qprof2(l)+qprof2(l-1))
      do 608 l=npsit+1,npsi+1
      adn(l) = fracn0*r0*vp(l)
      ade(l) = smallt*fracn0*r0*vpg(l)
      adp(l) = acoef(882)*ade(l)
      adi(l) = 2._R8/(qprof2(l)+qprof2(l-1))
  608 continue
      if(npsit .ge. (npsi-1)) go to 610
      call defcoef
      if(itype.eq.3) return
!
      if(ibalsw.eq.2 .and. kcycle.ge.0) call balloon
      if(ii66.ne.1) then
      if( numargs .lt. 1 ) then
         filename = "fort66"
      else
         filename = "fort66" // '.' // trim(suffix)
      end if
      open(66,file=trim(filename),status='unknown',iostat=ios66)
      ii66 = 1
      endif
!
      write(66,6012) (idn(j),j=1,50)
 6012 format(60i2)
      if(ibalsw.eq.1     .and.                                           &  
     &  (kcycle.eq.0 .or. iplt2 .le. nskip2+1-iskipsf)) call balloon
!
!.....calculate betapol,ali2,vol,uint
      vol = vary(npsit)
      arad = sqrt(vol/(tpi*pi*xmag))
      sum1 = 0._R8
      sum2 = 0._R8
      sum3 = 0._R8
      sum4 = 0._R8
      sum5 = 0._R8
      sum6 = 0._R8
      sum7 = 0._R8
      sum8 = 0._R8
      sum9 = 0._R8
      do 609 l=2,npsit
!
      call peval(xsv(l),2,pval,ppval,imag,jmag)
      call reval(xsv(l),idens,1,rval,rpval,imag,jmag)
      call geval(xsv(l),2,gval,gpval,gppval,imag,jmag)
      sum7 = sum7 + anhe(l)*vp(l)*dpsi
      sum8 = sum8 + alphade(l)*vp(l)*dpsi
      sum9 = sum9 + alphapr(l)*vp(l)*dpsi
      sum1 = sum1 + pval*vp(l)*dpsi
      sum3 = sum3 + vp(l)*dpsi
      sum4 = sum4 + rval*vp(l)*dpsi
!
      sum6 = sum6 + anhy(l)*vp(l)*dpsi
      sum5 = sum5 + .5_R8*gval**2*xmja(l)*dpsi/                          &  
     &     (tpi*.5_R8*(qprof2(l-1)+qprof2(l)))
  609 sum2 = sum2 + gxmja(l)*dpsi/(tpi*.5_R8*(qprof2(l-1)+qprof2(l)))
      sum2 = sum2 + gxmja2(npsit)*(psilim-psis1)
      uint = 1.5_R8*sum1*udsi
!
      enermjo = enermj
      enermj = uint*1.E-6_R8
      dtfluxv = times - toldfv
      if(toldfv .ge. times .or. enermjo .le. 0) go to 6123
      wdotmw = (enermj - enermjo)/dtfluxv
 6123 continue
      toldfv = times
!
!
      apld = tcurdtp*tpi
      if(rminor .le. 0._R8 .or. rmajor .le. 0._R8) then
      call shapenp( 4,xshape,fshape,1.E-8_R8,ifail)
      if(ifail .ne. 0)       ineg=15
      if(ineg.ne.0) return
      rminor = xshape(2)
      rmajor = xshape(1)
      if(rminor .ne. 0._R8) shape5 = ellmom/rminor
        endif
      betapol = 4._R8*pi**2*rminor**2*(1._R8+shape5**2)*sum1/(sum3*      &  
     & apld**2)
      ali2 = sum2/(rmajor*apld**2)
      betat = sum1/sum5
      sum10 = .5_R8*sum3*gzero**2/xplas**2
      betat0 = sum1/sum10
      rvolav = sum4/sum3*udsd
      rvolavi = sum6/sum3
      rvolavh = sum7/sum3
      alphades = sum8/sum3
      alphabes = sum9/sum10
      rtot = sum4*udsd
      aliga = sum2*perim**2/(sum3*apld**2)
!
!.....DEBUG
!     write(nout,8088) sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10
!8088 format(" sums from flxvol",1p10e12.4)
!     write(nout,8089) npsi,npsit
!8089 format("npsi,npsit=",2i5)
!     do j=1,npsit
!     write(nout,8090) j,adn(j),adp(j),ade(j),adi(j),vp(j),qprof2(j)
!     enddo
!8090 format(i4,1p10e12.4)
!
      if(npsit .lt. (npsi-1)) return
      npsi = npsi+1
      if(npsi .lt. ppsi-1) return
  610 write(nout,1069)
 1069 format(" error, toroidal flux domain too small, increase tfmult")
      ineg=16
!
      return
!
!.....appendix section to calculate q profile on equal toroidal
!.....flux grid for itype=1    . . .   used in initialization only
  603 continue
      tftemp(1) = 0._R8
      do 700 l=2,npsi
      qave = .5_R8*(qprof2(l-1)+qprof2(l))
      tftemp(l) = tftemp(l-1) + qave*(xsv2(l)-xsv2(l-1))*tpi
  700 continue
      if(tfmax.eq.0 .or. irst1.ne.2) tfmax = tftemp(npsi)*tfmult
      lval = 2
      dpsi = tfmax/(npsi-1)
      write(nout,6666) dpsi
 6666 format(" dpsi =",1pe12.4,"  webers")
!
      rdpsi = 1._R8/dpsi
      do 710 l=2,npsi
      tfcon = (l-1.5_R8)*dpsi
      do 720 ll=lval,npsi
      llsave = ll
      if(tftemp(ll).gt.tfcon) go to 730
  720 continue
      adi(l) = 1._R8/qprof2(npsi)
      go to 740
  730 continue
      lval = llsave
      adi(l) = 1._R8/(qprof2(llsave-1) + (tfcon-tftemp(llsave-1))        &  
     &        /(tftemp(llsave)-tftemp(llsave-1))                         &  
     &      *(qprof2(llsave)-qprof2(llsave-1)))
  740 continue
  710 continue
      adi(1) = adi(2)
!
!
!.....special stability calculation for FRC
      if(kcycle.gt.0 .or. acoef(894).le.0) return
      write(nterm,1900)
 1900 format("   i    psi         vp         vol         vpg",           &  
     &      "            p       vpg*p       d(vpg*p)")
      vol = 0._R8
      pvpgtempo = 0._R8
      call peval(psimin,1,paxis,ppaxis,imag,jmag)
      prel = acoef(894)
      pboundary = paxis*prel/(1._R8-prel)
      do 901 l=2,npsit
      ps = xsv2(l)
      vptemp = vp2(l)*tpi*qprof2(l)
      vol =vol +  (xsv2(l)-xsv2(l-1))*vptemp
      vpgtemp = vptemp**(5._R8/3._R8)
      call peval(ps,1,preturn,ppval,imag,jmag)
      pval = preturn + pboundary
      pvpgtemp = pval*vpgtemp
      dpvg = pvpgtemp - pvpgtempo
      pvpgtempo = pvpgtemp
      write(nterm,1901) l,ps,vptemp,vol,vpgtemp,pval,pvpgtemp,dpvg
 1901 format(1x,i3,1p7e12.4)
  901 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
