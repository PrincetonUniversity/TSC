      subroutine qcalc
!
!......this subroutine calculates the q-profile by doing a volume integr
!......of the toroidal flux inside a given poloidal flux surface.  it is
!......mainly used when isurf=0, ie, non-surface averaged mode.  it also
!......calculates beta-poloidal and li/2.
!
!
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j,n,i1,i2,i3,indx
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 psiinc,enpl,enbp,envol,ps,sumt,sumv,sumg,p1,p2,p3,p4
      REAL*8 x1,x2,x3,x4,z1,z2,z3,z4,at,at1save,sgn,area1,xa,zb,za
      REAL*8 xb,at2save,area2,area,tzone,psimid,pval,ppval,da,rmi
      REAL*8 dv,dpdxs,dpdzs,bsq,tcuroc
      REAL*8 sum
!============
!     dimension iskip(pnx,pnz)
!============      
      INTEGER :: istat = 0 
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iskip
!============      
      IF(.not.ALLOCATED(iskip)) ALLOCATE( iskip(pnx,pnz), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : qcalc  ' 
!============      
!
      psiinc = (psilim-psimin)/npts
      polflx(1) = tpi*psimin
      torflx(1) = 0._R8
      temflx(1) = 0._R8
      gemflx(1) = 0._R8
      volflx(1) = 0._R8
      enpl      = 0._R8
      enbp      = 0._R8
      envol     = 0._R8
      helic = 0
      do 490 i=iminn,imaxx
      do 490 j=jminn,jmaxx
  490 iskip(i,j) = 0
      do 500 n=1,npts
      ps = psimin + n*psiinc
      polflx(n+1) = tpi*ps
      sum = 0._R8
      sumt = 0._R8
      sumv = 0._R8
      sumg = 0._R8
      do 400 i=iminn,imaxx
      do 300 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 300
      if(iskip(i,j).eq.1) go to 210
!
!.....calculate area of triangular sections with psi lt ps
      p1 = psi(i+1,j+1)
      p2 = psi(i,j)
      p3 = psi(i+1,j)
      p4 = psi(i,j+1)
      if(p1.gt.ps.and.p2.gt.ps.and.p3.gt.ps.and.p4.gt.ps) go to 300
      x1 = xary(i+1)
      x2 = xary(i)
      x3 = xary(i+1)
      x4 = xary(i)
      z1 = zary(j+1)
      z2 = zary(j)
      z3 = zary(j)
      z4 = zary(j+1)
!
!
!.....area of first triangle
      at = .5_R8*(z1-z3)*(x3-x2)
      at1save = at
      sgn = -1._R8
      i1 = 0
      i2 = 0
      i3 = 0
      if(p1.gt.ps) i1 = 1
      if(p2.gt.ps) i2 = 1
      if(p3.gt.ps) i3 = 1
      indx = i1 + 2*i2 + 4*i3 + 1
      go to(10,20,30,35,40,25,15,50),indx
   10 continue
      area1 = at
      go to 100
   15 at = 0._R8
      sgn = 1._R8
   20 xa = ((p1-ps)*x2+(ps-p2)*x1)/(p1-p2)
      zb = ((p1-ps)*z3+(ps-p3)*z1)/(p1-p3)
      area1 = at+.5_R8*(x1-xa)*(z1-zb)*sgn
      go to 100
   25 at = 0
      sgn = 1._R8
   30 za = ((p1-ps)*z2+(ps-p2)*z1)/(p1-p2)
      xb = ((p3-ps)*x2+(ps-p2)*x3)/(p3-p2)
      area1 = at+.5_R8*(za-z2)*(xb-x2)*sgn
      go to 100
   35 at = 0
      sgn = 1._R8
   40 zb = ((p1-ps)*z3+(ps-p3)*z1)/(p1-p3)
      xa = ((p3-ps)*x2+(ps-p2)*x3)/(p3-p2)
      area1 = at+.5_R8*(zb-z3)*(x3-xa)*sgn
      go to 100
   50 continue
      area1 = 0._R8
  100 continue
!
!.....area of second triangle
      at = .5_R8*(z4-z2)*(x1-x4)
      at2save = at
      sgn = -1._R8
      i1 = 0._R8
      i2 = 0
      i3 = 0
      if(p2.gt.ps) i1 = 1
      if(p1.gt.ps) i2 = 1
      if(p4.gt.ps) i3 = 1
      indx = i1 + 2*i2 + 4*i3 + 1
      go to(60,70,80,85,90,75,65,95),indx
   60 area2 = at
      go to 200
   65 at = 0._R8
      sgn = 1._R8
   70 xa = ((p1-ps)*x2+(ps-p2)*x1)/(p1-p2)
      zb = ((p4-ps)*z2+(ps-p2)*z4)/(p4-p2)
      area2 = at+.5_R8*(zb-z2)*(xa-x2)*sgn
      go to 200
   75 at = 0._R8
      sgn = 1._R8
   80 za = ((p1-ps)*z2+(ps-p2)*z1)/(p1-p2)
      xb = ((p1-ps)*x4+(ps-p4)*x1)/(p1-p4)
      area2 = at+.5_R8*(x1-xb)*(z1-za)*sgn
      go to 200
   85 at = 0
      sgn = 1._R8
   90 xa = ((p1-ps)*x4+(ps-p4)*x1)/(p1-p4)
      zb = ((p4-ps)*z2+(ps-p2)*z4)/(p4-p2)
      area2 = at+.5_R8*(xa-x4)*(z4-zb)*sgn
      go to 200
   95 area2 = 0
  200 continue
      area = .50_R8*(area1/at1save + area2/at2save)
!
      if(area.gt. 0.999_R8) iskip(i,j)=1
      go to 211
  210 area=1._R8
      area1 = .5_R8*deex*deez
      area2 = .5_R8*deex*deez
  211 continue
      tzone = ((.5_R8*udsh)/udsd)*pr(i+1,j+1)*ajey(i+1)                  &  
     &        /roj(i+1,j+1)
      sumt = sumt+tzone*area
      sumv = sumv+area*ajey(i+1)
      sum = sum+area*g(i+1,j+1)
      sumg = sumg + g(i+1,j+1)*xsqoj(i+1)*ajey(i+1)*area
!
      if(n.ne.npts) goto 300
      psimid=.25_R8*(psi(i+1,j+1)+psi(i,j)+psi(i+1,j)+psi(i,j+1))
      call peval(psimid,1,pval,ppval,i,j)
      da=(area1+area2)
      rmi=xarh(i+1)
      dv=tpi*rmi*da
      enpl=enpl+pval*dv
      dpdxs=.5_R8*((psi(i+1,j)-psi(i,j))**2+(psi(i+1,j+1)-psi(i,j+1))**  &  
     & 2)
      dpdzs=.5_R8*((psi(i,j+1)-psi(i,j))**2+(psi(i+1,j+1)-psi(i+1,j))**  &  
     & 2)
      bsq=(dpdxs/deex**2+dpdzs/deez**2)/rmi**2
      enbp=enbp+bsq*dv
      envol=envol+dv
  300 continue
  400 continue
      torflx(n+1) = (1+isym)*sum
      qprof(n)=(torflx(n+1)-torflx(n))/(polflx(n+1)-polflx(n))
      temflx(n+1) = (1+isym)*sumt
      volflx(n+1) = (1+isym)*sumv
      if(volflx(n).ne.volflx(n+1))                                       &  
     &tprof(n) = (temflx(n+1)-temflx(n))/(volflx(n+1)-volflx(n))
      gemflx(n+1) = (1+isym)*sumg
      if(volflx(n).ne.volflx(n+1))                                       &  
     &gprof(n) = (gemflx(n+1)-gemflx(n))/(volflx(n+1)-volflx(n))
      helic = helic-(polflx(n)+polflx(n+1))*(torflx(n+1)-torflx(n))
  500 continue
      qprof(1) = qprof(3)
      qprof(2) = qprof(3)
!
      tcuroc = tcuro
      if(apl.ne.0) tcuroc = apl*usdi
      betapol = 4._R8*pi**2*rminor**2*(1._R8+shape5**2)*enpl/(envol*     &  
     & tcuroc**2)
      ali2 = (1+isym)*enbp/(xmag*tcuroc**2)
      vol = (1+isym)*envol
      uint = 1.5_R8*(1+isym)*enpl*udsi
      if(acoef(504).ne.0) uint = 0.5_R8*(1+isym)*enbp*udsi
      helic = helic + 2._R8*polflx(npts+1)*torflx(npts+1)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
      subroutine acalc(ps,area)
!
!......this subroutine calculates the area within the plasma
!......and the area within the halo region
!
!
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j,n,i1,i2,i3,indx
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 area,ps,p1,p2,p3,p4
      REAL*8 x1,x2,x3,x4,z1,z2,z3,z4,at,at1save,sgn,area1,xa,zb,za
      REAL*8 xb,area2,da
      area = 0.
      do 400 i=iminn,imaxx
      do 300 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 300
!
!.....calculate area of triangular sections with psi lt ps
      p1 = psi(i+1,j+1)
      p2 = psi(i,j)
      p3 = psi(i+1,j)
      p4 = psi(i,j+1)
      if(p1.gt.ps.and.p2.gt.ps.and.p3.gt.ps.and.p4.gt.ps) go to 300
      x1 = xary(i+1)
      x2 = xary(i)
      x3 = xary(i+1)
      x4 = xary(i)
      z1 = zary(j+1)
      z2 = zary(j)
      z3 = zary(j)
      z4 = zary(j+1)
!
!
!.....area of first triangle
      at = .5_R8*(z1-z3)*(x3-x2)
      sgn = -1._R8
      i1 = 0
      i2 = 0
      i3 = 0
      if(p1.gt.ps) i1 = 1
      if(p2.gt.ps) i2 = 1
      if(p3.gt.ps) i3 = 1
      indx = i1 + 2*i2 + 4*i3 + 1
      go to(10,20,30,35,40,25,15,50),indx
   10 continue
      area1 = at
      go to 100
   15 at = 0._R8
      sgn = 1._R8
   20 xa = ((p1-ps)*x2+(ps-p2)*x1)/(p1-p2)
      zb = ((p1-ps)*z3+(ps-p3)*z1)/(p1-p3)
      area1 = at+.5_R8*(x1-xa)*(z1-zb)*sgn
      go to 100
   25 at = 0
      sgn = 1._R8
   30 za = ((p1-ps)*z2+(ps-p2)*z1)/(p1-p2)
      xb = ((p3-ps)*x2+(ps-p2)*x3)/(p3-p2)
      area1 = at+.5_R8*(za-z2)*(xb-x2)*sgn
      go to 100
   35 at = 0
      sgn = 1._R8
   40 zb = ((p1-ps)*z3+(ps-p3)*z1)/(p1-p3)
      xa = ((p3-ps)*x2+(ps-p2)*x3)/(p3-p2)
      area1 = at+.5_R8*(zb-z3)*(x3-xa)*sgn
      go to 100
   50 continue
      area1 = 0._R8
  100 continue
!
!.....area of second triangle
      at = .5_R8*(z4-z2)*(x1-x4)
      sgn = -1._R8
      i1 = 0._R8
      i2 = 0
      i3 = 0
      if(p2.gt.ps) i1 = 1
      if(p1.gt.ps) i2 = 1
      if(p4.gt.ps) i3 = 1
      indx = i1 + 2*i2 + 4*i3 + 1
      go to(60,70,80,85,90,75,65,95),indx
   60 area2 = at
      go to 200
   65 at = 0._R8
      sgn = 1._R8
   70 xa = ((p1-ps)*x2+(ps-p2)*x1)/(p1-p2)
      zb = ((p4-ps)*z2+(ps-p2)*z4)/(p4-p2)
      area2 = at+.5_R8*(zb-z2)*(xa-x2)*sgn
      go to 200
   75 at = 0._R8
      sgn = 1._R8
   80 za = ((p1-ps)*z2+(ps-p2)*z1)/(p1-p2)
      xb = ((p1-ps)*x4+(ps-p4)*x1)/(p1-p4)
      area2 = at+.5_R8*(x1-xb)*(z1-za)*sgn
      go to 200
   85 at = 0
      sgn = 1._R8
   90 xa = ((p1-ps)*x4+(ps-p4)*x1)/(p1-p4)
      zb = ((p4-ps)*z2+(ps-p2)*z4)/(p4-p2)
      area2 = at+.5_R8*(xa-x4)*(z4-zb)*sgn
      go to 200
   95 area2 = 0
  200 continue
!
!
      da=(area1+area2)
      area = area + da
  300 continue
  400 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
