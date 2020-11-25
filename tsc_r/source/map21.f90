      subroutine map21(sary2d,sary1d)
!
!.....uses 2d density distribution SARY2D(PNX,PNZ)  (cell centered)
!     to define differential density flux aray SARY1D(PPSI) (cell centered)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,i,l,i1,i2,i3,indx
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sary1d,sary2d,sumold,sumvolo,ps,sumvol,p1
      REAL*8 p2,p3,p4,x1,x2,x3,x4,z1,z2,z3,z4,at,at1save,sgn,area1
      REAL*8 xa,zb,za,xb,at2save,area2,frac,dvol,volpls
      REAL*8 sum
!============
!     dimension iskip(pnx,pnz),
      dimension sary1d(ppsi),sary2d(pnx,pnz)
!     dimension xsvm(ppsi), xsv2m(ppsi)
!============      
      INTEGER :: istat = 0 
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iskip
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xsvm
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xsv2m
!============      
      IF(.not.ALLOCATED(iskip)) ALLOCATE( iskip(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(xsvm)) ALLOCATE( xsvm(ppsi), STAT=istat)
      IF(.not.ALLOCATED(xsv2m)) ALLOCATE( xsv2m(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : map21  ' 
!============      
!
      do 105 j = 1,npsi
      xsvm(j)  = xsv (j) - (xsv2(1)-psimin)
      xsv2m(j) = xsv2(j) - (xsv2(1)-psimin)
      sary1d(j)= 0._R8
  105 continue
!
      do 490 i=iminn,imaxx
      do 490 j=jminn,jmaxx
  490 iskip(i,j) = 0
      sumold = 0._R8
      sumvolo = 0._R8
      do 500 l=2,npsit
      ps = xsv2m(l)
      sum = 0._R8
      sumvol = 0._R8
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
      frac = .50_R8*(area1/at1save + area2/at2save)
!
      if(frac.gt. 0.999_R8) iskip(i,j)=1
      go to 211
  210 frac=1._R8
      area1 = .5_R8*deex*deez
      area2 = .5_R8*deex*deez
  211 continue
      dvol = tpi*xary(i)*deex*deez
      sum = sum+frac*sary2d(i,j)*dvol
      sumvol = sumvol + frac*dvol
  300 continue
  400 continue
      volpls = (1._R8+isym)*(sumvol-sumvolo)
      if(volpls.eq.0.0_R8) volpls = vp(l)*dpsi
      sumvolo = sumvol
      sary1d(l) = (1._R8+isym)*(sum-sumold)/volpls
      sumold = sum
  500 continue
      sary1d(1) = sary1d(2)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
