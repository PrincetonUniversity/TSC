!#include "f77_dcomplx.h"
      subroutine metricw
!.....9.10
!
!
      USE CLINAM
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,nthmax,j,i1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xmax2,zmax2,xmin2,zmin2,fac,xmaxn,xminn,zmaxn,zminn,xx
      REAL*8 zz,fmh,gmh,gpmh,gppval,xta,xpa,xxa,xtb,xpb,xxb,xt5,xp5
      REAL*8 xx5,zz5,zta,zpa,ztb,zpb,zt5,zp5,aja,ajb,bsqfc
      REAL*8 AREAL
!============
      nthep2 = nthe+2
      nthep3 = nthe+3
!
      if(kcycle.lt.0.or.iplt.le.nskipl+1-iskipsf.or.iwall.eq.0)go to     &  
     & 501
      xmax2=-1.E20_R8
      zmax2=-1.E20_R8
      xmin2=1.E20_R8
      zmin2=1.E20_R8
      do 16 i=3,nthe+2
      if(xw(i,npsit).gt.xmax2) xmax2 = xw(i,npsit)
      if(xw(i,npsit).lt.xmin2) xmin2 = xw(i,npsit)
      if(zw(i,npsit).gt.zmax2) zmax2 = zw(i,npsit)
      if(zw(i,npsit).lt.zmin2) zmin2 = zw(i,npsit)
   16 continue
      fac = (xmax2-xmin2)/(zmax2-zmin2)
      if(fac.gt.1) go to 20
      xmaxn = .5_R8*(xmax2+xmin2)+.5_R8*(xmax2-xmin2)/fac
      xminn = .5_R8*(xmax2+xmin2)-.5_R8*(xmax2-xmin2)/fac
      xmax2 = xmaxn
      xmin2 = xminn
      go to 30
   20 continue
      zmaxn = .5_R8*(zmax2+zmin2)+.5_R8*(zmax2-zmin2)*fac
      zminn = .5_R8*(zmax2+zmin2)-.5_R8*(zmax2-zmin2)*fac
      zmax2 = zmaxn
      zmin2 = zminn
   30 continue
      call maps(xmin2,xmax2,zmin2,zmax2,.142_R8,.858_R8,.285_R8,1._R8)
!
!...draw flux contours inside plasma
!
      nthmax = nthe+2
      do 100 j=2,npsit
      i = 2
      xx = xw(i,j)
      zz = zw(i,j)
      call setcrt(xx,zz)
      do 99 i=3,nthmax
      i1 = i
      xx = xw(i1,j)
      zz = zw(i1,j)
      call vector(xx,zz)
   99 continue
  100 continue
      do 200 i=2,nthe+1+isym
      j = 1
      i1 = i
      xx = xw(i1,j)
      zz = zw(i1,j)
      call setcrt(xx,zz)
      do 199 j=2,npsit
      xx = xw(i1,j)
      zz = zw(i1,j)
      call vector(xx,zz)
  199 continue
  200 continue
      call frscj(10)
  501 continue
      if(iwall .eq. 0 .and. irippl .eq. 0) return
!
!.....calculate metric coefficient arrays
      dthe = (2._R8-AREAL(isym))*pi/nthe
      do 400 j=2,npsit
!
      fmh = (xsv2(j)-xsv2(j-1))*rdpsi
      call geval(xsv(j),2,gmh,gpmh,gppval,imag,jmag)
      do 300 i=2,nthe+2
      xta = (xw(i+1,j)+xw(i+1,j-1)-xw(i-1,j)-xw(i-1,j-1))*.25_R8/dthe
      xpa = (xw(i,j)-xw(i,j-1))/dpsi
      xxa = (xw(i,j)+xw(i,j-1))*.5_R8
      xtb = (xw(i,j)-xw(i-1,j))/dthe
      xpb = (xw(i,j+1)+xw(i-1,j+1)-xw(i,j-1)-xw(i-1,j-1))*.25_R8/dpsi
      xxb = (xw(i,j)+xw(i-1,j))*.5_R8
      xt5 = (xw(i,j)+xw(i,j-1)-xw(i-1,j)-xw(i-1,j-1))*.5_R8/dthe
      xp5 = (xw(i,j)+xw(i-1,j)-xw(i,j-1)-xw(i-1,j-1))*.5_R8/dpsi
      xx5 = (xw(i,j)+xw(i,j-1)+xw(i-1,j)+xw(i-1,j-1))*.25_R8
      zz5 = (zw(i,j)+zw(i,j-1)+zw(i-1,j)+zw(i-1,j-1))*.25_R8
      zta = (zw(i+1,j)+zw(i+1,j-1)-zw(i-1,j)-zw(i-1,j-1))*.25_R8/dthe
      zpa = (zw(i,j)-zw(i,j-1))/dpsi
      ztb = (zw(i,j)-zw(i-1,j))/dthe
      zpb = (zw(i,j+1)+zw(i-1,j+1)-zw(i,j-1)-zw(i-1,j-1))*.25_R8/dpsi
      zt5 = (zw(i,j)+zw(i,j-1)-zw(i-1,j)-zw(i-1,j-1))*.5_R8/dthe
      zp5 = (zw(i,j)+zw(i-1,j)-zw(i,j-1)-zw(i-1,j-1))*.5_R8/dpsi
!
      aja =-xxa*(xpa*zta-xta*zpa)
      ajb =-xxb*(xpb*ztb-xtb*zpb)
      sqg(i,j) =-xx5*(xp5*zt5-xt5*zp5)
      if(sqg(i,j).le.0.0_R8) ineg=22
      if(aja.le.0.0_R8) ineg=22
      if(ajb.le.0.0_R8) ineg=22
       if(ineg.ne.0) go to 400
!
      bsqfc = (fmh**2*(xt5**2+zt5**2)*(xx5/sqg(i,j))**2                  &  
     &         + gmh**2)/xx5**2
      xcentfc(i,j) = xx5
      zcentfc(i,j) = zz5
      bmagfc(i,j)  = sqrt(bsqfc)
      if(iwall.eq.0) go to 300
!...note that vacuum toroidal field is used here
      g33(i,j) = aja/((xta**2+zta**2) + (aja*gzero/xxa)**2)
      g11(i,j) = xxb**2*(xtb**2+ztb**2)/ajb
      g22(i,j) = xxa**2*(xpa**2+zpa**2)/aja
      g21(i,j) = -xxa**2*(xta*xpa+zta*zpa)/aja
      g12(i,j) = -xxb**2*(xtb*xpb+ztb*zpb)/ajb
  300 continue
      if(isym.eq.1) go to 301
      sqg(nthe+3,j) = sqg(3,j)
      bmagfc(nthe+3,j) = bmagfc(3,j)
      xcentfc(nthe+3,j) = xcentfc(3,j)
      zcentfc(nthe+3,j) = zcentfc(3,j)
      if(iwall.eq.0) go to 301
      g33(nthe+3,j) = g33(3,j)
      g11(nthe+3,j) = g11(3,j)
      g22(nthe+3,j) = g22(3,j)
      g21(nthe+3,j) = g21(3,j)
      g12(nthe+3,j) = g12(3,j)
      go to 400
  301 continue
      sqg(nthe+3,j) = sqg(nthe+2,j)
      bmagfc(nthe+3,j) = bmagfc(nthe+2,j)
      xcentfc(nthe+3,j) = xcentfc(nthe+2,j)
      zcentfc(nthe+3,j) = zcentfc(nthe+2,j)
      if(iwall.eq.0) go to 400
      g33(nthe+3,j) = g33(nthe+1,j)
      g11(nthe+3,j) = g11(nthe+2,j)
      g22(nthe+3,j) = g22(nthe+1,j)
      g21(nthe+3,j) =-g21(nthe+1,j)
      g12(nthe+3,j) =-g12(nthe+2,j)
  400 continue
      do 450 i=2,nthe+3
      bmagfc(i,1) = bmagfc(i,2)
      bmagfc(i,npsit+1) = 2._R8*bmagfc(i,npsit)-bmagfc(i,npsit-1)
      xcentfc(i,1) = xmag
      zcentfc(i,1) = zmag
      xcentfc(i,npsit+1) = 2._R8*xcentfc(i,npsit)-xcentfc(i,npsit-1)
      zcentfc(i,npsit+1) = 2._R8*zcentfc(i,npsit)-zcentfc(i,npsit-1)
      sqg(i,1) = 0._R8
      sqg(i,npsit+1) = 2._R8*sqg(i,npsit)-sqg(i,npsit-1)
      if(sqg(i,npsit+1).le.0.0_R8) ineg=22
      if(iwall.eq.0) go to 450
      g33(i,1) = 0._R8
      g11(i,1) = 0._R8
      g22(i,1) = 0._R8
      g21(i,1) = 0._R8
      g12(i,1) = 0._R8
      g33(i,npsit+1) = 2._R8*g33(i,npsit)-g33(i,npsit-1)
      g11(i,npsit+1) = 2._R8*g11(i,npsit)-g11(i,npsit-1)
      g22(i,npsit+1) = 2._R8*g22(i,npsit)-g22(i,npsit-1)
      g21(i,npsit+1) = 2._R8*g21(i,npsit)-g21(i,npsit-1)
      g12(i,npsit+1) = 2._R8*g12(i,npsit)-g12(i,npsit-1)
  450 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
