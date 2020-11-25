!#include "f77_dcomplx.h"
      subroutine balloon
!
!     check all surfaces for mercier and balloon stability (ibal1.gt.0)
!
      USE CLINAM
      USE BALCLI
      USE WALLCL
      USE BALCOM1

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,i,ii,im
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rdthe2,rdpsi2,pval,gppval,count,fac,xe,xte
      REAL*8 zte,xpe,zpe,aje,xem,xtem,ztem,xpem,zpem,ajem,xep,xtep
      REAL*8 ztep,xpep,zpep,ajep,gps,qdp,dtp1,xsqdps,ajedp,gts,gpgt
      REAL*8 xsq,bsqd,coef1,coef2,term1,term2,term3,term4,coef3
      REAL*8 gamat,avg0,ava2i,avg1a2i,avg1,avg1sa2i
      REAL*8 AREAL
!============
!     common/balcom1/
!    .          bsqi(pnthe,ppsi),xbal(pnthe,ppsi),zbal(pnthe,ppsi),
!    1          delt(pnthe,ppsi),deltp(pnthe,ppsi),deltm(pnthe,ppsi),
!    2          delp1(pnthe,ppsi),delm1(pnthe,ppsi),delp(pnthe,ppsi),
!    3          del(pnthe,ppsi),
!    .          qpmh(ppsi),fpe(ppsi),fmh(ppsi),fph(ppsi),
!    1          fm3(ppsi),ppmh(ppsi),gmh(ppsi),gpmh(ppsi),
!    .          sum1(ppsi),sum2(ppsi),sum3(ppsi),sum4(ppsi),
!    1          sum5(ppsi),sum6(ppsi)
!
!
      if(nthe.le.0) then
      write(nterm,6001) nthe
      write(nout, 6001) nthe
 6001 format(" error in balloon, (type07 field 7)nthe= ",i5)
      ineg=31
      return
      endif
!
      dthe = (2._R8-AREAL(isym))*pi/nthe
      rdthe2 = .5_R8/dthe
      rdpsi2 = .5_R8*rdpsi
!
!.....initialize to zero
      do 80 j=1,npsi
      idi(j) = 0
      idr(j) = 0
      idn(j) = 0
      idf(j) = 0
      ibaloon(j) = 0
   80 continue
!
!
!     compute mercier criterion and store coefficients
!     for balloon euler equation
!
!
!.....computes averages over theta of things needed for di and dr
!     and similar things.
!
!.......first define new origin and sence for theta coordinate
      do 81 i=1,nthe+2+isym
      if(isym.eq.1) ii = nthe+3+isym-i
      if(isym.eq.0) ii = nthe/2 + 4 - i
      if(ii.le.0) ii = ii + nthe
      do 81 j=1,npsi
      xbal(i,j) = xw(ii,j)
      zbal(i,j) = zw(ii,j)
   81 continue
!
      do 82 j=3,npsit-2
      qpmh(j) = (qprof2(j+1)-qprof2(j-1))*rdpsi2
      fpe(j) = (xsv2(j+1)-2._R8*xsv2(j)+xsv2(j-1))*rdpsi**2
      fmh(j) = (xsv2(j+1)-xsv2(j-1))*rdpsi2
      fph(j) = (xsv2(j+1)-xsv2(j  ))*rdpsi
      fm3(j) = (xsv2(j  )-xsv2(j-1))*rdpsi
      call peval(xsv2(j),2,pval,ppmh(j),imag,jmag)
      call geval(xsv2(j),2,gmh(j),gpmh(j),gppval,imag,jmag)
   82 continue
!
      do 83 j=3,npsit-2
      sum1(j) = 0._R8
      sum2(j) = 0._R8
      sum3(j) = 0._R8
   83 continue
      count = 0._R8
      do 50 i=2,nthe+1+isym
      fac = 1._R8
      if(isym.eq.1 .and. (i.eq.2 .or. i.eq.nthe+2)) fac = .5_R8
      count = count + fac
      ii = i-1
      do 50 j=3,npsit-2
      xe = xbal(i,j)
      xte = (xbal(i+1,j) - xbal(i-1,j))*rdthe2
      zte = (zbal(i+1,j) - zbal(i-1,j))*rdthe2
      xpe = (xbal(i,j+1) - xbal(i,j-1))*rdpsi2
      zpe = (zbal(i,j+1) - zbal(i,j-1))*rdpsi2
      aje = -xe*(xte*zpe-xpe*zte)
      xem = .5_R8*(xbal(i,j-1) + xe)
      xtem = .5_R8*((xbal(i+1,j-1) - xbal(i-1,j-1))*rdthe2 + xte)
      ztem = .5_R8*((zbal(i+1,j-1) - zbal(i-1,j-1))*rdthe2 + zte)
      xpem = (xbal(i,j) - xbal(i,j-1))*rdpsi
      zpem = (zbal(i,j) - zbal(i,j-1))*rdpsi
      ajem =-xem*(xtem*zpem-xpem*ztem)
      xep = .5_R8*(xbal(i,j+1) + xe)
      xtep = .5_R8*((xbal(i+1,j+1) - xbal(i-1,j+1))*rdthe2 + xte)
      ztep = .5_R8*((zbal(i+1,j+1) - zbal(i-1,j+1))*rdthe2 + zte)
      xpep = (xbal(i,j+1) - xbal(i,j))*rdpsi
      zpep = (zbal(i,j+1) - zbal(i,j))*rdpsi
      ajep =-xep*(xtep*zpep-xpep*ztep)
      delt(ii,j) = aje/xe**2
      deltp(ii,j) = ajep/xep**2
      deltm(ii,j) = ajem/xem**2
      sum1(j) = sum1(j) + fac*delt(ii,j)
      sum2(j) = sum2(j) + fac*deltp(ii,j)
      sum3(j) = sum3(j) + fac*deltm(ii,j)
      gps = (xte**2 + zte**2)*(xe/aje)**2
      bsqi(i,j) = xe**2/(fmh(j)**2*gps + gmh(j)**2)
   50 continue
      if(isym.eq.0) go to 51
      do 84 j=3,npsit-2
      bsqi(1,j) = bsqi(3,j)
      bsqi(nthe+3,j) = bsqi(nthe+1,j)
   84 continue
      go to 52
   51 continue
      do 85 j=3,npsit-2
      bsqi(nthe+2,j) = bsqi(2,j)
      bsqi(1,j) = bsqi(nthe+1,j)
   85 continue
   52 continue
!
      do 86 j=3,npsit-2
      sum1(j) = sum1(j)/count
      sum2(j) = sum2(j)/count
      sum3(j) = sum3(j)/count
      del(1,j) = 0._R8
      delp1(1,j) = 0._R8
      delm1(1,j) = 0._R8
   86 continue
      do 60 i=2,nthe+1+isym
      ii = i-1
      im = ii-1
      do 87 j=3,npsit-2
      delt(ii,j) = delt(ii,j)/sum1(j) - 1._R8
      deltp(ii,j) = deltp(ii,j)/sum2(j) - 1._R8
      deltm(ii,j) = deltm(ii,j)/sum3(j) - 1._R8
   87 continue
      if(ii.eq.1) go to 61
      do 88 j=3,npsit-2
      del(ii,j) = del(im,j) + .5_R8*dthe*(delt(im,j)+delt(ii,j))
      delp1(ii,j) = delp1(im,j) + .5_R8*dthe*(deltp(im,j)+deltp(ii,j))
      delm1(ii,j) = delm1(im,j) + .5_R8*dthe*(deltm(im,j)+deltm(ii,j))
   88 continue
   61 continue
      do 60 j=3,npsit-2
      delp(ii,j) = (delp1(ii,j)-delm1(ii,j))*rdpsi
   60 continue
!
      do 89 j=3,npsit-2
      sum1(j) = 0._R8
      sum2(j) = 0._R8
      sum3(j) = 0._R8
      sum4(j) = 0._R8
      sum5(j) = 0._R8
   89 continue
!
      count = 0._R8
      do 100 i=2,nthe+1+isym
      fac = 1._R8
      if(isym.eq.1 .and. (i.eq.2 .or. i.eq.nthe+2)) fac = .5_R8
      count = count + fac
      ii = i-1
      do 91 j=3,npsit-2
      xe = xbal(i,j)
      xte = (xbal(i+1,j) - xbal(i-1,j))*rdthe2
      zte = (zbal(i+1,j) - zbal(i-1,j))*rdthe2
      xpe = (xbal(i,j+1) - xbal(i,j-1))*rdpsi2
      zpe = (zbal(i,j+1) - zbal(i,j-1))*rdpsi2
      aje =-xe*(xte*zpe-xpe*zte)
      xem = .5_R8*(xbal(i,j-1) + xe)
      xtem = .5_R8*((xbal(i+1,j-1) - xbal(i-1,j-1))*rdthe2 + xte)
      ztem = .5_R8*((zbal(i+1,j-1) - zbal(i-1,j-1))*rdthe2 + zte)
      xpem = (xbal(i,j) - xbal(i,j-1))*rdpsi
      zpem = (zbal(i,j) - zbal(i,j-1))*rdpsi
      ajem =-xem*(xtem*zpem-xpem*ztem)
      xep = .5_R8*(xbal(i,j+1) + xe)
      xtep = .5_R8*((xbal(i+1,j+1) - xbal(i-1,j+1))*rdthe2 + xte)
      ztep = .5_R8*((zbal(i+1,j+1) - zbal(i-1,j+1))*rdthe2 + zte)
      xpep = (xbal(i,j+1) - xbal(i,j))*rdpsi
      zpep = (zbal(i,j+1) - zbal(i,j))*rdpsi
      ajep =-xep*(xtep*zpep-xpep*ztep)
      qdp = qpmh(j)*del(ii,j) + qprof2(j)*delp(ii,j)
      dtp1 = delt(ii,j) + 1._R8
      xsqdps = 2._R8*xpe/xe
!
      ajedp = (ajep/fph(j) - ajem/fm3(j))*rdpsi/fmh(j)
!
      gps = (xte**2 + zte**2)*(xe/aje)**2
      gts = (xpe**2 + zpe**2)*(xe/aje)**2
      gpgt = -(xte*xpe + zte*zpe)*(xe/aje)**2
      xsq = xe**2
      bsqd = (bsqi(i+1,j)-bsqi(i-1,j))*rdthe2
 
      coef1 = fmh(j)*bsqi(i,j)/aje
!
      alfav(ii,j) = coef1*(1._R8/xsq+qprof2(j)**2*gts*dtp1**2+qdp**2*    &  
     & gps                                                               &  
     &                   + 2._R8*qdp*qprof2(j)*dtp1*gpgt)
      alfa1(ii,j) = coef1*(2._R8*qprof2(j)*qpmh(j)*gpgt*dtp1             &  
     &                     + 2._R8*qpmh(j)*qdp*gps)
      alfa2(ii,j) = coef1*(qpmh(j)**2*gps)
!
      coef2 = ppmh(j)*bsqi(i,j)
      betav(ii,j) = coef2*(gpgt*aje/xsq)
      beta1(ii,j) = 0._R8
!
      term1 = -ajedp*fmh(j)**2*gps*bsqi(i,j)/xsq
      term2 =  aje*bsqi(i,j)/fmh(j)*(ppmh(j)+gmh(j)*gpmh(j)/xsq)
      term3 = -xsqdps*aje*bsqi(i,j)*(gmh(j)**2/xsq)/fmh(j)**2
      term4 = gmh(j)*bsqd/fmh(j)*qdp
      coef3 = ppmh(j)
      gamat = coef3*gmh(j)*qpmh(j)/fmh(j)*bsqi(i,j)
      gamav(ii,j) = coef3*(term1+term2+term3+term4)
      gama1(ii,j) = coef3*gmh(j)*qpmh(j)*bsqd/fmh(j)
!
!
!
      sum1(j) = sum1(j) + fac*gamav(ii,j)
      sum2(j) = sum2(j) + fac/alfa2(ii,j)
      sum3(j) = sum3(j) + fac*gamat/alfa2(ii,j)
      sum4(j) = sum4(j) + fac*gamat
      sum5(j) = sum5(j) + fac*gamat**2/alfa2(ii,j)
   91 continue
!
!
  100 continue
!
!
      do 92 j=3,npsit-2
      avg0 = sum1(j)/count
      ava2i = sum2(j)/count
      avg1a2i = sum3(j)/count
      avg1 = sum4(j)/count
      avg1sa2i = sum5(j)/count
!
      di(j)=avg0*ava2i + avg1a2i-avg1*ava2i + ava2i*avg1sa2i-avg1a2i**2  &  
     &   - 0.25_R8
!
!
   92 continue
!
      ifcount = 0
!
!     integrate euler equation ten times around torus
!
      call intgrte
      do 110 j=3,npsit-2
!....no of nodes in ballooning eigenfunction
      idn(j) = node1(j)
!....no of functional evaluations
      idf(j) = ifcount
      if(di(j).gt.0) idi(j) = 1
!
      ibaloon(j) = 0
!
!.....mark integer q surfaces
      if(qprof2(j-1) .le. 1 .and. qprof2(j) .gt. 1) ibaloon(j) = 6
      if(qprof2(j-1) .le. 2 .and. qprof2(j) .gt. 2) ibaloon(j) = 7
      if(qprof2(j-1) .le. 3 .and. qprof2(j) .gt. 3) ibaloon(j) = 8
      if(idn(j) .gt. 0) ibaloon(j) = 1
      if(idi(j) .gt. 0) ibaloon(j) = ibaloon(j) + 2
  110 continue
      ibaloon(npsit) = 5
      do 120 j=1,npsit-2
      idn(j) = idn(j) + idi(j)
  120 continue
!
!...adjust values around orign
      if(idn(5).gt.idn(4)) idn(4) = idn(5)
      if(idn(4).gt.idn(3)) idn(3) = idn(4)
      idn(2) = idn(3)
      idn(1) = idn(3)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
