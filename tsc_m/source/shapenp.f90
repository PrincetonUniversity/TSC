      subroutine shapenp(nshape,xshape,fshape,tol,ifail)
!.....3.94 shapenp
!....neil'l replacement for the nag routine c05nbf to obtain shape parameters
!
!
      USE CLINAM
      USE GEOPAR

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.... common block shared with rawmeas:
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ifail,nshape,kk,kkmax,i,itop,ibot
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xshape,fshape,tol,xra,xla
      REAL*8 xtopa,ztopa,xbota,zbota,xlb,denom,denom2,denom3,am2
      REAL*8 am3,ams,bcoef,ccoef,disc,xcrit,ycrit,ydenom,ynorm,f1
      REAL*8 xlb99,elong,center,hwidth,arg,triang,f2,xl95,xlb95
      REAL*8 xr95,xbot95,xtop95,zbot95,ztop95,cent95,hw95,arg95
      REAL*8 xl90,xlb90,xr90,xbot90,xtop90,zbot90,ztop90,cent90
      REAL*8 hw90,arg90
!============
!     common /geopar/ xbot,xtop,zbot,ztop,xr,xl
!
!
      dimension xshape(nshape),fshape(nshape)
      dimension xra(6),xla(6),xtopa(6),ztopa(6),xbota(6),zbota(6),xlb(6)
!============      
      INTEGER :: istat = 0 
!============      
      ifail = 0
!
!.....look for surfaces bordering 95% surface
!
!.....calculate shape values at 6 flux surfaces
      do 301 kk=1,6
!     write(nout,9911) kk, kcycle
 9911 format(" shapenp, kk, kcycle",2i5)
      kkmax = kmax
      if(kk.eq.2) kkmax = kmaxo
      if(kk.eq.3) kkmax = kkm95m
      if(kk.eq.4) kkmax = kkm95
      if(kk.eq.5) kkmax = kkm90m
      if(kk.eq.6) kkmax = kkm90
      if(isym.eq.0) go to 1
      xra(kk)=xplot(kk,kkmax)
      xla(kk)=xplot(kk,2)
      go to 3
    1 continue
      xra(kk) = xplot(kk,1)
      xla(kk) = xplot(kk,2)
      do 4 i=1,kkmax
      xra(kk) = max(xra(kk),xplot(kk,i))
!     write (nout,1299) kk, kkmax
 1299 format (" shapenp, kk, kkmax",2i5)
    4 continue
    3 continue
      xlb(kk) = xplot(kk,2)
      do 6 i=1,kkmax
    6 xlb(kk) = min(xlb(kk),xplot(kk,i))
!
      ztopa(kk)=-1000._R8
      do 2 i=2,kkmax
      if(ztopa(kk).gt.zplot(kk,i)) go to 2
      ztopa(kk)=zplot(kk,i)
      xtopa(kk)=xplot(kk,i)
      itop = i
!     write (nout,8111) kk,kkmax
 8111 format ("shapenp, kk, kmax",2i5)
2     continue
!
!.....use cubic interpolation to find maximum point
      if(zplot(kk,itop-1) .gt. zplot(kk,itop+1)) itop = itop - 1
      denom = xplot(kk,itop+1)-xplot(kk,itop)
      if(denom .eq.0._R8) go to 100
      denom2 = xplot(kk,itop+2) - xplot(kk,itop)
      if(denom2.eq.0._R8) go to 100
      denom3 = xplot(kk,itop+1) - xplot(kk,itop-1)
      if(denom3.eq.0._R8) go to 100
      am2 = (zplot(kk,itop+1)-zplot(kk,itop-1))/denom3
      am3 = (zplot(kk,itop+2)-zplot(kk,itop  ))/denom2
      ams = (zplot(kk,itop+1)-zplot(kk,itop  ))/denom
      bcoef = (3._R8*ams - 2._R8*am2 - am3)/denom
      ccoef = (-2._R8*ams + am2 + am3)/denom**2
      disc = (3._R8*ams - am3 - am2)**2 - am2*am3
      if(disc.lt.0._R8) go to 100
      if(ccoef .eq. 0._R8)  go to 100
      xcrit = (-3._R8*ams + 2._R8*am2 + am3 - sqrt(disc))/(3._R8*ccoef*  &  
     & denom)
      if(xcrit/denom.lt.-1.E-6_R8.or. xcrit/denom.gt.1._R8+1.E-6_R8) go   &  
     & to 100
      ycrit = am2*xcrit + bcoef*xcrit**2 + ccoef*xcrit**3
      ydenom = zplot(kk,itop+1)-zplot(kk,itop)
      ynorm = 0._R8
      if(ydenom .ne.0) ynorm = ycrit/ydenom
      ztopa(kk) = ycrit + zplot(kk,itop)
      if(ycrit .le. 0) ztopa(kk) = zplot(kk,itop)
      xtopa(kk) = xcrit + xplot(kk,itop)
!
      if(isym.eq.0) go to 200
      zbota(kk) = -ztopa(kk)
      xbota(kk) = xtopa(kk)
      go to 300
  200 continue
      zbota(kk)=+1000._R8
      do 5 i=2,kkmax
      if(zbota(kk).lt.zplot(kk,i)) go to 5
      zbota(kk)=zplot(kk,i)
      xbota(kk)=xplot(kk,i)
      ibot = i
!     write (nout,5555) kk,kkmax
 5555 format ("junk, kk, kkmax",2i5)
5     continue
!
!.....use cubic interpolation to find minimum point
      if(zplot(kk,ibot-1) .lt. zplot(kk,ibot+1)) ibot = ibot - 1
      denom = xplot(kk,ibot+1)-xplot(kk,ibot)
      if(denom .eq.0) go to 101
      denom2 = xplot(kk,ibot+2) - xplot(kk,ibot)
      if(denom2.eq.0) go to 101
      denom3 = xplot(kk,ibot+1) - xplot(kk,ibot-1)
      if(denom3.eq.0) go to 101
      am2 = (zplot(kk,ibot+1)-zplot(kk,ibot-1))/denom3
      am3 = (zplot(kk,ibot+2)-zplot(kk,ibot  ))/denom2
      ams = (zplot(kk,ibot+1)-zplot(kk,ibot  ))/denom
      bcoef = (3._R8*ams - 2._R8*am2 - am3)/denom
      ccoef = (-2._R8*ams + am2 + am3)/denom**2
      disc = (3._R8*ams - am3 - am2)**2 - am2*am3
      if(disc.lt.0._R8) go to 101
      if(ccoef .eq. 0._R8)  go to 101
      xcrit = (-3._R8*ams + 2._R8*am2 + am3 - sqrt(disc))/(3._R8*ccoef*  &  
     & denom)
      if(xcrit/denom.lt.-1.E-6_R8.or. xcrit/denom.gt.1._R8+1.E-6_R8) go   &  
     & to 101
      ycrit = am2*xcrit + bcoef*xcrit**2 + ccoef*xcrit**3
      ydenom = zplot(kk,ibot+1)-zplot(kk,ibot)
      ynorm = 0._R8
      if(ydenom .ne.0) ynorm = ycrit/ydenom
      zbota(kk) = ycrit + zplot(kk,ibot)
      if(ycrit .ge. 0) zbota(kk) = zplot(kk,ibot)
      xbota(kk) = xcrit + xplot(kk,ibot)
!
!
! 990 write (nout,9901) denom,denom2,denom3,am2,am3,
!    1ams,bcoef,ccoef,disc,xcrit,ycrit,
!    1ydenom,zbota(kk),xbota(kk)
 9901 format ("Junkntermputa",1p4e12.4)
  300 continue
  301 continue
      f1 = (psis3-psis1)/(psis1-psis2)
      if(f1.gt.1._R8) f1 = 1._R8
      if(f1.lt.0)  f1 = 0._R8
      xl = (1._R8+f1)*xla(1) - f1*xla(2)
      xlb99 = (1._R8+f1)*xlb(1) - f1*xlb(2)
      xr = (1._R8+f1)*xra(1) - f1*xra(2)
      xbot = (1._R8+f1)*xbota(1) - f1*xbota(2)
      xtop = (1._R8+f1)*xtopa(1) - f1*xtopa(2)
      zbot = (1._R8+f1)*zbota(1) - f1*zbota(2)
      ztop = (1._R8+f1)*ztopa(1) - f1*ztopa(2)
!
!.....special diagnostic printout
!     write(nout,1212) psis1,psis2,psis3,f1
!     write(nout,1213) xra(1),xra(2),xr
!     write(nout,1214) xla(1),xla(2),xl
!     write(nout,1215) xtopa(1),xtopa(2),xtop
!     write(nout,1216) ztopa(1),ztopa(2),ztop
 1212 format(" psis1,psis2,psis3,f1",1p4e12.4)
 1213 format(" xr ",1p3e12.4)
 1214 format(" xl ",1p3e12.4)
 1215 format(" xtop",1p3e12.4)
 1216 format(" ztop",1p3e12.4)
!
      elong=(ztop-zbot)/(xr-xl)
      center=(xr+xl)/2._R8
      hwidth=(xr-center)
!     write(nout,1511) hwidth
 1511 format(" hwidth =",e12.4)
      arg = (center-.5_R8*(xtop+xbot))/hwidth
      if(abs(arg) .le. 1.0_R8) triang = asin(arg)
!
      xshape(1) = .5_R8*(xr+xl)
      xshape(2) = .5_R8*(xr-xl)
      xshape(3) = triang
      xshape(4) = 0._R8
      ellmom = elong*xshape(2)
!
!
!.....indentation parameters
      ain99 = (xl  - xlb99)/(xr  - xlb99)
!.......95% flux surface shape parameters
!
!
      el95 = elong
      del95 = triang
      ain95 = ain99
      if(i95 .le. 0) go to 1399
!     write(nout,1512) i95,psi95,xsv2(i95),xsv2(i95-1)
 1512 format("i95...",i5,1p3e12.4)
      f2 = (psi95-xsv2(i95-1))/(xsv2(i95)-xsv2(i95-1))
      xl95 = (1._R8-f2)*xla(3) + f2*xla(4)
      xlb95 = (1._R8-f2)*xlb(3) + f2*xlb(4)
      xr95 = (1._R8-f2)*xra(3) + f2*xra(4)
      xbot95 = (1._R8-f2)*xbota(3) + f2*xbota(4)
      xtop95 = (1._R8-f2)*xtopa(3) + f2*xtopa(4)
      zbot95 = (1._R8-f2)*zbota(3) + f2*zbota(4)
      ztop95 = (1._R8-f2)*ztopa(3) + f2*ztopa(4)
      el95 = (ztop95-zbot95)/(xr95-xl95)
!
!.....special diagnostic printout
!     write(nout,1312) psis1,psis2,psis3,f2
!     write(nout,1313) xra(1),xra(2),xr
!     write(nout,1314) xla(1),xla(2),xl
!     write(nout,1315) xtopa(1),xtopa(2),xtop95
!     write(nout,1316) ztopa(1),ztopa(2),ztop95
 1312 format(" psis1,psis2,psis3,f2",1p4e12.4)
 1313 format(" xr ",1p3e12.4)
 1314 format(" xl ",1p3e12.4)
 1315 format(" xtop95",1p3e12.4)
 1316 format(" ztop95",1p3e12.4)
!
      cent95 = (xr95+xl95)/2._R8
      hw95 = xr95-cent95
      arg95 = (cent95-0.5_R8*(xtop95+xbot95))/hw95
      if(abs(arg95).le.1.0_R8) del95 = asin(arg95)
      ain95 = (xl95- xlb95)/(xr95- xlb95)
!
 1399 continue
!.......90% flux surface shape parameters
!
!
      el90 = elong
      del90 = triang
      ain90 = ain99
      if(i90 .le. 0) go to 1398
      f2 = (psi90-xsv2(i90-1))/(xsv2(i90)-xsv2(i90-1))
      xl90 = (1._R8-f2)*xla(5) + f2*xla(6)
      xlb90 = (1._R8-f2)*xlb(5) + f2*xlb(6)
      xr90 = (1._R8-f2)*xra(5) + f2*xra(6)
      xbot90 = (1._R8-f2)*xbota(5) + f2*xbota(6)
      xtop90 = (1._R8-f2)*xtopa(5) + f2*xtopa(6)
      zbot90 = (1._R8-f2)*zbota(5) + f2*zbota(6)
      ztop90 = (1._R8-f2)*ztopa(5) + f2*ztopa(6)
!
!.....special diagnostic printout
!     write(nout,1412) psis1,psis2,psis3,f1
!     write(nout,1413) xra(1),xra(2),xr
!     write(nout,1414) xla(1),xla(2),xl
!     write(nout,1415) xtopa(1),xtopa(2),xtop90
!     write(nout,1416) ztopa(1),ztopa(2),ztop90
 1412 format(" psis1,psis2,psis3,f1",1p4e12.4)
 1413 format(" xr ",1p3e12.4)
 1414 format(" xl ",1p3e12.4)
 1415 format(" xtop90",1p3e12.4)
 1416 format(" ztop90",1p3e12.4)
!
      el90 = (ztop90-zbot90)/(xr90-xl90)
      cent90 = (xr90+xl90)/2._R8
      hw90 = xr90-cent90
      arg90 = (cent90-0.5_R8*(xtop90+xbot90))/hw90
      if(abs(arg90).le.1.0_R8) del90 = asin(arg90)
      ain90 = (xl90- xlb90)/(xr90- xlb90)
 1398 continue
      return
  100 write(nout,1001) itop,kkmax,xplot(kk,itop-1),xplot(kk,itop),       &  
     &xplot(kk,itop+1),xplot(kk,itop+2),zplot(kk,itop-1),zplot(kk,itop),  &  
!    &                                                                   &  
     &zplot(kk,itop+1),zplot(kk,itop+2),denom,denom2,denom3,am2,am3,     &  
     &ams,bcoef,ccoef,disc,xcrit,ycrit,ynorm,ztop,xtop
 1001 format(" error in shapenp ...itop,kkmax= ",2i5,/,                  &  
     &  1p10e12.4,/,1p10e12.4,/,1p10e12.4)
      ifail = 1
      return
  101 write(nout,1011) ibot,kkmax,xplot(kk,ibot-1),xplot(kk,ibot),       &  
     &xplot(kk,ibot+1),xplot(kk,ibot+2),zplot(kk,ibot-1),zplot(kk,ibot),  &  
!    &                                                                   &  
     &zplot(kk,ibot+1),zplot(kk,ibot+2),denom,denom2,denom3,am2,am3,     &  
     &ams,bcoef,ccoef,disc,xcrit,ycrit,ynorm,zbot,xbot
 1011 format(" error in shapenp ...ibot,kkmax= ",2i5,/,                  &  
     &       1p10e12.4,/,1p10e12.4,/,1p10e12.4)
      ifail = 1
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
