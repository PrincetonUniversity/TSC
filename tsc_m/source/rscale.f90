      subroutine rscale
!
!......6.93 savprof
!
!.....define dimensionless quantities and vacuum temperature
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.....first check for errors
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,l,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 allam,tedge,dx,dz,xx,rsmax,amag,dpdz,dpdx,ajmid,gmid
      REAL*8 dfdz,xmidsq,dfdx,denom
!============
      if(udsd .le. 0 .or. ffac .le. 0 .or. tevv.le.0) then
      write(nout,1011) udsd,ffac,tevv
 1011 format(" * * * error * * *..udsd,ffac,tevv=",1p3e12.4)
      ineg=21
      return
      endif
!
      udst  = sqrt((4._R8*pi*1.E-7_R8)*udsd*(1.67E-27_R8))*ffac
      udsv  = 1._R8/udst
      udsr  = (4._R8*pi*1.E-7_R8)/udst
      udse  = 1._R8/udst
!
      usdt  = 1._R8/udst
      usdv  = 1._R8/udsv
      usde  = 1._R8/udse
      usdr  = 1._R8/udsr
!
      dtmin   = dtmins*usdt
      dtmax   = dtmaxs*usdt
!
!.....vacuum resistivity
      allam   = 24._R8-log(1.E-3_R8*sqrt(udsd*r0*fracn0)/tevv)
      etav    = (0.5_R8*1.03E-4_R8*allam)*tevv**(-1.5_R8)*usdr*zeff
      tedge = tevv
      if(whalos.gt.0 .or. teflat_time.gt.0) tedge = max(tevv,thalos)
      if(acoef(880) .gt. 0) tedge=acoef(880)
      smallt  = tedge*(usdh/usdd)
!
!.....halo resistivity
      etah = etav
      if(thalos.gt.0) then
        allam   = 24._R8-log(1.E-3_R8*sqrt(udsd*r0*fracn0)/thalos)
        etah=(0.5_R8*1.03E-4_R8*allam)*thalos**(-1.5_R8)*usdr*zeff
      endif
!
      dx = (alx-ccon)/(nx-1)
      dz = 2._R8*alz/((nz-1)*(1+isym))
      do 518 n=1,ncoil
  518 rscoil(n) = rscoils(n)*usdr
      do 517 n=1,nwire
      xx = xwire(n)
      rsmax = etav*tpi*xx/(dx*dz)
      if(rswire(n).lt.rsmax) go to 604
      rswire(n) = .9_R8*rsmax
      rswires(n) = .9_R8*rsmax*udsr
      write(nout,1603) n,rswires(n)
 1603 format("  wire",i3," resistivity decreased to max of",1pe12.4)
!
  604 continue
  517 rswire(n) = rswires(n)*usdr
      time = times*usdt
      do 520 l=1,ptpts
  520 tpro(l) = tpros(l)*usdt
      do 521 l=1,numfb
      tfbon(l) = tfbons(l)*usdt
      tfbof(l) = tfbofs(l)*usdt
  521 continue
!
!
      if(kcycle.le.0 .or. lrswtch.gt.0) return
      resid = 0._R8
      amag = 0._R8
      do 333 i=iminn,imaxx
      do 333 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 333
      if(psi(i,j).gt.psilim) go to 333
      dpdz = .25_R8*(psi(i,j+1)+psi(i-1,j+1)-psi(i,j-1)-psi(i-1,j-1))
      dpdx = .25_R8*(psi(i+1,j)+psi(i+1,j-1)-psi(i-1,j)-psi(i-1,j-1))
      ajmid = .5_R8*(ajphi(i,j)+ajphi(i-1,j))*xarh(i)
      gmid = .5_R8*(g(i,j+1)+g(i,j))*xsqoj(i)
      dfdz = (ajmid*dpdz                                                 &  
     &         + (gmid)*(g(i,j+1)-g(i,j))*xsqoj(i)                       &  
     &         + xarh(i)**2*(pr(i,j+1)-pr(i,j)))/deez
      ajmid = .5_R8*(ajphi(i,j)+ajphi(i,j-1))*xary(i)
      gmid = .5_R8*(g(i+1,j)*xsqoj(i+1)+g(i,j)*xsqoj(i))
      xmidsq = xary(i)**2
      dfdx = (ajmid*dpdx                                                 &  
     &         + (gmid)*(g(i+1,j)*xsqoj(i+1)-g(i,j)*xsqoj(i))            &  
     &         + xmidsq*(pr(i+1,j)-pr(i,j)))/deex
      denom = (dpdz**2 + dpdx**2) + 1.E-12_R8
      resid = resid + (dpdz*dfdz+dpdx*dfdx)**2/denom
      amag = amag + (ajmid*(dpdz/deez + dpdx/deex))**2
  333 continue
      if(amag.le.0) return
      resid = sqrt(resid/amag)
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
