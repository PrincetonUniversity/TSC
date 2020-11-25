      subroutine xptcalc
!
      USE CLINAM
      USE SVDCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER icount,icount2,ifstsvd,i,i1,i2,ip,im,j,ic,jc,j2,j4,k1
      INTEGER k3,mrows,ncol,ierr,null,k,iunique
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rco,zco,ph,relerr,rc,zc,dis,r1,r3,z2,z4,pc,pinterp,p1
      REAL*8 p2,p3,p4,tau,sigmax,sigmin,sigrat,rsq,rmsnorm,denom
      REAL*8 drsep,dzsep,rxx1,rzz1,rxz1,rxx2,rzz2,rxz2,drsep2
      REAL*8 dzsep2,c1,c2,c3,dr10,rl,zl,rr,zr,rt,zt,pval,dr20,rinb
      REAL*8 zinb,rine,zine,smax,theta,r0s,z0s,costh,sinth,fac,fac2
      REAL*8 fac3,sumd,routb,zoutb,route,zoute,y1,y2
      REAL*8 sum
!============
      dimension rco(nptsp),zco(nptsp),ph(7,4)
!
!.....coordinates of experimental flux loops
      data rco/ 1.143_R8, 1.168_R8, 1.244_R8, 1.346_R8, 1.447_R8,        &  
     & 1.549_R8,                                                         &  
     &          1.651_R8, 1.727_R8, 1.793_R8/
      data zco/ 1.250_R8, 1.346_R8, 1.432_R8, 1.519_R8, 1.574_R8,        &  
     & 1.582_R8,                                                         &  
     &          1.536_R8, 1.455_R8, 1.371_R8/
      data icount,icount2/10,0/
      data ifstsvd/1/
!
      if(xsep(1) .le. 0) return
      xsepcal = xsep(1)
      zsepcal = zsep(1)
      psepcal = psisep
      if(isvd .le. 0) return
!
!.....set up design matrix for 7 basis functions expanding around (rs,zs)
      rs = 1.46_R8
      zs = 1.27_R8
      relerr = 1.E-5_R8
!
      rsi = 1._R8/rs
      do 100 i=1,nptsp
      rcocom(i) = rco(i)
      zcocom(i) = zco(i)
      i1 = 2*(i-1) + 1
      i2 = 2*(i-1) + 2
      rc = rco(i)
      zc = zco(i)
!
!.....compute normal vector
      ip = i+1
      im = i-1
      if(i.eq.1) im = i
      if(i.eq.nptsp) ip = i
      dis = sqrt((rco(ip)-rco(im))**2+(zco(ip)-zco(im))**2)
      rnormv(i) = (zco(ip)-zco(im))/dis
      znormv(i) =-(rco(ip)-rco(im))/dis
!
!
      if(ifstsvd.eq.0) go to 201
      do 200 j=1,7
!
!.....flux basis vectors at point i
      am(i1,1) = 1._R8
      am(i1,2) = (rc-rs) + 0.5_R8*rsi*(zc-zs)**2*(1._R8-rsi*(rc-rs))
      am(i1,3) = (zc-zs)
      am(i1,4) = (rc-rs)**2 - (zc-zs)**2 + rsi*(rc-rs)*(zc-zs)**2
      am(i1,5) = (rc-rs)*(zc-zs)*(1._R8+0.5_R8*rsi*(rc-rs))
      am(i1,6) = (zc-zs)*((zc-zs)**2 - 3._R8*(rc-rs)**2)
      am(i1,7) = (rc-rs)*((rc-rs)**2 - 3._R8*(zc-zs)**2)
!
!.....normal derivative at point I
      am(i2,1) = 0._R8
      am(i2,2) = rnormv(i)*(1._R8- 0.5_R8*rsi**2*(zc-zs)**2)             &  
     &         + znormv(i)*(rsi*(zc-zs)*(1._R8-rsi*(rc-rs)))
      am(i2,3) = 0._R8                                                   &  
     &         + znormv(i)
      am(i2,4) = rnormv(i)*(2._R8*(rc-rs) + rsi*(zc-zs)**2)              &  
     &         + znormv(i)*(-2._R8*(zc-zs) + 2._R8*rsi*(rc-rs)*(zc-zs))
      am(i2,5) = rnormv(i)*((zc-zs)*(1._R8+rsi*(rc-rs)))                 &  
     &         + znormv(i)*((rc-rs)*(1._R8+0.5_R8*rsi*(rc-rs)))
      am(i2,6) = rnormv(i)*(-6._R8*(rc-rs)*(zc-zs))                      &  
     &         + znormv(i)*(3._R8*((zc-zs)**2 - (rc-rs)**2))
      am(i2,7) = rnormv(i)*(3._R8*((rc-rs)**2 - (zc-zs)**2))             &  
     &         + znormv(i)*(-6._R8*(rc-rs)*(zc-zs))
!
  200 continue
  201 continue
!
!.....flux and derivatives at point I from Data
      r1 = rc + .01_R8
      r3 = rc - .01_R8
      z2 = zc + .01_R8
      z4 = zc - .01_R8
      ic = (rc-ccon)/deex + 2
      jc = (zc-zzero)/deez + 2
      pc = pinterp(rc,zc,ic,jc)
      j2 = (z2-zzero)/deez + 2
      j4 = (z4-zzero)/deez + 2
      k1 = (r1-ccon)/deex + 2
      k3 = (r3-ccon)/deex + 2
      p1 = pinterp(r1,zc,k1,jc)
      p2 = pinterp(rc,z2,ic,j2)
      p3 = pinterp(r3,zc,k3,jc)
      p4 = pinterp(rc,z4,ic,j4)
      datav(i1) = pc
      datav(i2) = rnormv(i)*(p1-p3)/(2._R8*.01_R8)                       &  
     &          + znormv(i)*(p2-p4)/(2._R8*.01_R8)
  100 continue
      mrows = 2*nptsp
      ncol = 7
!.....singular value decomposition
      call svd1(mrows,ncol,1,1,ierr,relerr,tau,sigmax,sigmin,sigrat,     &  
     &         rsq,rmsnorm,null,ifstsvd)
      ifstsvd = 0
      if(ierr.eq.0) go to 500
      write(nout,1001)
 1001 format(" error return from svd1")
      ineg=28
      return
  500 continue
!
!.....calculate new position of x-point
      denom = -4._R8*coef(4)**2 + 2._R8*rsi*coef(4)*coef(2)              &  
     &      - coef(5)**2
      if(denom.eq.0) go to 700
      drsep = (2._R8*coef(4)*coef(2)-rsi*coef(2)**2+coef(5)*coef(3))/    &  
     & denom
      dzsep = (coef(5)*coef(2)-2._R8*coef(4)*coef(3))/denom
!
      rxx1 = (-3._R8*coef(7)*(-2._R8*coef(4)+coef(2)*rsi)                &  
     &     - coef(5)*(3._R8*coef(6)-.5_R8*coef(5)*rsi))/denom
      rzz1 = ((-2._R8*coef(4)+coef(2)*rsi)*(.5_R8*coef(2)*rsi**2         &  
     &     - coef(4)*rsi+3._R8*coef(7))+3._R8*coef(5)*coef(6))/denom
      rxz1 =((-2._R8*coef(4)+coef(2)*rsi)*(-coef(5)*rsi+6._R8*coef(6))   &  
     &     - coef(5)*(coef(2)*rsi**2-2._R8*coef(4)*rsi+6._R8*coef(7)))/  &  
     & denom
      rxx2 = (3._R8*coef(5)*coef(7)+2._R8*coef(4)*(3._R8*coef(6)-.5_R8*  &  
     & coef(5)*rsi))                                                     &  
     &     /denom
      rzz2 = (-coef(5)*(.5_R8*coef(2)*rsi**2-coef(4)*rsi+3._R8*coef(7))  &  
     &     - 6._R8*coef(4)*coef(6))/denom
      rxz2 = (-coef(5)*(-coef(5)*rsi+6._R8*coef(6))                      &  
     &    + 2._R8*coef(4)*(coef(2)*rsi**2-2._R8*coef(4)*rsi+6._R8*       &  
     & coef(7)))/denom
      drsep2 = rxx1*drsep**2 + rzz1*dzsep**2 + rxz1*dzsep*drsep
      dzsep2 = rxx2*drsep**2 + rzz2*dzsep**2 + rxz2*dzsep*drsep
!
      xsepcal = rs + drsep + drsep2
      zsepcal = zs + dzsep + dzsep2
!.....calculate flux at xsepcal,zsepcal
      rc = xsepcal
      zc = zsepcal
      psepcal =                                                          &  
     & + coef(1)                                                         &  
     & + coef(2)*((rc-rs) + 0.5_R8*rsi*(zc-zs)**2*(1._R8-rsi*(rc-rs)))   &  
     & + coef(3)*((zc-zs))                                               &  
     & + coef(4)*((rc-rs)**2 - (zc-zs)**2 + rsi*(rc-rs)*(zc-zs)**2)      &  
     & + coef(5)*((rc-rs)*(zc-zs)*(1._R8+0.5_R8*rsi*(rc-rs)))            &  
     & + coef(6)*((zc-zs)*((zc-zs)**2 - 3._R8*(rc-rs)**2))               &  
     & + coef(7)*((rc-rs)*((rc-rs)**2 - 3._R8*(zc-zs)**2))
!
!
!.....calculate error parameters here
!
!.....error eps1...outer midplane..expand about r=2.35,z=0
      eps10 = -.05_R8
      rc = 2.35_R8
      zc = 0.0_R8
      r1 = rc + .01_R8
      ic = (rc-ccon)/deex + 2
      jc = (zc-zzero)/deez + 2
      c1 = pinterp(rc,zc,ic,jc)
      p1 = pinterp(r1,zc,ic,jc)
      c2 = (p1-c1)/(.01_R8)
      z2 = .05_R8
      j2 = (z2-zzero)/deez + 2
      p2 = pinterp(rc,z2,ic,j2)
      c3 = ((c1-p2) + c2*(.05_R8)**2*.5_R8/rc)/(.05_R8)**2
      dr10 = (psepcal-c1)/c2
      eps1c = dr10 - (c3/c2)*dr10**2
!
!.....binary search to find actual interface
      rl = xmag
      zl = 0._R8
      rr = rc
      zr = 0._R8
      rt = rr
      zt = zr
      pval = psisep
      call binary(rl,zl,rr,zr,pval,rt,zt)
      eps1a = rt-rc
!
      eps1a = min(eps1a,0.05_R8)
      eps1c = min(eps1c,0.05_R8)
      eps10 = min(eps10,0.05_R8)
      eps1a = max(eps1a,-.15_R8)
      eps1c = max(eps1c,-.15_R8)
      eps10 = max(eps10,-.15_R8)
!
!
!.....error eps2...inner midplane..expand about r=1.15,z=0
      eps20 = +.05_R8
      rc = 1.15_R8
      zc = 0.0_R8
      r1 = rc + .01_R8
      ic = (rc-ccon)/deex + 2
      jc = (zc-zzero)/deez + 2
      c1 = pinterp(rc,zc,ic,jc)
      p1 = pinterp(r1,zc,ic,jc)
      c2 = (p1-c1)/(.01_R8)
      z2 = .05_R8
      j2 = (z2-zzero)/deez + 2
      p2 = pinterp(rc,z2,ic,j2)
      c3 = ((c1-p2) + c2*(.05_R8)**2*.5_R8/rc)/(.05_R8)**2
      dr20 = (psepcal-c1)/c2
      eps2c = dr20 - (c3/c2)*dr20**2
!
!
!
!.....binary search to find actual interface
      rl = xmag
      zl = 0._R8
      rr = rc
      zr = 0._R8
      rt = rr
      zt = zr
      pval = psisep
      call binary(rl,zl,rr,zr,pval,rt,zt)
      eps2a = rt-rc
!
      eps2a = min(eps2a,0.15_R8)
      eps2c = min(eps2c,0.15_R8)
      eps20 = min(eps20,0.15_R8)
      eps2a = max(eps2a,-.05_R8)
      eps2c = max(eps2c,-.05_R8)
      eps20 = max(eps20,-.05_R8)
!
!.....error eps3 ... inner strike-point
      rinb = 1.280_R8
      zinb = 1.492_R8
      rine = 1.202_R8
      zine = 1.306_R8
      smax = sqrt((rinb-rine)**2 + (zinb-zine)**2)
      theta = atan2(zinb-zine , rinb-rine)
      r0s = rine-rs
      z0s = zine-zs
      costh = cos(theta)
      sinth = sin(theta)
      fac = 1._R8-r0s/rs
      fac2 = 1._R8+ 0.5_R8*r0s/rs
      fac3 = 1._R8+     r0s/rs
!
      ph(1,1) = 1._R8
      ph(1,2) = 0._R8
      ph(1,3) = 0._R8
      ph(1,4) = 0._R8
!
      ph(2,1) = r0s + .5_R8/rs*z0s**2*fac
      ph(2,2) = costh*(1._R8-0.5_R8/rs**2*z0s**2) + z0s/rs*sinth*fac
      ph(2,3) = sinth*(-z0s*costh/rs**2 + 0.5_R8*sinth*fac/rs)
      ph(2,4) = -sinth**2*costh*0.5_R8/rs**2
!
      ph(3,1) = z0s
      ph(3,2) = sinth
      ph(3,3) = 0._R8
      ph(3,4) = 0._R8
!
      ph(4,1) = r0s**2 - z0s**2 + r0s*z0s**2/rs
      ph(4,2) =2._R8*r0s*costh-2._R8*z0s*sinth+(2._R8*r0s*z0s*sinth+     &  
     & z0s**2*costh)                                                     &  
     &        /rs
      ph(4,3) = costh**2-sinth**2+(r0s*sinth**2 + 2._R8*costh*sinth*z0s)  &  
     & /rs
      ph(4,4) = costh*sinth**2/rs
!
      ph(5,1) = r0s*z0s*fac2
      ph(5,2) = costh*z0s*fac3 + sinth*r0s*fac2
      ph(5,3) = costh*(sinth*fac3 + costh*z0s*0.5_R8/rs)
      ph(5,4) = sinth*costh**2*0.5_R8/rs
!
      ph(6,1) = z0s*(z0s**2 - 3._R8*r0s**2)
      ph(6,2) = z0s*(3._R8*z0s*sinth - 6._R8*r0s*costh) - 3._R8*r0s**2*  &  
     & sinth
      ph(6,3) = -3._R8*z0s*costh**2 + 3._R8*z0s*sinth**2 - 6._R8*r0s*    &  
     & sinth*costh
      ph(6,4) = sinth*(sinth**2 - 3._R8*costh**2)
!
      ph(7,1) = r0s*(r0s**2 - 3._R8*z0s**2)
      ph(7,2) = r0s*(3._R8*r0s*costh - 6._R8*z0s*sinth) - 3._R8*z0s**2*  &  
     & costh
      ph(7,3) = 3._R8*r0s*(costh**2 - sinth**2) - 6._R8*sinth*costh*z0s
      ph(7,4) = costh*(costh**2 - 3._R8*sinth**2)
!
!.....apply newtons method about old time value
      sum = 0._R8
      sumd = 0._R8
      if(abs(eps3co) .le. 1.E-6_R8) eps3co = 1.E-6_R8
      if(eps3co .gt. 1.5_R8*smax) eps3co = 1.5_R8*smax
      if(eps3co .lt.-0.5_R8*smax) eps3co =-0.5_R8*smax
      do 803 i=1,7
      do 802 k=1,4
      if(k.le.1) go to 802
      sumd = sumd + coef(i)*ph(i,k)*(k-1)*eps3co**(k-2)
  802 sum  = sum  + coef(i)*ph(i,k)*eps3co**(k-1)
  803 continue
      eps3c = eps3co + (psepcal-sum)/sumd
      eps3co = eps3c
      rt = rinb
      zt = zinb
      call binary(rinb,zinb,rine,zine,psisep,rt,zt)
      eps3a = sqrt((rt-rine)**2 + (zt-zine)**2)
      if(rt.lt.rine) eps3a = -eps3a
      if(times .le. 5.8_R8) eps30 = .5_R8*smax
      if(times .ge.10.8_R8) eps30 = 0._R8
      if(times .gt. 5.8_R8.and. times .lt. 10.8_R8)                      &  
     & eps30 = .5_R8*smax*(1._R8-(times-5.8_R8)/5._R8)
!
      eps3c = min(1.5_R8*smax,eps3c)
      eps30 = min(1.5_R8*smax,eps30)
      eps3a = min(1.5_R8*smax,eps3a)
      eps3c = max(-0.5_R8*smax,eps3c)
      eps30 = max(-0.5_R8*smax,eps30)
      eps3a = max(-0.5_R8*smax,eps3a)
!
!
!.....error eps4 ... outer strike-point
!     routb = 1.593
!     zoutb = 1.485
      routb = 1.782_R8
      zoutb = 1.419_R8
      route = 1.404_R8
      zoute = 1.551_R8
      smax = sqrt((routb-route)**2 + (zoutb-zoute)**2)
      theta = atan2(zoutb-zoute , routb-route)
      r0s = route-rs
      z0s = zoute-zs
      costh = cos(theta)
      sinth = sin(theta)
      fac = 1._R8-r0s/rs
      fac2 = 1._R8+ 0.5_R8*r0s/rs
      fac3 = 1._R8+     r0s/rs
!
      ph(1,1) = 1._R8
      ph(1,2) = 0._R8
      ph(1,3) = 0._R8
      ph(1,4) = 0._R8
!
      ph(2,1) = r0s + .5_R8/rs*z0s**2*fac
      ph(2,2) = costh*(1._R8-0.5_R8/rs**2*z0s**2) + z0s/rs*sinth*fac
      ph(2,3) = sinth*(-z0s*costh/rs**2 + 0.5_R8*sinth*fac/rs)
      ph(2,4) = -sinth**2*costh*0.5_R8/rs**2
!
      ph(3,1) = z0s
      ph(3,2) = sinth
      ph(3,3) = 0._R8
      ph(3,4) = 0._R8
!
      ph(4,1) = r0s**2 - z0s**2 + r0s*z0s**2/rs
      ph(4,2) =2._R8*r0s*costh-2._R8*z0s*sinth+(2._R8*r0s*z0s*sinth+     &  
     & z0s**2*costh)                                                     &  
     &        /rs
      ph(4,3) = costh**2-sinth**2+(r0s*sinth**2 + 2._R8*costh*sinth*z0s)  &  
     & /rs
      ph(4,4) = costh*sinth**2/rs
!
      ph(5,1) = r0s*z0s*fac2
      ph(5,2) = costh*z0s*fac3 + sinth*r0s*fac2
      ph(5,3) = costh*(sinth*fac3 + costh*z0s*0.5_R8/rs)
      ph(5,4) = sinth*costh**2*0.5_R8/rs
!
      ph(6,1) = z0s*(z0s**2 - 3._R8*r0s**2)
      ph(6,2) = z0s*(3._R8*z0s*sinth - 6._R8*r0s*costh) - 3._R8*r0s**2*  &  
     & sinth
      ph(6,3) = -3._R8*z0s*costh**2 + 3._R8*z0s*sinth**2 - 6._R8*r0s*    &  
     & sinth*costh
      ph(6,4) = sinth*(sinth**2 - 3._R8*costh**2)
!
      ph(7,1) = r0s*(r0s**2 - 3._R8*z0s**2)
      ph(7,2) = r0s*(3._R8*r0s*costh - 6._R8*z0s*sinth) - 3._R8*z0s**2*  &  
     & costh
      ph(7,3) = 3._R8*r0s*(costh**2 - sinth**2) - 6._R8*sinth*costh*z0s
      ph(7,4) = costh*(costh**2 - 3._R8*sinth**2)
!
!.....apply newtons method about old time value
      sum = 0._R8
      sumd = 0._R8
      if(abs(eps4co) .le. 1.E-6_R8) eps4co = 1.E-6_R8
      if(eps4co .gt. 1.5_R8*smax) eps4co = 1.5_R8*smax
      if(eps4co .lt.-0.5_R8*smax) eps4co =-0.5_R8*smax
      do 903 i=1,7
      do 902 k=1,4
      if(k.le.1) go to 902
      sumd = sumd + coef(i)*ph(i,k)*(k-1)*eps4co**(k-2)
  902 sum  = sum  + coef(i)*ph(i,k)*eps4co**(k-1)
  903 continue
      eps4c = eps4co + (psepcal-sum)/sumd
      eps4co = eps4c
      rt = routb
      zt = zoutb
      call binary(routb,zoutb,route,zoute,psisep,rt,zt)
      eps4a = sqrt((rt-route)**2 + (zt-zoute)**2)
      if(rt.lt.route) eps4a = -eps4a
      if(times .le. 5.8_R8) eps40 = .5_R8*smax
      if(times .ge.10.8_R8) eps40 = 0._R8
      if(times .gt. 5.8_R8.and. times .lt. 10.8_R8)                      &  
     & eps40 = .5_R8*smax*(1._R8-(times-5.8_R8)/5._R8)
!
      eps4c = min(1.5_R8*smax,eps4c)
      eps40 = min(1.5_R8*smax,eps40)
      eps4a = min(1.5_R8*smax,eps4a)
      eps4c = max(-0.5_R8*smax,eps4c)
      eps40 = max(-0.5_R8*smax,eps40)
      eps4a = max(-0.5_R8*smax,eps4a)
!
!
!
      if(kcycle.le.0 .or. iprnt.le.nskipr) return
      icount2 = icount2 + 1
      if(icount2 .gt. 25) return
!
      write(nout,1200)
 1200 format(1h1,"    * * * *  results from SVD analysis * * * * ")
      write(nout,1201) sigmin,sigmax,tau,relerr,sigrat,rsq,rmsnorm,null
 1201 format(1x,//," minimum eigenvalue ",1pe12.4,                       &  
     &           /," maximum eigenvalue ",1pe12.4,                       &  
     &           /," cuttoff eigenvalue ",1pe12.4,                       &  
     &           /," relative error tol ",1pe12.4,                       &  
     &           /," ratio: min to max  ",1pe12.4,                       &  
     &           /," rms normed residual",1pe12.4,                       &  
     &           /," rms norm of data   ",1pe12.4,                       &  
     &           /," count of sv .lt.tau",i3 )
!
!......evaluate model at data points and compare with data
      write(nout,1600)
 1600 format(/////,"  comparison of fit with data ",                     &  
     &          //,10x,"  i         fit        data",                    &  
     &   "            i         fit        data")
      do 600 i=1,nptsp
      i1 = 2*(i-1) + 1
      i2 = 2*(i-1) + 2
      y1 = 0._R8
      y2 = 0._R8
      do 610 j=1,nordp
      y1 = y1 + am(i1,j)*coef(j)
      y2 = y2 + am(i2,j)*coef(j)
  610 continue
      write(nout,1610) i1,y1,datav(i1),i2,y2,datav(i2)
 1610 format(1x,2(10x,i3,1p2e12.4) )
  600 continue
!
      write(nout,1300)
 1300 format(//,"  singular values ")
      iunique = 0
      do 300 j=1,ncol
      write(nout,1103) j,sigma(j)
 1103 format(" eigenvalue no ",i3," is ",1pe12.4)
      if(sigma(j).gt.tau) go to 300
      iunique = 1
      write(nout,1101)
      write(nout,1111) (vv(i,j),i=1,ncol)
  300 continue
      if(iunique.eq.0) write(nout,1102)
 1101 format(" null coefficients follow:")
 1111 format( 1x,1p10e12.4)
 1102 format(" coefficients are unique" )
!
      write(nout,1500) psisep,psepcal,                                   &  
     &              xsep(1),xsepcal,drsep,drsep2,                        &  
     &              zsep(1),zsepcal,dzsep,dzsep2
 1500 format(///,"  psisep,  psepcal         ",1p2e12.4,                 &  
     &         /,"  xsep(1), xsepcal, drsep  ",1p4e12.4,                 &  
     &         /,"  zsep(1), zsepcal, dzsep  ",1p4e12.4 )
!
      return
!
!.....error exit
  700 continue
      write(nout,1700)
 1700 format(" denom zero in x-point calculation ")
      ineg=28
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
