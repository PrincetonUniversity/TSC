      subroutine vplot(itype)
!......6.20 vplot
!
!....NOTE:  changed sign on option 4 on 5/13/2010
!           other options likely need sign change also
!
!
!.....makes vector plots
!       ITYPE   PLOT
!         1     poloidal magnetic field
!         2     incompressible velocity field
!         3     forces
!         4     poloidal current
!         5     compressible velocity field
!         6     total velocity field
!         7     total E poloidal
!
      USE CLINAM
      USE SCRATCH
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,lref,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 speedm,xe,area,ai,aj,oi,oj,vecx1,vecz1,vecx2,vecz2,g5
      REAL*8 g6,g8,g9,etaa_local,etab_local,                             &  
     &       etac,etad,gi,gj,vecx3,vecz3,gis,gjs
      REAL*8 speed,vmax,scale,x858,xco,x142,zmin,xmax,zmax,xcosv,xx
      REAL*8 zz,vx,vz
      REAL*8 sign
!============
      lref = 0
      sign = 1._R8
      speedm = 0
!                                       Start of loop
      do 200 i=2,nxp
      do 100 j=2,nzp
      vecx(i,j) = 0._R8
      vecz(i,j) = 0._R8
      xe = xary(i)
      area = dxdz*xary(i)
      go to(10,20,30,40,50,60,70),itype
   10 continue
      ai = (psi(i+1,j)-psi(i-1,j))*.5_R8
      aj = (psi(i,j+1)-psi(i,j-1))*.5_R8
      vecx(i,j) = (aj*deex)/area
      vecz(i,j) = (-ai*deez)/area
      go to 99
   20 continue
      ai = (abig(i+1,j)-abig(i-1,j))*.5_R8/deex
      aj = (abig(i,j+1)-abig(i,j-1))*.5_R8/deez
      oi = .5_R8*(omeg(i+1,j+1)+omeg(i+1,j)-omeg(i,j+1)-omeg(i,j))/deex
      oj = .5_R8*(omeg(i+1,j+1)+omeg(i,j+1)-omeg(i+1,j)-omeg(i,j))/deez
      vecx(i,j) = -(ai/xary(i)**2 - oj/xary(i) )*gs(i,j)*udsv
      vecz(i,j) = -(aj/xary(i)**2 + oi/xary(i) )*gs(i,j)*udsv
!
      go to 99
   30 continue
      ai = .25_R8*( (psi(i+1,j)-psi(i,j))*                               &  
     &           (w(i+1,j)/ajey(i+1)+w(i+1,j+1)/ajey(i+1))               &  
     &          +(psi(i,j)-psi(i-1,j))*                                  &  
     &           (w(i  ,j)/ajey(i  )+w(i  ,j+1)/ajey(i  )) )/deex
      aj = .25_R8*( (psi(i,j+1)-psi(i,j))*                               &  
     &           (w(i,j+1)/ajey(i  )+w(i+1,j+1)/ajey(i+1))               &  
     &          +(psi(i,j)-psi(i,j-1))*                                  &  
     &           (w(i,j  )/ajey(i  )+w(i+1,j  )/ajey(i+1)) )/deez
      vecx(i,j) =  ai/xary(i)**2*udsv
      vecz(i,j) =  aj/xary(i)**2*udsv
      go to 99
   70 continue
      ai = (abig(i+1,j)-abig(i-1,j))*.5_R8/deex
      aj = (abig(i,j+1)-abig(i,j-1))*.5_R8/deez
      oi = .5_R8*(omeg(i+1,j+1)+omeg(i+1,j)-omeg(i,j+1)-omeg(i,j))/deex
      oj = .5_R8*(omeg(i+1,j+1)+omeg(i,j+1)-omeg(i+1,j)-omeg(i,j))/deez
      vecx1     = -(ai/xary(i)**2 - oj/xary(i) )*gs(i,j)*udsv
      vecz1     = -(aj/xary(i)**2 + oi/xary(i) )*gs(i,j)*udsv
      ai = .25_R8*( (psi(i+1,j)-psi(i,j))*                               &  
     &           (w(i+1,j)/ajey(i+1)+w(i+1,j+1)/ajey(i+1))               &  
     &          +(psi(i,j)-psi(i-1,j))*                                  &  
     &           (w(i  ,j)/ajey(i  )+w(i  ,j+1)/ajey(i  )) )/deex
      aj = .25_R8*( (psi(i,j+1)-psi(i,j))*                               &  
     &           (w(i,j+1)/ajey(i  )+w(i+1,j+1)/ajey(i+1))               &  
     &          +(psi(i,j)-psi(i,j-1))*                                  &  
     &           (w(i,j  )/ajey(i  )+w(i+1,j  )/ajey(i+1)) )/deez
      vecx2     =  ai/xary(i)**2*udsv
      vecz2     =  aj/xary(i)**2*udsv
      g5 = g(i,j)*xsqoj(i)
      g6 = g(i+1,j)*xsqoj(i+1)
      g8 = g(i,j+1)*xsqoj(i)
      g9 = g(i+1,j+1)*xsqoj(i+1)
      etaa_local = .5_R8*(etay(i,j)+etay(i,j-1))
      etab_local = .5_R8*(etay(i,j)+etay(i+1,j))
      etac = .5_R8*(etay(i,j)+etay(i,j+1))
      etad = .5_R8*(etay(i,j)+etay(i-1,j))
      gi = .5_R8*( (g9-g8)*etac + (g6-g5)*etaa_local )/deex
      gj = .5_R8*( (g8-g5)*etad + (g9-g6)*etab_local )/deez
      vecx3     = -gj/xary(i)*udsv
      vecz3     =  gi/xary(i)*udsv
!
      vecx(i,j) = vecx1+vecx2+vecx3
      vecz(i,j) = vecz1+vecz2+vecz3
      go to 99
   40 continue
      g5 = g(i,j)*xsqoj(i)
      g6 = g(i+1,j)*xsqoj(i+1)
      g8 = g(i,j+1)*xsqoj(i)
      g9 = g(i+1,j+1)*xsqoj(i+1)
      ai = .5_R8*(g9+g6-g8-g5)
      aj = .5_R8*(g8+g9-g5-g6)
!
!.....changed sign on 5/13/2010
      vecx(i,j) = -(aj*deex)/area
      vecz(i,j) = (ai*deez)/area
      go to 99
   50 continue
      g5 = g(i,j)*xsqoj(i)
      g6 = g(i+1,j)*xsqoj(i+1)
      g8 = g(i,j+1)*xsqoj(i)
      g9 = g(i+1,j+1)*xsqoj(i+1)
      etaa_local = .5_R8*(etay(i,j)+etay(i,j-1))
      etab_local = .5_R8*(etay(i,j)+etay(i+1,j))
      etac = .5_R8*(etay(i,j)+etay(i,j+1))
      etad = .5_R8*(etay(i,j)+etay(i-1,j))
      gi = .5_R8*( (g9-g8)*etac + (g6-g5)*etaa_local )/deex
      gj = .5_R8*( (g8-g5)*etad + (g9-g6)*etab_local )/deez
      vecx(i,j) = -gj/xary(i)*udsv
      vecz(i,j) =  gi/xary(i)*udsv
      go to 99
   60 continue
      if(psi(i,j).gt.psilim) go to 99
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 99
      gis = (deez**2)*(xe/area)**2
      gjs = (deex**2)*(xe/area)**2
      oi = .5_R8*(omeg(i+1,j+1)+omeg(i+1,j)-omeg(i,j+1)-omeg(i,j))
      oj = .5_R8*(omeg(i+1,j+1)+omeg(i,j+1)-omeg(i+1,j)-omeg(i,j))
      ai = (abig(i+1,j)-abig(i-1,j))*.5_R8
      aj = (abig(i,j+1)-abig(i,j-1))*.5_R8
      vecx(i,j) = deex*oi*gis/udst
      vecz(i,j) =  deez*oj*gjs/udst
      vecx(i,j) = vecx(i,j) +(aj*deex)/area/udst
      vecz(i,j) = vecz(i,j) +(-ai*deez)/area/udst
   99 continue
      speed = vecx(i,j)**2 + vecz(i,j)**2
      if(speed.gt.speedm) speedm = speed
  100 continue
  200 continue
!                                       End of loop
      vmax = sqrt(speedm)
      if(vmax.lt.small) go to 500
      scale = 8.0_R8* deex/vmax
      x858 = .858_R8
      if(alx-ccon .gt. 2._R8*alz*.716_R8/.55_R8) go to 1212
      xco = .45_R8
      x142 = .858_R8- .55_R8*(alx-ccon)/(2._R8*alz)
      go to 1213
 1212 continue
      x142 = .142_R8
      xco = 1._R8- .716_R8*(2._R8*alz/(alx-ccon))
 1213 continue
      zmin = -alz
      xmax = alx
      zmax = alz
      xcosv = xco
!
      call maps(ccon,xmax,zmin,zmax,x142,x858,xco,1._R8)
!
!.....draw arrows
      do 400 i=2,nxp
      do 300 j=2,nzp
      xx = xary(i)
      zz = zary(j)
      speed = vecx(i,j)**2 + vecz(i,j)**2
      if(speed .lt. (.001_R8)*speedm) go to 300
      vx = vecx(i,j)*scale
      vz = vecz(i,j)*scale
      call arrow(xx,zz,vx,vz)
  300 continue
  400 continue
!
      if (isym.ne.0)   then
      do 401 i=2,nxp
      do 301 j=2,nzp
      xx = xary(i)
      zz = -zary(j)
      speed = vecx(i,j)**2 + vecz(i,j)**2
      if(speed.lt. (.001_R8)*speedm) go to 301
      vx = vecx(i,j)*scale
      vz = vecz(i,j)*scale
      go to (303,303,302,303,303,302,303),itype
  302 vz = - vz
      go to 304
  303 vx = - vx
  304 continue
      call arrow(xx,zz,vx,vz)
  301 continue
  401 continue
                       endif
!
!
      call setld(1._R8,40._R8,1,0,1,1)
      go to(101,102,103,104,105,106,107),itype
  101 write(s100,1001)vmax,times,kcycle
      call gtextm(s100,80,0,1,2)
 1001 format(1x," poloidal field,max=",1pe12.4,/,"  time=",1pe12.4       &  
     &      ,"  cycle=",i7)
      write(nsc1,1011) kcycle
 1011 format(" poloidal field  cycle=",i7)
      go to 110
  102 write(s100,1002)vmax,times,kcycle
      call gtextm(s100,80,0,1,2)
 1002 format(1x," Vp X Bt, max=",1pe12.4,/,"  time=",1pe12.4             &  
     &      ,"  cycle=",i7)
      write(nsc1,1012) kcycle
 1012 format(" Vp X Bt  cycle=",i7)
      go to 110
  103 write(s100,1003)vmax,times,kcycle
      call gtextm(s100,80,0,1,2)
 1003 format(1x,"  vT X Bp   ,max=",1pe12.4,/,"   time=",1pe12.4,        &  
     &     "  cycle=",i7)
      write(nsc1,1013) kcycle
 1013 format(" vT X Bp   cycle=",i7)
!
      go to 110
  104 write(s100,1004)vmax,times,kcycle
      call gtextm(s100,80,0,1,2)
 1004 format(1x,"pol cur ,max=",1pe12.4,/,"  time=",1pe12.4,             &  
     &     " cycle=",i7)
      write(nsc1,1014) kcycle
 1014 format(" pol current   cycle=",i7)
      go to 110
  105 write(s100,1005)vmax,times,kcycle
      call gtextm(s100,80,0,1,2)
 1005 format(1x," eta*j-poloidal, max=",1pe12.4,/,"  time=",1pe12.4      &  
     &      ,"  cycle=",i7)
      write(nsc1,1015) kcycle
 1015 format(" eta*j-poloidal  cycle=",i7)
      go to 110
  106 write(s100,1006)vmax,times,kcycle
      call gtextm(s100,80,0,1,2)
 1006 format(1x," pol velocity(total), max=",1pe12.4,/,"  time=",        &  
     & 1pe12.4                                                           &  
     &      ,"  cycle=",i7)
      write(nsc1,1016) kcycle
 1016 format(" pol velocity(total)cycle=",i7)
      go to 110
  107 write(s100,1007)vmax,times,kcycle
      call gtextm(s100,80,0,1,2)
 1007 format(1x," total E-poloidal, max=",1pe12.4,/,"  time=",1pe12.4    &  
     &      ,"  cycle=",i7)
      write(nsc1,1017) kcycle
 1017 format(" total E-poloidal  cycle=",i7)
  110 continue
  499 call frscj(6)
  500 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
