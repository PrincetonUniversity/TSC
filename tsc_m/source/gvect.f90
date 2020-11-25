      subroutine gvect(r,z,xi,zi,n,g,gr,gz,grz,gzz,grr,nmult,ineg)
!.....4.90 gvect
!
      USE PARAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!......calculates derivatives wrt first argument
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,nmult,ineg,i,imult,imultp
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 z,xi,zi,g,gr,gz,grz,gzz,grr,r,a0,a1,a2,a3,a4,b0,b1,b2
      REAL*8 b3,b4,c1,c2,c3,c4,d1,d2,d3,d4,pi,tpi,rpxi,rxi,zmzi
      REAL*8 rksq,rk,sqrxi,x,ce,ck,term1,term2,rz,co
!============
      dimension r(n),z(n),xi(n),zi(n),g(n),gr(n),gz(n),                  &  
     &          grz(n),gzz(n),grr(n)
!
      data a0,a1,a2,a3,a4/1.38629436112_R8,9.666344259E-2_R8,            &  
     &3.590092383E-2_R8,3.742563713E-2_R8,1.451196212E-2_R8/
      data b0,b1,b2,b3,b4/.5_R8,.12498593597_R8,6.880248576E-2_R8,       &  
     &3.328355346E-2_R8,4.41787012E-3_R8/
      data c1,c2,c3,c4/.44325141463_R8,6.260601220E-2_R8,                &  
     &4.757383546E-2_R8,1.736506451E-2_R8/
      data d1,d2,d3,d4/.24998368310_R8,9.200180037E-2_R8,                &  
     &4.069697526E-2_R8,5.26449639E-3_R8/
      data pi,tpi/3.1415926535_R8,6.283185308_R8/
!============      
!
      do 100 i=1,n
      rpxi=r(i)+xi(i)
      rxi=r(i)*xi(i)
      zmzi=z(i)-zi(i)
      rksq=4._R8*rxi/(rpxi**2+zmzi**2)
      rk=sqrt(rksq)
      sqrxi=sqrt(rxi)
      x=1._R8-rksq
      ce=1._R8+x*(c1+x*(c2+x*(c3+x*c4)))+                                &  
     &x*(d1+x*(d2+x*(d3+x*d4)))*(-log(x))
      ck=a0+x*(a1+x*(a2+x*(a3+x*a4)))+                                   &  
     &(b0+x*(b1+x*(b2+x*(b3+x*b4))))*(-log(x))
 
      term1=2._R8*ck-2._R8*ce-ce*rksq/x
      term2=2._R8*xi(i)-rksq*rpxi
 
      g(i) =- sqrxi*(2._R8*ck-2._R8*ce-ck*rksq)/rk
      gr(i)=-rk*0.25_R8/sqrxi*(rpxi*term1                                &  
     &+2._R8*xi(i)*(ce/x-ck))
      gz(i)=-rk*0.25_R8*zmzi/sqrxi*term1
      grz(i)=0.0625_R8*zmzi*(rk/sqrxi)**3*(rpxi*term1                    &  
     &+(ce-ck+2._R8*ce*rksq/x)*                                          &  
     &(term2)/x)
      gzz(i)=-rk*0.25_R8/sqrxi*(term1*                                   &  
     &(1._R8-rksq*zmzi**2/(4._R8*rxi))+zmzi**2*rksq**2/(4._R8*rxi*x)     &  
     &*(ce-ck+2._R8*ce*rksq/x))
      grr(i)=-rk*0.25_R8/sqrxi*(-rksq*rpxi/(4._R8*rxi)*                  &  
     &(rpxi*term1+2._R8*xi(i)*(ce/x-ck))+term1-                          &  
     &rksq*rpxi/(4._R8*rxi*x)*(ce-ck+2._R8*ce*rksq/x)*                   &  
     &(term2)+rksq/(2._R8*r(i)*x)*(2._R8*ce/x-ck)*term2)
  100 continue
!
!.....check for multipolar coils
      if(nmult .le. 0) return
      do 200 i=1,n
      if(xi(i) .lt. 100._R8) go to 200
      rz = zi(i)
      imult = int(xi(i) - 100._R8)
      if(imult .lt. 0 .or. imult.gt.10) go to 250
      imultp = imult + 1
      go to(10,11,12,13,14,15,16,17,18,19,20),imultp
   10 continue
!
!....even nullapole
      g(i) = tpi*rz**2
      gr(i) = 0._R8
      gz(i) = 0._R8
      grz(i) = 0._R8
      gzz(i) = 0._R8
      grr(i) = 0._R8
      go to 200
   11 continue
!
!....odd nullapole
      g(i) = 0._R8
      gr(i) = 0._R8
      gz(i) = 0._R8
      grz(i) = 0._R8
      gzz(i) = 0._R8
      grr(i) = 0._R8
      go to 200
   12 continue
!
!....even dipole
      g(i) = tpi*(r(i)**2 - rz**2)/2._R8
      gr(i) = tpi*r(i)
      gz(i) = 0._R8
      grz(i) = 0._R8
      gzz(i) = 0._R8
      grr(i) = tpi
      go to 200
   13 continue
!
!....odd dipole
      co=tpi/rz
      g(i) = co*(r(i)**2*z(i))
      gr(i) = co*(2._R8*r(i)*z(i))
      gz(i) = co*(r(i)**2)
      grz(i) = co*2._R8*r(i)
      gzz(i) = 0._R8
      grr(i) = co*2*z(i)
      go to 200
   14 continue
!
!....even quadrapole
      co=pi/(4._R8*rz**2)
      g(i) = co*(r(i)**4-4._R8*r(i)**2*z(i)**2 - 2._R8*r(i)**2*rz**2+    &  
     & rz**4)
      gr(i) = co*(4._R8*r(i)**3-8._R8*r(i)*z(i)**2-4._R8*r(i)*rz**2)
      gz(i) = co*(-8._R8*r(i)**2*z(i))
      grz(i) = co*(-16._R8*r(i)*z(i))
      gzz(i) = co*(-8._R8*r(i)**2)
      grr(i) = co*(12._R8*r(i)**2 - 8._R8*z(i)**2 - 4._R8*rz**2)
      go to 200
   15 continue
!
!....odd quadrapole
      co=pi/(3._R8*rz**3)
      g(i) = co*r(i)**2*z(i)*(3._R8*r(i)**2-4._R8*z(i)**2-3._R8*rz**2)
      gr(i) = co*(12._R8*r(i)**3*z(i)-8._R8*r(i)*z(i)**3-6._R8*r(i)*z(i)  &  
     & *rz**2)
      gz(i) = co*(3._R8*r(i)**4 - 12._R8*r(i)**2*z(i)**2 - 3._R8*r(i)**  &  
     & 2*rz**2)
      grz(i) = co*(12._R8*r(i)**3-24._R8*r(i)*z(i)**2 - 6._R8*r(i)*rz**  &  
     & 2)
      gzz(i) = co*(-24._R8*r(i)**2*z(i))
      grr(i) = co*(36._R8*r(i)**2*z(i)-8._R8*z(i)**3-6._R8*z(i)*rz**2)
      go to 200
   16 continue
!
!.....even hexapole
      co=pi/(12._R8*rz**4)
      g(i) = co*(r(i)**6 - 12._R8*r(i)**4*z(i)**2 - 3._R8*r(i)**4*rz**2  &  
     &       + 8._R8*r(i)**2*z(i)**4 + 12._R8*r(i)**2*z(i)**2*rz**2      &  
     &       + 3._R8*r(i)**2*rz**4 - rz**6 )
      gr(i)= co*(6._R8*r(i)**5 - 48._R8*r(i)**3*z(i)**2                  &  
     &       - 12._R8*r(i)**3*rz**2 + 16._R8*r(i)*z(i)**4                &  
     &       + 24._R8*r(i)*z(i)**2*rz**2 + 6._R8*r(i)*rz**4 )
      gz(i)= co*(-24._R8*r(i)**4*z(i) + 32._R8*r(i)**2*z(i)**3           &  
     &       + 24._R8*r(i)**2*z(i)*rz**2)
      grz(i)=co*(-96._R8*r(i)**3*z(i)+64._R8*r(i)*z(i)**3                &  
     &       + 48._R8*r(i)*z(i)*rz**2 )
      grr(i)=co*(30._R8*r(i)**4-144._R8*r(i)**2*z(i)**2-36._R8*r(i)**2*  &  
     & rz**2                                                             &  
     &       + 16._R8*z(i)**4 + 24._R8*z(i)**2*rz**2 + 6._R8*rz**4)
      gzz(i)=co*(-24._R8*r(i)**4 + 96._R8*r(i)**2*z(i)**2 + 24._R8*r(i)  &  
     & **2*rz**2)
      go to 200
   17 continue
!
!.....odd hexapole
      co=pi/(30._R8*rz**5)
      g(i) = co*(15._R8*r(i)**6*z(i) - 60._R8*r(i)**4*z(i)**3            &  
     &       - 30._R8*r(i)**4*z(i)*rz**2 + 24._R8*r(i)**2*z(i)**5        &  
     &       + 40._R8*r(i)**2*z(i)**3*rz**2 + 15._R8*r(i)**2*z(i)*rz**4)     
       
      gr(i)= co*(90._R8*r(i)**5*z(i) - 240._R8*r(i)**3*z(i)**3           &  
     &       - 120._R8*r(i)**3*z(i)*rz**2 + 48._R8*r(i)*z(i)**5          &  
     &       + 80._R8*r(i)*z(i)**3*rz**2 + 30._R8*r(i)*z(i)*rz**4)
      gz(i)= co*(15._R8*r(i)**6 - 180._R8*r(i)**4*z(i)**2                &  
     &       - 30._R8*r(i)**4*rz**2 + 120._R8*r(i)**2*z(i)**4            &  
     &       +120._R8*r(i)**2*z(i)**2*rz**2 + 15._R8*r(i)**2*rz**4)
      grz(i)=co*(90._R8*r(i)**5 - 720._R8*r(i)**3*z(i)**2                &  
     &       - 120._R8*r(i)**3*rz**2 + 240._R8*r(i)*z(i)**4              &  
     &       + 240._R8*r(i)*z(i)**2*rz**2 + 30._R8*r(i)*rz**4)
      gzz(i)=co*(-360._R8*r(i)**4*z(i) + 480._R8*r(i)**2*z(i)**3         &  
     &       + 240._R8*r(i)**2*z(i)*rz**2)
      grr(i)=co*(450._R8*r(i)**4*z(i) - 720._R8*r(i)**2*z(i)**3          &  
     &       - 360._R8*r(i)**2*z(i)*rz**2 + 48._R8*z(i)**5               &  
     &       + 80._R8*z(i)**3*rz**2 + 30._R8*z(i)*rz**4)
      go to 200
   18 continue
!
!....even octapole
      co=pi/(160._R8*rz**6)
      g(i) = co*(5._R8*r(i)**8 - 120._R8*r(i)**6*z(i)**2 - 20._R8*r(i)**  &  
     & 6*rz**2                                                           &  
     &       + 240._R8*r(i)**4*z(i)**4 + 240._R8*r(i)**4*z(i)**2*rz**2   &  
     &       + 30._R8*r(i)**4*rz**4 - 64._R8*r(i)**2*z(i)**6             &  
     &       - 160._R8*r(i)**2*z(i)**4*rz**2 - 120._R8*r(i)**2*z(i)**2*  &  
     & rz**4                                                             &  
     &       - 20._R8*r(i)**2*rz**6 + 5._R8*rz**8)
      gr(i)= co*(40._R8*r(i)**7 - 720._R8*r(i)**5*z(i)**2 - 120._R8*r(i)  &  
     & **5*rz**2                                                         &  
     &       + 960._R8*r(i)**3*z(i)**4 + 960._R8*r(i)**3*z(i)**2*rz**2   &  
     &       + 120._R8*r(i)**3*rz**4 - 128._R8*r(i)*z(i)**6              &  
     &       - 320._R8*r(i)*z(i)**4*rz**2 - 240._R8*r(i)*z(i)**2*rz**4   &  
     &       - 40._R8*r(i)*rz**6)
      gz(i)= co*(-240._R8*r(i)**6*z(i) + 960._R8*r(i)**4*z(i)**3         &  
     &       + 480._R8*r(i)**4*z(i)*rz**2 - 384._R8*r(i)**2*z(i)**5      &  
     &       - 640._R8*r(i)**2*z(i)**3*rz**2 - 240._R8*r(i)**2*z(i)*rz**  &  
     & 4)
      grz(i)=co*(-1440._R8*r(i)**5*z(i) + 3840._R8*r(i)**3*z(i)**3       &  
     &       + 1920._R8*r(i)**3*z(i)*rz**2 - 768._R8*r(i)*z(i)**5        &  
     &       - 1280._R8*r(i)*z(i)**3*rz**2 - 480._R8*r(i)*z(i)*rz**4)
      gzz(i)=co*(-240._R8*r(i)**6 + 2880._R8*r(i)**4*z(i)**2             &  
     &       + 480._R8*r(i)**4*rz**2 - 1920._R8*r(i)**2*z(i)**4          &  
     &       - 1920._R8*r(i)**2*z(i)**2*rz**2 - 240._R8*r(i)**2*rz**4)
      grr(i)=co*(280._R8*r(i)**6 - 3600._R8*r(i)**4*z(i)**2              &  
     &       - 600._R8*r(i)**4*rz**2 + 2880._R8*r(i)**2*z(i)**4          &  
     &       + 2880._R8*r(i)**2*z(i)**2*rz**2                            &  
     &       + 360._R8*r(i)**2*rz**4 - 128._R8*z(i)**6                   &  
     &       - 320._R8*z(i)**4*rz**2 - 240._R8*z(i)**2*rz**4 - 40._R8*   &  
     & rz**6)
      go to 200
   19 continue
!
!....odd octapole
      co=pi/(140._R8*rz**7)
      g(i) = co*r(i)**2*z(i)*(35._R8*r(i)**6 - 280._R8*r(i)**4*z(i)**2   &  
     &       - 105._R8*r(i)**4*rz**2 + 336._R8*r(i)**2*z(i)**4           &  
     &       + 420._R8*r(i)**2*z(i)**2*rz**2 + 105._R8*r(i)**2*rz**4     &  
     &       - 64._R8*z(i)**6 - 168._R8*z(i)**4*rz**2                    &  
     &       - 140._R8*z(i)**2*rz**4 - 35._R8*rz**6)
      gr(i)= co*(280._R8*r(i)**7*z(i) - 1680._R8*r(i)**5*z(i)**3         &  
     &       - 630._R8*r(i)**5*z(i)*rz**2 + 1344._R8*r(i)**3*z(i)**5     &  
     &       + 1680._R8*r(i)**3*z(i)**3*rz**2 + 420._R8*r(i)**3*z(i)*    &  
     & rz**4                                                             &  
     &       - 128._R8*r(i)*z(i)**7 - 336._R8*r(i)*z(i)**5*rz**2         &  
     &       - 280._R8*r(i)*z(i)**3*rz**4 - 70._R8*r(i)*z(i)*rz**6)
      gz(i)= co*(35._R8*r(i)**8-840._R8*r(i)**6*z(i)**2-105._R8*r(i)**6*  &  
     & rz**2                                                             &  
     &       + 1680._R8*r(i)**4*z(i)**4 + 1260._R8*r(i)**4*z(i)**2*rz**  &  
     & 2                                                                 &  
     &       + 105._R8*r(i)**4*rz**4 - 448._R8*r(i)**2*z(i)**6           &  
     &       - 840._R8*r(i)**2*z(i)**4*rz**2 - 420._R8*r(i)**2*z(i)**2*  &  
     & rz**4                                                             &  
     &       - 35._R8*r(i)**2*rz**6)
      grz(i)=co*(280._R8*r(i)**7 - 5040._R8*r(i)**5*z(i)**2              &  
     &       - 630._R8*r(i)**5*rz**2 + 6720._R8*r(i)**3*z(i)**4          &  
     &       + 5040._R8*r(i)**3*z(i)**2*rz**2 + 420._R8*r(i)**3*rz**4    &  
     &       - 896._R8*r(i)*z(i)**6 - 1680._R8*r(i)*z(i)**4*rz**2        &  
     &       - 840._R8*r(i)*z(i)**2*rz**4 - 70._R8*r(i)*rz**6)
      gzz(i)=co*(-1680._R8*r(i)**6*z(i) + 6720._R8*r(i)**4*z(i)**3       &  
     &       + 2520._R8*r(i)**4*z(i)*rz**2 - 2688._R8*r(i)**2*z(i)**5    &  
     &       - 3360._R8*r(i)**2*z(i)**3*rz**2 - 840._R8*r(i)**2*z(i)*    &  
     & rz**4)
      grr(i)=co*(1960._R8*r(i)**6*z(i) - 8400._R8*r(i)**4*z(i)**3        &  
     &       - 3150._R8*r(i)**4*z(i)*rz**2 + 4032._R8*r(i)**2*z(i)**5    &  
     &       + 5040._R8*r(i)**2*z(i)**3*rz**2 + 1260._R8*r(i)**2*z(i)*   &  
     & rz**4                                                             &  
     &       - 128._R8*z(i)**7 - 336._R8*z(i)**5*rz**2                   &  
     &       - 280._R8*z(i)**3*rz**4 - 70._R8*z(i)*rz**6)
      go to 200
   20 continue
!
!....even decapole
      co=pi/(560._R8*rz**8)
      g(i) = co*(7._R8*r(i)**10 - 280._R8*r(i)**8*z(i)**2 - 35._R8*r(i)  &  
     & **8*rz**2                                                         &  
     &       + 1120._R8*r(i)**6*z(i)**4 + 840._R8*r(i)**6*z(i)**2*rz**2  &  
     &       + 70._R8*r(i)**6*rz**4 - 896._R8*r(i)**4*z(i)**6            &  
     &       - 1680._R8*r(i)**4*z(i)**4*rz**2 - 840._R8*r(i)**4*z(i)**2*  &  
     & rz**4                                                             &  
     &       - 70._R8*r(i)**4*rz**6 + 128._R8*r(i)**2*z(i)**8            &  
     &       + 448._R8*r(i)**2*z(i)**6*rz**2 + 560._R8*r(i)**2*z(i)**4*  &  
     & rz**4                                                             &  
     &       + 280._R8*r(i)**2*z(i)**2*rz**6 + 35._R8*r(i)**2*rz**8      &  
     &       - 7._R8*rz**10)
      gr(i)= co*(70._R8*r(i)**9 - 2240._R8*r(i)**7*z(i)**2               &  
     &       - 280._R8*r(i)**7*rz**2 + 6720._R8*r(i)**5*z(i)**4          &  
     &       + 5040._R8*r(i)**5*z(i)**2*rz**2 + 420._R8*r(i)**5*rz**4    &  
     &       - 3584._R8*r(i)**3*z(i)**6 - 6720._R8*r(i)**3*z(i)**4*rz**  &  
     & 2                                                                 &  
     &       - 3360._R8*r(i)**3*z(i)**2*rz**4 - 280._R8*r(i)**3*rz**6    &  
     &       + 256._R8*r(i)*z(i)**8 + 896._R8*r(i)*z(i)**6*rz**2         &  
     &       + 1120._R8*r(i)*z(i)**4*rz**4 + 560._R8*r(i)*z(i)**2*rz**6  &  
     &       + 70._R8*r(i)*rz**8)
      gz(i)= co*(-560._R8*r(i)**8*z(i) + 4480._R8*r(i)**6*z(i)**3        &  
     &       + 1680._R8*r(i)**6*z(i)*rz**2 - 5376._R8*r(i)**4*z(i)**5    &  
     &       - 6720._R8*r(i)**4*z(i)**3*rz**2 - 1680._R8*r(i)**4*z(i)*   &  
     & rz**4                                                             &  
     &       + 1024._R8*r(i)**2*z(i)**7 + 2688._R8*r(i)**2*z(i)**5*rz**  &  
     & 2                                                                 &  
     &       + 2240._R8*r(i)**2*z(i)**3*rz**4 + 560._R8*r(i)**2*z(i)*    &  
     & rz**6)
      grz(i)=co*(-4480._R8*r(i)**7*z(i) + 26880._R8*r(i)**5*z(i)**3      &  
     &       + 10080._R8*r(i)**5*z(i)*rz**2 - 21504._R8*r(i)**3*z(i)**5  &  
     &       - 26880._R8*r(i)**3*z(i)**3*rz**2 - 6720._R8*r(i)**3*z(i)*  &  
     & rz**4                                                             &  
     &       + 2048._R8*r(i)*z(i)**7 + 5376._R8*r(i)*z(i)**5*rz**2       &  
     &       + 4480._R8*r(i)*z(i)**3*rz**4 + 1120._R8*r(i)*z(i)*rz**6)
      gzz(i)=co*(-560._R8*r(i)**8 + 13440._R8*r(i)**6*z(i)**2            &  
     &       + 1680._R8*r(i)**6*rz**2 - 26880._R8*r(i)**4*z(i)**4        &  
     &       - 20160._R8*r(i)**4*z(i)**2*rz**2 - 1680._R8*r(i)**4*rz**4  &  
     &       + 7168._R8*r(i)**2*z(i)**6 + 13440._R8*r(i)**2*z(i)**4*rz**  &  
     & 2                                                                 &  
     &       + 6720._R8*r(i)**2*z(i)**2*rz**4 + 560._R8*r(i)**2*rz**6)
      grr(i)=co*(630._R8*r(i)**8 - 15680._R8*r(i)**6*z(i)**2             &  
     &       - 1960._R8*r(i)**6*rz**2 + 33600*r(i)**4*z(i)**4            &  
     &       + 25200._R8*r(i)**4*z(i)**2*rz**2 + 2100._R8*r(i)**4*rz**4  &  
     &       - 10752._R8*r(i)**2*z(i)**6 - 20160._R8*r(i)**2*z(i)**4*    &  
     & rz**2                                                             &  
     &       - 10080._R8*r(i)**2*z(i)**2*rz**4 - 840._R8*r(i)**2*rz**6   &  
     &       + 256._R8*z(i)**8 + 896._R8*z(i)**6*rz**2                   &  
     &       + 1120._R8*z(i)**4*rz**4 + 560._R8*z(i)**2*rz**6            &  
     &       + 70._R8*rz**8)
      go to 200
  200 continue
      return
  250 continue
!......error
      ineg=39
 
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
