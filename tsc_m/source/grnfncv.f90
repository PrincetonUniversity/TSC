      subroutine grnfncv(ineg,nmult,xt,zt,xs,zs,n,ans)
!.....4.71 grnfncv
!
!......calculates poloidal flux at point (xt,zt) due to current
!......at location (xs,zs)  . . . returns value in ans
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nmult,n,ineg,i,imult,imultp
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xt,zt,xs,zs,ans,q0,q1,q2,q3,q4,q5,r0,r1,r2,r3,r4,r5,pi
      REAL*8 tpi,zr,qk2,qn,qlg,bk,rz,co
!============
      dimension xt(n),zt(n),xs(n),zs(n),ans(n)
      data q0,q1,q2,q3,q4,q5/-.30685281944_R8, .29822748722_R8,          &  
     &   .00367617100_R8, -.01091055499_R8, .00860373511_R8,             &  
     &  .00725598106_R8/
      data r0,r1,r2,r3,r4,r5/ .25_R8, .06250928488_R8,                   &  
     &  .00489241049_R8, .01034604435_R8, .01358621540_R8,               &  
     &  .00220893506_R8/
      data pi,tpi/3.1415926535_R8,6.283185308_R8/
!============      
      INTEGER :: istat = 0 
!============      
!
      do 100 i=1,n
      zr = zs(i) - zt(i)
      qk2 = 4._R8*xs(i)*xt(i)/((xs(i)+xt(i))**2 + zr**2)
      qn = 1.0_R8- qk2
      qlg = -log(qn)
      bk = q0 + qn*(q1+qn*(q2+qn*(q3+qn*(q4+qn*q5))))                    &  
     &   + (r0+qn*(r1+qn*(r2+qn*(r3+qn*(r4+qn*r5)))))*qlg
      ans(i) = -sqrt(xt(i)*xs(i)/qk2)*bk*2._R8
  100 continue
!.....check for multipolar coils
      if(nmult .le. 0) return
      do 200 i=1,n
      if(xs(i) .lt. 100._R8) go to 200
      rz = zs(i)
      imult = int(xs(i) - 100._R8)
      if(imult .lt. 0 .or. imult.gt.10) go to 250
      imultp = imult + 1
      go to(10,11,12,13,14,15,16,17,18,19,20),imultp
   10 continue
!
!....even nullapole
      ans(i) = tpi*rz**2
      go to 200
   11 continue
!
!....odd nullapole
      ans(i) = 0._R8
      go to 200
   12 continue
!
!....even dipole
      ans(i) = tpi*(xt(i)**2 - rz**2)/2._R8
      go to 200
   13 continue
!
!....odd dipole
      ans(i) = tpi*(xt(i)**2*zt(i))/rz
      go to 200
   14 continue
!
!....even quadrapole
      ans(i)=pi*(xt(i)**4-4._R8*xt(i)**2*zt(i)**2-2._R8*xt(i)**2*rz**2+  &  
     & rz**4)/                                                           &  
     &          (4*rz**2)
      go to 200
   15 continue
!
!....odd quadrapole
      ans(i)=pi*xt(i)**2*zt(i)*(3._R8*xt(i)**2-4*zt(i)**2-3._R8*rz**2)   &  
     &      /(3*rz**3)
      go to 200
   16 continue
!
!.....even hexapole
      ans(i)=(tpi/3._R8)*(xt(i)**6-12._R8*xt(i)**4*zt(i)**2-3*xt(i)**4*  &  
     & rz**2                                                             &  
     &               + 8*xt(i)**2*zt(i)**4 + 12*xt(i)**2*zt(i)**2*rz**2  &  
     &               + 3._R8*xt(i)**2*rz**4 - rz**6 )/(8._R8*rz**4)
      go to 200
   17 continue
!
!.....odd hexapole
      co=pi/(30._R8*rz**5)
      ans(i)=co*(15._R8*xt(i)**6*zt(i) - 60._R8*xt(i)**4*zt(i)**3        &  
     &       - 30._R8*xt(i)**4*zt(i)*rz**2 + 24._R8*xt(i)**2*zt(i)**5    &  
     &       + 40._R8*xt(i)**2*zt(i)**3*rz**2 + 15._R8*xt(i)**2*zt(i)*   &  
     & rz**4)
      go to 200
   18 continue
!
!....even octapole
      co=pi/(160._R8*rz**6)
      ans(i)=co*(5._R8*xt(i)**8-120._R8*xt(i)**6*zt(i)**2-20._R8*xt(i)**  &  
     & 6*rz**2                                                           &  
     &       + 240._R8*xt(i)**4*zt(i)**4 + 240._R8*xt(i)**4*zt(i)**2*    &  
     & rz**2                                                             &  
     &       + 30._R8*xt(i)**4*rz**4 - 64._R8*xt(i)**2*zt(i)**6          &  
     &      -160._R8*xt(i)**2*zt(i)**4*rz**2-120._R8*xt(i)**2*zt(i)**2*  &  
     & rz**4                                                             &  
     &       - 20._R8*xt(i)**2*rz**6 + 5._R8*rz**8)
      go to 200
   19 continue
!
!....odd octapole
      co=pi/(140._R8*rz**7)
      ans(i)=co*xt(i)**2*zt(i)*(35._R8*xt(i)**6 - 280._R8*xt(i)**4*zt(i)  &  
     & **2                                                               &  
     &       - 105._R8*xt(i)**4*rz**2 + 336._R8*xt(i)**2*zt(i)**4        &  
     &       + 420._R8*xt(i)**2*zt(i)**2*rz**2 + 105._R8*xt(i)**2*rz**4  &  
     &       - 64._R8*zt(i)**6 - 168._R8*zt(i)**4*rz**2                  &  
     &       - 140._R8*zt(i)**2*rz**4 - 35._R8*rz**6)
      go to 200
   20 continue
!
!....even decapole
      co=pi/(560._R8*rz**8)
      ans(i)=co*(7._R8*xt(i)**10-280._R8*xt(i)**8*zt(i)**2-35._R8*xt(i)  &  
     & **8*rz**2                                                         &  
     &       + 1120._R8*xt(i)**6*zt(i)**4 + 840._R8*xt(i)**6*zt(i)**2*   &  
     & rz**2                                                             &  
     &       + 70._R8*xt(i)**6*rz**4 - 896._R8*xt(i)**4*zt(i)**6         &  
     &      -1680._R8*xt(i)**4*zt(i)**4*rz**2-840._R8*xt(i)**4*zt(i)**2*  &  
     & rz**4                                                             &  
     &       - 70._R8*xt(i)**4*rz**6 + 128._R8*xt(i)**2*zt(i)**8         &  
     &      +448._R8*xt(i)**2*zt(i)**6*rz**2+560._R8*xt(i)**2*zt(i)**4*  &  
     & rz**4                                                             &  
     &       + 280._R8*xt(i)**2*zt(i)**2*rz**6 + 35._R8*xt(i)**2*rz**8   &  
     &       - 7._R8*rz**10)
      go to 200
  200 continue
      return
  250 continue
!......error
      ineg=39
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
