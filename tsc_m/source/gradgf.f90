      subroutine gradgf(ineg,nmult,xt,zt,xs,zs,gradx,gradz)
!.....4.80 gradgf
!
!.....calculates gradient of the poloidal flux at
!.....point (xt,zt) due to current at location (xs,zs)
!
!
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nmult,ineg
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xt,zt,xs,zs,gradx,gradz,pi,ak0,ak1,ak2,ak3,ak4,bk0,bk1
      REAL*8 bk2,bk3,bk4,ae1,ae2,ae3,ae4,be1,be2,be3,be4,tpi,zd,ck2
      REAL*8 x1,x2,x3,x4,ck,elipk,elipe,term1,term2,fac,ans
!============
      data pi/3.1415926535897_R8/,                                       &  
     &    ak0 /1.38629436112_R8/,                                        &  
     &    ak1 /0.09666344259_R8/,                                        &  
     &    ak2 /0.03590092383_R8/,                                        &  
     &    ak3 /0.03742563713_R8/,                                        &  
     &    ak4 /0.01451196212_R8/,                                        &  
     &    bk0 /0.5_R8/,                                                  &  
     &    bk1 /0.12498593597_R8/,                                        &  
     &    bk2 /0.06880248576_R8/,                                        &  
     &    bk3 /0.03328355346_R8/,                                        &  
     &    bk4 /0.00441787012_R8/,                                        &  
     &    ae1 /0.44325141463_R8/,                                        &  
     &    ae2 /0.0626060122_R8/,                                         &  
     &    ae3 /0.04757383546_R8/,                                        &  
     &    ae4 /0.01736506451_R8/,                                        &  
     &    be1 /0.2499836831_R8/,                                         &  
     &    be2 /0.09200180037_R8/,                                        &  
     &    be3 /0.04069697526_R8/,                                        &  
     &    be4 /0.00526449639_R8/
!
!
      tpi = 2._R8*pi
      zd = zs - zt
      ck2 = 4._R8*xs*xt/((xs+xt)**2 + zd**2)
      x1 = 1._R8- ck2
      x2 = x1*x1
      x3 = x2*x1
      x4 = x3*x1
      ck = sqrt(ck2)
!
      elipk = ak0+ak1*x1+ak2*x2+ak3*x3+ak4*x4                            &  
     &             - (bk0+bk1*x1+bk2*x2+bk3*x3+bk4*x4)*log(x1)
!
      elipe=1.0_R8
      if(abs(x1) .gt. 1.0E-6_R8)                                         &  
     & elipe=1.0_R8+ae1*x1+ae2*x2+ae3*x3+ae4*x4                          &  
     &         - (be1*x1+be2*x2+be3*x3+be4*x4)*log(x1)
!
      term1 = -2._R8*elipk
      term2 = (2._R8-ck2)/(1._R8-ck2)*elipe
      fac = 1._R8/(4._R8*sqrt(xs*xt))/ck*(term1+term2)
      call gf(ineg,nmult,xt,zt,xs,zs,ans)
      gradx = ans/2._R8/xt - fac*(2._R8*xs - ck2*(xt+xs))
      gradz = -fac*ck2*(zs-zt)
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
