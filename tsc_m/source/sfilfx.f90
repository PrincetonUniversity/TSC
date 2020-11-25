      subroutine sfilfx(rf,zf,at,rp,zp,fx)
!
!
!
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 zf,at,rp,zp,fx,rf,c1,c2,z,r,p2,p,a2,h,v,eka,ekb,eea
      REAL*8 eeb,y,ce,ck
!============
  733 c1=12.56637E-7_R8
  742 if(rp) 746,744,746
  744 fx=0.0_R8
      return
!
  746 continue
      if(at.eq.0.0_R8) go to 775
      c2=c1*at
  750 z=(zp-zf)/rf
  752 r=rp/rf
  754 p2=(1.0_R8+r)*(1.0_R8+r)+z*z
  756 p=sqrt (p2)
  758 a2=p2-4.0_R8*r
      if(a2.le.0.0_R8) go to 773
!
  760 h=4.0_R8*r/p2
  761 v=a2/p2
!
 1552 eka=(((0.0145119621_R8*v+0.0374256371_R8)*v+0.0359009238_R8)*v+    &  
     & 0.0966634426_R8                                                   &  
     &)*v+1.38629436_R8
 1554 ekb=(((0.00441787012_R8*v+0.0332835535_R8)*v+0.0688024858_R8)*v+   &  
     & 0.1249859360_R8                                                   &  
     &)*v+0.5_R8
 1556 eea=(((0.0173650645_R8*v+0.0475738355_R8)*v+0.0626060122_R8)*v+    &  
     & 0.443251415_R8)                                                   &  
     &*v+1.0_R8
 1558 eeb=(((0.00526449639_R8*v+0.0406969753_R8)*v+0.0920018004_R8)*v+   &  
     & 0.249983683_R8                                                    &  
     &)*v
!
  765 y=log(v)
  766 ce=eea-y*eeb
  768 ck=eka-y*ekb
  772 fx=p*((1.0_R8-0.5_R8*h)*ck-ce)*c2*rf
      go to 775
  773 fx=sign(1.0E20_R8,at)
  775  continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
