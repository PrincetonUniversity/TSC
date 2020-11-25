      SUBROUTINE PULLn (VPULL,KCHOP,PF,NCHOP,VC,DC,VPS,V01)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER kchop,nchop
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pf,vc,dc,vps,v01,vpull,aipull,f,v02,c,pk2,vp,dc1,t1,t2
      REAL*8 t3
!============
      aipull=abs(PF)*1.E3_R8
!
      if(KCHOP.eq.1)then
!
!------------/DC modulation ratio/ -----
      DC=0.694_R8+2.012E-2_R8*VC-2.755E-3_R8*VC**2+                      &  
     &  3.136E-4_R8*VC**3+1.434E-5_R8*VC**4-2.238E-6_R8*VC**5
!-------------------------------------------
!
!--------------/ f=frequency/----
      f=2944-76.97_R8*VC**2+2.97_R8*VC**4-0.0779_R8*VC**6+               &  
     &  8.98E-4_R8*VC**8-3.61E-6_R8*VC**10
!-----------------
!
!----------/ V01=avearge push voltage/----
      V01=0.9584_R8*VPS*( 1._R8+0.146_R8*(1._R8-VPS/600._R8)*(aipull+    &  
     & 1500._R8)/                                                        &  
     &  (aipull+300._R8) )
!
!----------/ VPULL=avearge pull voltage/----
      v02=-94.5_R8*(aipull/NCHOP)**0.25_R8
      c=120.E-6_R8
      pk2=0.9_R8
      vp=0.5_R8*pk2*c*v02*NCHOP*f/((1._R8-DC)*aipull)
!     print *,' pk2 *c *v02 *NCHOP *f',pk2,c,v02,NCHOP,f
      dc1=1-DC
!     print *,' vc  VPS',vc,VPS
!     print *,'dc1 aipull',dc1,aipull
!ccc  VPULL=v02*(1.+(0.5*pk2*c*v02*NCHOP*f)/((1.-DC)*aipull))
!ccc  if(VPULL.gt.0.)VPULL=0.
      t1=DC/f
      t2=-pk2*c*v02*NCHOP/aipull
      t3=(1._R8-DC)/f-t2
      if(t3.ge.0)VPULL=v02*(0.5_R8*t2+t3)/(t2+t3)
      if(t3.lt.0._R8)VPULL=0.5_R8*v02*(1._R8-DC)/(f*t2)
!     print *,' v02 DC f aipull VPULL vp',v02,DC,f,aipull,VPULL,vp
!---o-----------------
!     pause
      end if
!000
      if(KCHOP.eq.2)then
!------------/DC modulation ratio/ -----
      DC=0.628_R8+1.77E-2_R8*VC-1.965E-3_R8*VC**2+                       &  
     &  4.253E-4_R8*VC**3+1.027E-5_R8*VC**4-2.901E-6_R8*VC**5
!-------------------------------------------
!
!--------------/ f=frequency/----
      f=3039.1_R8-69.33_R8*VC**2+2.44_R8*VC**4-0.0619_R8*VC**6+          &  
     &  6.9E-4_R8*VC**8-2.69E-6_R8*VC**10
!-----------------
!
!----------/ V01=avearge push voltage/----
      V01=VPS
!
!----------/ VPULL=avearge pull voltage/----
      if(abs(aipull/NCHOP).lt.700._R8)then
      v02=-0.75_R8*(aipull/NCHOP)
      else
      v02=-(1575._R8*aipull/NCHOP-6.E5_R8)**0.476_R8
      end if
      c=60.E-6_R8
      pk2=1.8_R8
!ccc  VPULL=v02*(1.+(0.55*pk2*c*v02*NCHOP*f)/((1.-DC)*aipull))
!ccc  if(VPULL.gt.0.)VPULL=0.
      t1=DC/f
      t2=-pk2*c*v02*NCHOP/aipull
      t3=(1._R8-DC)/f-t2
      if(t3.ge.0)VPULL=v02*(0.5_R8*t2+t3)/(t2+t3)
      if(t3.lt.0._R8)VPULL=0.5_R8*v02*(1._R8-DC)/(f*t2)
!--------------------
      end if
      RETURN
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
