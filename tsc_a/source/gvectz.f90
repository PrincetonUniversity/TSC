      subroutine gvectz(rr,z,xi,zi,n,gz,gzz)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!......calculates derivatives wrt first argument
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 z,xi,zi,gz,gzz,rr,a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,c1,c2
      REAL*8 c3,c4,d1,d2,d3,d4,rpxi,rxi,zmzi,rksq,rk,sqrxi,x,ce,ck
      REAL*8 term1,term2
      REAL*8, intrinsic :: log
!============
      dimension xi(n),zi(n),gz(n),gzz(n)
!
      data a0,a1,a2,a3,a4/1.38629436112_R8,9.666344259E-2_R8,            &  
     &3.590092383e-2_R8,3.742563713e-2_R8,1.451196212e-2_R8/
      data b0,b1,b2,b3,b4/.5_R8,.12498593597_R8,6.880248576E-2_R8,       &  
     &3.328355346e-2_R8,4.41787012e-3_R8/
      data c1,c2,c3,c4/.44325141463_R8,6.260601220E-2_R8,                &  
     &4.757383546e-2_R8,1.736506451e-2_R8/
      data d1,d2,d3,d4/.24998368310_R8,9.200180037E-2_R8,                &  
     &4.069697526e-2_R8,5.26449639e-3_R8/
!
      do 100 i=1,n
      rpxi= rr  +xi(i)
      rxi=  rr  *xi(i)
      zmzi= z   -zi(i)
      rksq=4._R8*rxi/(rpxi**2+zmzi**2)
      rk=sqrt(rksq)
      sqrxi=sqrt(rxi)
      x=1._R8-rksq
      ce=1._R8+x*(c1+x*(c2+x*(c3+x*c4)))+                                &  
     &x*(d1+x*(d2+x*(d3+x*d4)))*(-log(x))
      ck=a0+x*(a1+x*(a2+x*(a3+x*a4)))+                                   &  
     &(b0+x*(b1+x*(b2+x*(b3+x*b4))))*(-log(x))
!
      term1=2._R8*ck-2._R8*ce-ce*rksq/x
      term2=2._R8*xi(i)-rksq*rpxi
!
      gz(i)=-rk*0.25_R8*zmzi/sqrxi*term1
      gzz(i)=-rk*0.25_R8/sqrxi*(term1*                                   &  
     &(1.-rksq*zmzi**2/(4.*rxi))+zmzi**2*rksq**2/(4.*rxi*x)              &  
     &*(ce-ck+2._R8*ce*rksq/x))
  100 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
