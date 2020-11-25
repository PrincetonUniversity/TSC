      subroutine outpr1
!......6.41 outpr1
!
!.....output dg=g(i,nh+j)-g(i,nh-j+1) etc. to check symmetry of all
!.....variables.
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nxp1,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 denom,dg,dpr,dr,dajey,dq,dxsqoj,dbmagy,dx,dz,detay
      REAL*8 dpsia,dajphi,dgs,domeg,dw,dabig,du,db
!============
  102 format (1x,i3,11(1x,1pe10.3))
  120 format (1h1,1x,"i= ",i4,2x,"cycle=",i7,2x,"time=",1pe12.4,         &  
     &        2x,"dt=",1pe12.4)
  122 format (/1x,"  j    dg         dpr        dr         dajey"        &  
     &       ,"      dq         dxsqoj     dbmagy")
  124 format (/1x,"  j    dx         dz         detay      dpsi",        &  
     &        "       dajphi     dgs        domeg       dw         ",    &  
     &        "dabig       du         db")
!
      if(isym.ne.0) return
      nxp1=nxp+1
!
      do 280 i=1,nxp1
      write (nout,120) i,kcycle,time,dt
      write (nout,122)
      do 250 j=1,nh
      denom=max(abs(g(i,nh+j)+g(i,nh-j+1)),1.0E-12_R8)
      dg=(g(i,nh+j)-g(i,nh-j+1))/denom
      denom=max(abs(pr(i,nh+j)+pr(i,nh-j+1)),1.0E-12_R8)
      dpr=(pr(i,nh+j)-pr(i,nh-j+1))/denom
      denom=max(abs(r(i,nh+j)+r(i,nh-j+1)),1.0E-12_R8)
      dr=(r(i,nh+j)-r(i,nh-j+1))/denom
      denom=max(abs(ajey(i)+ajey(i)),1.0E-12_R8)
      dajey=(ajey(i)-ajey(i))/denom
      denom=max(abs(q(i,nh+j)+q(i,nh-j+1)),1.0E-12_R8)
      dq=(q(i,nh+j)-q(i,nh-j+1))/denom
      denom=max(abs(xsqoj(i)+xsqoj(i)),1.0E-12_R8)
      dxsqoj=(xsqoj(i)-xsqoj(i))/denom
      denom=max(abs(bmagy(i,nh+j)+bmagy(i,nh-j+1)),1.0E-12_R8)
      dbmagy=(bmagy(i,nh+j)-bmagy(i,nh-j+1))/denom
      write (nout,102) j,dg,dpr,dr,dajey,dq,dxsqoj,dbmagy
  250 continue
!
      write (nout,124)
      do 260 j=1,nh-1
      denom=max(abs(xary(i)+xary(i)),1.0E-12_R8)
      dx=(xary(i)-xary(i))/denom
      denom=max(abs(zary(nh+j)-zary(nh-j)),1.0E-12_R8)
      dz=(zary(nh+j)+zary(nh-j))/denom
      denom=max(abs(etay(i,nh+j)+etay(i,nh-j)),1.0E-12_R8)
      detay=(etay(i,nh+j)-etay(i,nh-j))/denom
      denom=max(abs(psi(i,nh+j)+psi(i,nh-j)),1.0E-12_R8)
      dpsia=(psi(i,nh+j)-psi(i,nh-j))/denom
      denom=max(abs(ajphi(i,nh+j)+ajphi(i,nh-j)),1.0E-12_R8)
      dajphi=(ajphi(i,nh+j)-ajphi(i,nh-j))/denom
      denom=max(abs(gs(i,nh+j)+gs(i,nh-j)),1.0E-12_R8)
      dgs=(gs(i,nh+j)-gs(i,nh-j))/denom
      denom=max(abs(omeg(i,nh+j)+omeg(i,nh-j+1)),1.0E-12_R8)
      domeg=(omeg(i,nh+j)-omeg(i,nh-j+1))/denom
      denom=max(abs(w(i,nh+j)-w(i,nh-j+1)),1.0E-12_R8)
      dw=(w(i,nh+j)+w(i,nh-j+1))/denom
      denom=max(abs(abig(i,nh+j)-abig(i,nh-j)),1.0E-12_R8)
      dabig=(abig(i,nh+j)+abig(i,nh-j))/denom
      denom=max(abs(u(i,nh+j)+u(i,nh-j+1)),1.0E-12_R8)
      du=(u(i,nh+j)-u(i,nh-j+1))/denom
      denom=max(abs(b(i,nh+j)-b(i,nh-j)),1.0E-12_R8)
      db=(b(i,nh+j)+b(i,nh-j))/denom
      write (nout,102) j,dx,dz,detay,dpsia,dajphi,dgs,domeg,dw,          &  
     &              dabig,du,db
  260 continue
  280 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
