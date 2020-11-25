      subroutine intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER neqn,kold,ki,kip1,i,j,jm1,limit1,l
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 y,xout,yout,ypout,phi,psi,x,g,w,rho,hi,temp1,term
      REAL*8 psijm1,gamma,eta,temp2,temp3
!============
      dimension y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
      dimension g(13),w(13),rho(13)
      data g(1)/1.0_R8/,rho(1)/1.0_R8/
!============      
      hi=xout-x
      ki=kold+1
      kip1=ki+1
      do 5 i=1,ki
      temp1=i
5     w(i)=1.0_R8/temp1
      term=0.0_R8
      do 15 j=2,ki
      jm1=j-1
      psijm1=psi(jm1)
      gamma=(hi+term)/psijm1
      eta=hi/psijm1
      limit1=kip1-j
      do 10 i=1,limit1
10    w(i)=gamma*w(i)-eta*w(i+1)
      g(j)=w(1)
      rho(j)=gamma*rho(jm1)
15    term=psijm1
      do 20 l=1,neqn
      ypout(l)=0.0_R8
20    yout(l)=0.0_R8
      do 30 j=1,ki
      i=kip1-j
      temp2=g(i)
      temp3=rho(i)
      do 25 l=1,neqn
      yout(l)=yout(l)+temp2*phi(l,i)
25    ypout(l)=ypout(l)+temp3*phi(l,i)
30    continue
      do 35 l=1,neqn
35    yout(l)=y(l)+hi*yout(l)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
