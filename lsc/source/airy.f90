!c
!         compute airy functions
      SUBROUTINE airy(ai,bi,aip,bip,x)
!         x < 0
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nm, m, n
      REAL*8    ai, bi, aip, bip, x, a1, a2, x3, xx, f, ft, g, gt
      REAL*8    c, cn, sn, c0, c1, c2, f1, f2
      REAL*8    aa(2),bb(2)
      DATA    a1/0.355028_R8/,a2/0.2588194_R8/
      if(x.gt.7._R8) x=7._R8
      if(abs(x).gt.7._R8) go to 100
!         small argument expansion
      nm=3+2.5_R8*abs(x)
      do 20 m=1,2
      xx=x+0.01_R8*(m-1.5_R8)
      x3=xx**3
      f=1._R8
      ft=f
      g=xx
      gt=g
      do 10 n=1,nm
      f=f*x3/(3._R8*n*(3._R8*n-1._R8))
      g=g*x3/(3._R8*n*(3._R8*n+1._R8))
      ft=ft+f
   10 gt=gt+g
      aa(m)=a1*ft-a2*gt
   20 bb(m)=sqrt(3._R8)*(a1*ft+a2*gt)
      go to 500
!         asymptotic expansion
  100 do 120 m=1,2
      xx=-x-0.01_R8*(m-1.5_R8)
      c=2._R8/3._R8*xx**1.5_R8
      cn=cos(c)
      sn=sin(c)
      c0=0.3989423_R8*xx**(-0.25_R8)
      c1=0.0694444_R8/c
      c2=0.0371335_R8/c**2
      f1=1._R8-c1-c2
      f2=1._R8+c1-c2
      aa(m)=c0*(cn*f1+sn*f2)
  120 bb(m)=c0*(cn*f2-sn*f1)
  500 ai=0.5_R8*(aa(1)+aa(2))
      bi=0.5_R8*(bb(1)+bb(2))
      aip=100._R8*(aa(2)-aa(1))
      bip=100._R8*(bb(2)-bb(1))
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
