!
!----------------------------------------------------------/invert
!         invert 'aa' - apply same operations to unit matrix as to 'aa'
      SUBROUTINE invert(idim,nguide,aa,ai)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER idim, nguide, l,k,l1
      COMPLEX*16 f,f0,aa(idim,1),ai(idim,1)
!         set 'ai' equal to unit matrix
      do 100 l=1,nguide
      do 100 k=1,nguide
      ai(l,k)=(0._R8,0._R8)
  100 ai(l,l)=(1._R8,0._R8)
!
      do 500 l=1,nguide
!
      if(l.eq.nguide) go to 180
      do 150 l1=l+1,nguide
      f=aa(l,l)+aa(l1,l)
      f0=(1._R8,0._R8)
      if(abs(f).lt.abs(aa(l,l))) f0=(-1._R8,0._R8)
      do 120 k=1,nguide
      aa(l,k)=aa(l,k)+aa(l1,k)*f0
  120 ai(l,k)=ai(l,k)+ai(l1,k)*f0
  150 continue
!         set aa(l,l) = 1
  180 f=aa(l,l)
      do 200 k=1,nguide
      ai(l,k)=ai(l,k)/f
  200 aa(l,k)=aa(l,k)/f
!         set aa(l1,l) = 0  ,  l1 /= l
      do 400 l1=1,nguide
      if(l1.eq.l) go to 400
      f=aa(l1,l)
      do 300 k=1,nguide
      ai(l1,k)=ai(l1,k)-f*ai(l,k)
      aa(l1,k)=aa(l1,k)-f*aa(l,k)
  300 continue
  400 continue
  500 continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
