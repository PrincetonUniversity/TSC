      subroutine der ( t, y, yp ,neqn)
!
!.....10.8 der
!
      USE CLINAM
      USE BALCLI
      USE DERCOM

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER neqn,nord,ing,mth,ll,l,ith,n,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 y,yp,t,t0
      REAL*8 t1,t2,t3,sym,tt,term
!============
      dimension y(neqn),yp(neqn)
!============      
!     common/dercom/factb(4),
!    1          alfafb(4,ppsi),betafb(4,ppsi),gamafb(4,ppsi),
!    2          alfa0(ppsi),beta0(ppsi),gama0(ppsi),ingold
      if(ifcount.le.0) ingold=0
!
!..... nearest grid.
!
      nord = 2*(npsit-4)
      ing = t/dthe
      if(ing.eq.ingold .and. ifcount.ne.0) go to 13
      mth = (1+isym)*nthe
      t0 = (ing-1)*dthe
      t1 = (ing  )*dthe
      t2 = (ing+1)*dthe
      t3 = (ing+2)*dthe
      do 10 ll=0,3
      l = ll+1
      ith = mod(ing+ll,mth)
      if(ith.le.0) ith = ith+mth
      sym = 1._R8
      if(isym.eq.1 .and. ith.gt.nthe+1) go to 11
      go to 12
   11 ith = 2*(nthe+1)-ith
      sym = -1._R8
   12 continue
!
      tt = (ing-1+ll)*dthe
      do 10 n=1,nord,2
      j = (n+1)/2+2
      alfafb(l,j) = alfav(ith,j)+sym*alfa1(ith,j)*tt                     &  
     &          + alfa2(ith,j)*tt**2
      betafb(l,j) = sym*betav(ith,j)+beta1(ith,j)*tt
      gamafb(l,j) = gamav(ith,j)+sym*gama1(ith,j)*tt
   10 continue
   13 ingold = ing
!
!
!.....define hermite cubic interpolations coefficients
      factb(4) = (t-t1)**2*(t-t2)*.5_R8/dthe**3
      factb(1) = -(t-t2)**2*(t-t1)*.5_R8/dthe**3
      term = (t-t1)**2*(3._R8*t2-t1-2*t)/dthe**3
      factb(2) = 1._R8-term-factb(4)
      factb(3) = term-factb(1)
!
      do 19 j=3,npsit-2
      alfa0(j)= 0._R8
      beta0(j)= 0._R8
      gama0(j)= 0._R8
   19 continue
      do 20 l=1,4
      do 20 j=3,npsit-2
      alfa0(j)= alfa0(j)+ factb(l)*alfafb(l,j)
      beta0(j)= beta0(j)+ factb(l)*betafb(l,j)
      gama0(j)= gama0(j)+ factb(l)*gamafb(l,j)
   20 continue
!
!
!
      do 30 n=1,nord,2
      j = (n+1)/2 + 2
      yp(n) = ( beta0(j)*y(n) + y(n+1) ) / alfa0(j)
!
      yp(n+1) = -((alfa0(j)*gama0(j) + beta0(j)**2)*y(n)                 &  
     &        + beta0(j)*y(n+1))/alfa0(j)
   30 continue
!
!
      ifcount = ifcount + 4
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
