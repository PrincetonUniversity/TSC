!
!     -----------------------------------------------------------------
!
      SUBROUTINE RungeLSC(BoundsEr)
      USE dielec
      USE params
      USE PlPr
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ftion, prtout
 
      INTEGER ia,ib,ic,  i,k
      INTEGER BoundsEr
      REAL*8 q(NEQSP1),a(4),b(4),c(4)
      REAL*8 dum
      DATA (a(i),i=1,4)/ 0.5_R8, 0.29289322_R8, 1.70710675_R8,           &  
     &                   0.1666667_R8/
      DATA (b(i),i=1,4)/ 2._R8, 1._R8, 1._R8, 2._R8/
      DATA (c(i),i=1,4)/ 0.5_R8, 0.29289322_R8, 1.70710675_R8, 0.5_R8/
      f(NEQSP1) = 1._R8
      q(NEQSP1) =0._R8
      call ftion(BoundsEr)
        if(BoundsEr .ne. 0) go to 100
      do 10 i=1,NEQS
      y1(i) = y(i)
      f1(i) = f(i)
10    q(i)=0._R8
 
      do 40 ia=2,4
      do 13 ib=1,4
12    call ftion(BoundsEr)
        if(BoundsEr .ne. 0) go to 100
      do 13 ic=1,NEQSP1
      dum = a(ib)*(f(ic)-b(ib)*q(ic))
      y(ic) = y(ic) + HstpLH*dum
      q(ic) = q(ic) + 3._R8*dum-c(ib)*f(ic)
13    continue
      k = ia-1
      call prtout(1)
15    go to (20,30,40) k
20    do 22 ic= 1,NEQS
      f2(ic) = f(ic)
      y2(ic) = y(ic)
22    continue
      go to 40
30    do 35 ic=1,NEQS
      y3(ic) = y(ic)
      f3(ic) = f(ic)
35    continue
40    continue
      return
 100  continue
      BoundsEr = 1
      call LSCwarn(' Boundary error in RungeLSC')
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
