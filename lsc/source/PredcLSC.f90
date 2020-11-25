!#include "f77_dcomplx.h"
!
!     PredcLSC begins                       ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE PredcLSC
      USE dielec
      USE params
      USE PlPr
      USE RayWrk
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL RungeLSC, ftion, prtout
!     These integer variables allow stopping the
!     integration when some condition is met.
!     Lstop    is meant for serious conditions
!              such as going out of bounds.
!     The integer with which PrtOut is called ... call PrtOut(Lprint)...
!     Lprint   tells whether to force printing, or wait till nfreq
!              steps go by.  Lprint=1  --> force print
!                            Lprint=0  --> advance counter and decide
!     The include statements are to catch: begin,h,nstep,lstop,
!     and especially NEQS and NEQSP1
      INTEGER i,j, jstart, BoundsEr, iBndsErr, nBndsErr
      REAL*8    yok(NEQSP1), thi, tho
      REAL*8 pc(NEQS), c(NEQS),p(NEQS),dumy(NEQS)
      REAL*8 k1,k2,k3,k4,k5,k6,k7,k8,k9,k10
      DATA      k1,         k2,        k3,         k4,          k5,      &  
     &          k6,         k7,        k8,         k9,          k10 /    &  
     &  2.65277776_R8, -1.4861111_R8, 1.5138888_R8,  -.34722222_R8,      &  
     &   .94266666_R8,                                                   &  
     &   .34722222_R8,  1.2638888_R8,  .59722222_R8,  .125_R8,           &  
     &   .05733333_R8/
      DATA  iBndsErr, nBndsErr / 0 , 25 /
      REAL*8 AREAL

      jstart = 4
      iscatplt = 1
 10   continue
      do 20 i = 1, NEQS
      pc(i) = 0.0_R8
 20   continue
      BoundsEr = 0
      lstop = 0
      call RungeLSC(BoundsEr)
      if (BoundsEr .ne. 0 ) then
        call LSCendr(' Cant recover in PredcLSC')
        return
      endif
 
      call KeepGood(y, yok, NEQSP1)
 
      do 80 j = jstart, nstep
      call KeepGood(y, yok, NEQSP1)
      if ( iError .ge. 1 ) return
      if ( lstop .eq. 1) go to 90
      y(NEQSP1) = begin + AREAL(j)*HstpLH
      do 30 i=1,NEQS
      dumy(i) = (y3(i)+y3(i) + y2(i))/3._R8
      p(i)=dumy(i)+HstpLH*(k1*f(i)+k2*f3(i)+k3*f2(i)+k4*f1(i))
      y1(i)=y2(i)
      y2(i)=y3(i)
      y3(i)=y(i)
! modify
      y(i)=p(i)-k5*pc(i)
      f1(i) = f2(i)
      f2(i) = f3(i)
      f3(i) = f(i)
30    continue
      call ftion(BoundsEr)
        if(BoundsEr .ne. 0) go to 100
      do 40 i = 1,NEQS
! correct
      c(i)= dumy(i)+HstpLH*(k6*f(i)+k7*f3(i)+k8*f2(i)+k9*f1(i))
      pc(i) = p(i)-c(i)
! final
      y(i) = c(i)+k10*pc(i)
40    continue
 
60    call prtout(0)
70    call ftion(BoundsEr)
        if(BoundsEr .ne. 0) go to 100
80    continue
      return
90    call prtout(1)
      return
!     -                                 Attempt to recover from ray error,
!     -                                 most likely running out of bounds,
!     -                                 by specular reflection off the last
!     -                                 good set of values.
 100  continue
      jstart = j
!      if (EnhcNpar .eq. TRUE) then
!        call BounShft(yok, y, NEQSP1, NparEnhc)
!      else
!        call BounceIt(yok, y, NEQSP1)
!      endif
      if ( ScatKdeg .gt. 0.001_R8) then
        call BounShft(yok, y, NEQSP1, ScatKdeg, thi, tho)
        inciThet(iscatplt) = thi/3.1416_R8*180._R8
        scatThet(iscatplt) = tho/3.1416_R8*180._R8
        iscatplt=iscatplt+1
        if (iscatplt .gt. NPLTDIM) iscatplt=NPLTDIM
 
      else
        call BounceIt(yok, y, NEQSP1)
      endif
 
      iBndsErr = iBndsErr + 1
      if (iBndsErr .ge. nBndsErr) then
        iBndsErr = 0
        write(nLSCcomm, '(                                               &  
     &  '' Rays out of bounds; recovered'',i4,'' times'')') nBndsErr
      endif
      goto 10
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
