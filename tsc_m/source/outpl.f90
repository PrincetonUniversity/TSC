!#include "f77_dcomplx.h"
      subroutine outpl
!......6.80 outpl
!
!.....generates plots of the solution
!
      USE CLINAM
      USE SCR1
      USE SCR11
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER irs1sw,iig,ig,i,iabs,igp,lendsk,isd,idisk,ios17
      INTEGER icount,irec,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xmin,xmax,ymin,ymax,yval,x1,y1,x2,y2,psione
      REAL*8 AREAL, sum
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: apsi
      INTEGER :: istat
!============
      data irs1sw/1/
!
      if(.NOT.ALLOCATED(apsi)) ALLOCATE (apsi(pnx,pnz), STAT=istat)
      if (istat .ne. 0 ) stop "Allocation Error : outpl" 

      if((imovie.gt.0 .and. imovie.lt.10)) go to 489
!..rxw:
      if(noplot(11).gt.0) goto 1236
      if( isurf.ne.0 .and. lrswtch.eq.0.and.ineg.ge.0) call surfplot
      if( isurf.ne.0 .and. lrswtch.eq.0.and.ineg.ge.0) call surfplot3
      if( isurf.ne.0 .and. lrswtch.eq.0.and.ineg.ge.0                    &  
     &   .and.noplot(89).eq.0._R8) call surfplt2
      if( isurf.ne.0 .and. lrswtch.eq.0.and. (idens.eq.0 .or. idens.eq.3)) call densplot
      if( isurf.ne.0 .and. lrswtch.eq.0.and. ialpha.eq.1) call           &  
     & alphaplot
      if( isurf.ne.0 .and. lrswtch.eq.0 .and. ineg.eq.0 .and. iimp.gt.1)  &  
!    &                                                                   &  
     &        call impplot
      if(irippl.gt.0 .and. lrswtch.eq.0.and.ineg.eq.0) call ripplot
 1236 continue
!
      if(ncoil.le.0) go to 489
!
!......output values of coil currents
!
!
      write(nout,1086) kcycle,times
 1086 format(" coil curr (ka)  at cycle= ",i7, " time ",1pe13.5,// ,     &  
     &  "  group   actual sum   prepro sum  voltage   prepro v ",        &  
     &  "  gcur0(ka)   gcurfb(ka) ")
!      if(acoef(901).ge.1.)write(nterm,10860)
!10860 format(2x,"group",7x,"I_g(kA)",7x,
!     1   "J_c(kA/cm**2)")
      call groupcur(grsum,grsum0,gvsum,gvsum0,gcur0ka,gcurfka)
      sum=0._R8
      do 187 iig=1,ngroupt
      ig=nogroupt(iig)
      if(grsum(ig)+gcur0ka(ig)+gcurfka(ig).eq.0) go to 187
      write(nout,1087) ig,grsum(ig),grsum0(ig),gvsum(ig),gvsum0(ig),     &  
     &                 gcur0ka(ig),gcurfka(ig)
 1087 format(i5,1p6e12.4)
!      if(acoef(901).ge.1.) then
!      ccdens=abs(grsum(ig))/(1.+isym)*1.e-4/(dxcoil(ig)*dzcoil(ig))
!      write(nterm,10870) ig,grsum(ig)/(1.+isym),ccdens
!10870 format(i5,10x,1p1e12.4,15x,1pe12.2)
!10870 format(i5,4x,0pf12.2,5x,0pf12.2)
!      sum=sum+(grsum(ig)/(1.+isym))**2
!      endif
  187 continue
!      if(acoef(901).ge.1.) then
!      cnorm=sqrt(sum)
!      write(nterm,10871) cnorm
!10871 format(" norm of coil current vector [kA]=",f12.2)
!      endif
!
!
!.....  new coil currents plot
!..rxw:
      if(noplot(5).gt.0) goto 1231
!
      xmin = 0._R8
      xmax = AREAL(nwire+1)
      if(xmax .le. xmin) go to 489
      ymin = 0._R8
      ymax = 1._R8
      do 700 i=1,nwire
      yval = AREAL(iabs(igroupw(i)))
      ymax = max(ymax,yval)
  700 continue
      ymax = ymax+1
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.656_R8,.801_R8,1.0_R8)
      call setcrt(0.5_R8,ymin)
      do 701 i=1,nwire
      x1 = AREAL(i)-0.5_R8
      y1 = AREAL(iabs(igroupw(i)))
      call vector(x1,y1)
      x2 = AREAL(i)+0.5_R8
      y2 = AREAL(iabs(igroupw(i)))
      call vector(x2,y2)
      call vector(x2,ymin)
  701 continue
      call setold(xmin-(xmax-xmin)*.15_R8,ymin,1,0,1,1)
      write(s100,1701)
      call gtext(s100,80,0)
 1701 format("     group no ")
!
      call setold(xmax,ymax - .05_R8*(ymax-ymin),1,0,1,0)
      write(s100,2701)
      call gtext(s100,80,0)
 2701 format("  group currents(ka) ")
      do 706 ig=1,pngroup-1,2
      igp =ig + 1
      if(grsum(ig).eq.0 .and. grsum(igp).eq.0) go to 706
      write(s100,2706) ig,grsum(ig),igp,grsum(igp)
      call gtext(s100,80,0)
 2706 format(i3,f9.2,i5,f9.2)
  706 continue
!
!
      ymin = 0._R8
      ymax = 0._R8
      do 702 i=1,nwire
      yval = ccoil(ncoil-nwire+i)*udsi*.001_R8
      ymax = max(ymax,yval)
      ymin = min(ymin,yval)
  702 continue
      if(ymax.le.ymin) go to 488
      ymax = 1.05_R8*ymax
      ymin = 1.05_R8*ymin
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.656_R8,.543_R8,.741_R8)
      call setcrt(0.5_R8,0.0_R8)
      do 703 i=1,nwire
      x1 = AREAL(i)-0.5_R8
      y1 = ccoil(ncoil-nwire+i)*udsi*.001_R8
      call vector(x1,y1)
      call vector(x1+1._R8,y1)
  703 continue
      call setold(xmin-(xmax-xmin)*.15_R8,ymin,1,0,1,1)
      write(s100,1702)
      call gtext(s100,80,0)
 1702 format(" actual current(ka)")
!
      ymin = 0._R8
      ymax = 0._R8
      do 705 i=1,nwire
      yval = (ccoil(ncoil-nwire+i)-cwire0(i))*udsi*.001_R8
      ymax = max(ymax,yval)
      ymin = min(ymin,yval)
  705 continue
      if(ymax.le.ymin) go to 488
      ymax = ymax*1.05_R8
      ymin = ymin*1.05_R8
      call mapg(xmin,xmax,ymin,ymax,.142_R8,.656_R8,.285_R8,.483_R8)
      call setcrt(0.5_R8,0.0_R8)
      do 704 i=1,nwire
      x1 = AREAL(i)-0.5_R8
      y1 = (ccoil(ncoil-nwire+i)-cwire0(i))*udsi*.001_R8
      call vector(x1,y1)
      call vector(x1+1._R8,y1)
  704 continue
      call setold(xmin-(xmax-xmin)*.15_R8,ymin,1,0,1,1)
      write(s100,1705)
      call gtext(s100,80,0)
 1705 format("induced current(ka)")
!
  488 continue
      call frscj(7)
  489 continue
!..rxw:
 1231 continue
!
!....special plot of terms in ohms law
!.....December, 1987  KML
      if(irfp.ne.1 .or. lrswtch.ne.0) go to 493
      if(irs1sw.eq.1 .and. irst1.eq.1) go to 493
      if(kcycle.le.0) go to 493
      call lingplt
      call lingpl2
      call lingpl3
!
 493  continue
!
!.....plot poloidal flux function psi
!..rxw:
      if(noplot(6).gt.0) goto 1232
      call cplot(psi,1)
      if((imovie.gt.0 .and. imovie.lt.10)) go to 498
!
!.....plot toroidal current jphi/x
!..rxw:
      if(noplot(91).gt.0) goto 1233
      call cplot(ajphi,3)
 1233 continue
      if(icplgf.ne.0) call cplot(g,8)
      if(lrswtch .ne. 0) go to 495
!
!.....write disk file for summary plot
      if(ifrst(4).eq.0) go to 100
      ifrst(4) = 0
      lendsk = ((ncycle-kcycle)/nskipl+2)*pnx*pnz
      isd = 0
      idisk = 1
!     ifiles2(1:6) = 'spplsc'
!     ifiles2(7:7) = isuffix(1:1)
      if( numargs .lt. 1 ) then
         ifiles2 = 'spplsc' // isuffix(1:1)
      else
         ifiles2 = 'spplsc' // '.' // trim(suffix)
      end if
      open(nsc2,file=trim(ifiles2),status='unknown',form='unformatted',  &  
     &     iostat=ios17)
      nplrec = 0
      icount = 0
  100 continue
      nplrec = nplrec+1
      if(nplrec .ge. 250) ineg=40
      iminy(nplrec) = iminn
      imaxy(nplrec) = imaxx
      jminy(nplrec) = jminn
      jmaxy(nplrec) = jmaxx
      if(iplim.lt.0) jminy(nplrec) = 2
      if(iplim.lt.0) jmaxy(nplrec) = nzp
      irec = pnx*pnz
      psione = psi(1,1)
      psi(1,1) = psilim
      psi(1,2) = iplim
      call bufout(nsc2,psi(1,1),psi(pnx,pnz))
      psi(1,1) = psione
      icount = icount+irec
!
!
!.....special x-point plot
!..rxw:
      if(noplot(8).gt.0) goto 1234
      if(icplxp .gt. 0) call xptplt
 1234 continue
!.....special divertor heat distribution plot
!..rxw:
      if(noplot(9).gt.0) goto 1235
      if(idiv.ne.0) call dplplt
 1235 continue
  495 continue
!
!
!.....optional plots
!.....December,1987         KML
      if(irfp.eq.1) go to 496
      if(icplet.ne.0) call cplot(etay,2)
!     if(icplet.ne.0) call cplot(etay,10)
      if(icplwf.ne.0 .and. kcycle.gt.0) call cplot(w,7)
      if(icplpr.ne.0) call cplot(pr,5)
      if(icplbv.ne.0 .and. kcycle.gt.0) call cplot(b,4)
      if(icpluv.ne.0 .and. kcycle.gt.0) call cplot(u,6)
      go to 497
 496  continue
      if(irs1sw.eq.1 .and. irst1.eq.1) go to 497
      if(icplet.ne.0 .and. kcycle.gt.0) call cplot(pling1,2)
      if(icplwf.ne.0 .and. kcycle.gt.0 ) call cplot(pling4,7)
      if(icplpr.ne.0) call cplot(pr,5)
      if(icplbv.ne.0 .and. kcycle.gt.0 ) call cplot(pling2,4)
      if(icpluv.ne.0 .and. kcycle.gt.0 ) call cplot(pling3,6)
!
 497  continue
      if(kcycle.lt.0) go to 498
      if(ivplbp.ne.0) call vplot(1)
      if(ivplvi.ne.0) call vplot(2)
      if(ivplfr.ne.0) call vplot(3)
      if(ivpljp.ne.0) call vplot(4)
      if(ivplvc.ne.0) call vplot(5)
      if(ivplvt.ne.0) call vplot(6)
      if(ivplvi.ne.0.or.ivplfr.ne.0.or.ivplvc.ne.0) call vplot(7)
!
!
  498 continue
!
      if(isym.eq.1) go to 1232
      do 2232 i=1,nxp
      apsi(i,nh) = 0._R8
      do 2232 j=1,nh-1
      apsi(i,nh+j) = (psi(i,nh+j) - psi(i,nh-j))
      apsi(i,nh-j) = (psi(i,nh-j) - psi(i,nh+j))
 2232 continue
      call cplot(apsi,9)
 1232 continue
      irs1sw = 0
!     DEALLOCATE (apsi)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
