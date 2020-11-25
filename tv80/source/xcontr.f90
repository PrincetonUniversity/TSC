!#include "f77_dcomplx.h"
      subroutine xcontr (ii, jjj1, bc, jjj2, ba, max, imin, imax         &  
     & , istop, jmin, jmax, jstop, axarr, ayarr)
 
!**************************************************************************
!
!  xcontr
!
!  calls        gqcntn
!               gqnt
!               setcrt
!               q7plot
!
!  description  This routine was designed to be called by the contour routines
!               only.
!
!**************************************************************************
 
!======================================================================
! 12/28/2001 ler pppl - change arguments to "call setcrt" to double
!                        precision
!=======================================================================
!
!c author: f. n. fritsch
!
!c notes:
!
! 1) it is assumed that a(i,j) = fcn(axarr(i),ayarr(j)).
!
!c variable declarations:
!
!============
      USE q7quad
      USE tvgxx1
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER jjj1,jjj2,imin,imax,istop,jmin,jmax,jstop,ii
      INTEGER i,i1m1,ie,nt,j,istep,ima
      INTEGER jstep,jma,nstep,nmin,itemz
      INTEGER max
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   vp,wd,czmax2,czmin2
      REAL   REAL
!============
!     common /tvgxx1/ ksent, kclip, krastr, ktyp, cxfact, cxcons
!    * , cyfact, cycons, cxmi, cxma, cymi, cyma, klx, kly, kpol
!    * , knamsb, knerr, crasmx, karrsz, cxarr(256)
!    * , cyarr(256), cxlast, cylast, kspace, kcfx, kcfy
!
!     dimension ba(2), bc(2), axarr(2), ayarr(2)
      dimension ba(*), bc(*), axarr(*), ayarr(*)
      REAL*8 ba, bc, axarr, ayarr
       dimension wd(4),vp(4)
!
!     common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
!    * , int, k01, k02, c01, cdf, czmin, czmax
!
!     equivalence (i,i1m1)
      save
 

#ifndef NCAR_DUMMY
!
! set min and max NDC window
!
        call gqcntn (ie,nt)
        call gqnt (nt,ie,wd,vp)
        cxmi = vp(1)
        cxma = vp(2)
        cymi = vp(3)
        cyma = vp(4)
!
!c program statements:
!
      call setcrt (0.0d0 , 0.0d0 )
      kp = 0
      k02 = jjj2
      k01 = jjj1
      if (k01) 1, 50, 3
!
!  get contours for negative k01 option...
    1 k01 = -k01
      bc(k01) = bc(2)
      i1m1 = k01-1
      i = i1m1
      cdf = (bc(k01)-bc(1))/REAL(i1m1)
      do 2 j = 2, i1m1
         bc(j) = bc(j-1)+cdf
    2 continue
!
!  set up addresses for q7plot..
    3 continue
!
!  begin logic for sweeping through grid...
    4 continue
      istep = istop
      ima = imax-istep
      jstep = jstop
      jma = jmax-jstep
      nstep = jstep*max
      nmin = (jmin-1)*max
      do 40  i = imin, ima, istep
      if (ii .eq. 1) go to 10
      cxc(1) = axarr(i)
      cxc(2) = axarr(i)
      cxc(3) = axarr(i+istep)
      cxc(4) = axarr(i+istep)
      go to 15
!
   10 cxc(1) = i
      cxc(2) = i
      cxc(3) = i+istep
      cxc(4) = i+istep
!
   15 continue
!
! element (i,jmin)
      czc(1) = ba(nmin+i)
!
! element (i+istep,jmin)
      czc(4) = ba(nmin+i+istep)
      if (czc(4) .ge. czc(1)) go to 22
      czmax = czc(1)
      czmin = czc(4)
      go to 23
!
   22 czmax = czc(4)
      czmin = czc(1)
   23 continue
      itemz = nmin+i
      do 35 j = jmin, jma, jstep
      if (ii .eq. 1) go to 2350
      cyc(1) = ayarr(j)
      cyc(2) = ayarr(j+jstep)
      cyc(3) = ayarr(j+jstep)
      cyc(4) = ayarr(j)
      go to 2375
!
 2350 cyc(1) = j
      cyc(4) = j
      cyc(2) = j+jstep
      cyc(3) = j+jstep
 2375 continue
      itemz = itemz+nstep
!
! element (i,j+jstep)
      czc(2) = ba(itemz)
!
! element (i+istep,j+jstep)
      czc(3) = ba(itemz+istep)
      if (czc(3) .ge. czc(2)) go to 24
      czmax2 = czc(2)
      czmin2 = czc(3)
      go to 25
!
   24 czmax2 = czc(3)
      czmin2 = czc(2)
!
   25 if (czmax .ge. czmax2) go to 26
      czmax = czmax2
!
   26 if (czmin2 .ge. czmin) go to 27
      czmin = czmin2
!
   27 call q7plot (bc(1))
      if (kp .eq. 1) go to 9999
!
   34 czc(1) = czc(2)
      czc(4) = czc(3)
      czmax = czmax2
      czmin = czmin2
   35 continue
   40 continue
!
!     plotting completed...
   30 continue
      go to 9999
!
!     set up for k01 = 0 option...
   50 c01 = bc(1)
      cdf = bc(2)
      if (cdf .gt. 0) go to 4
!
!     error exit if  cdf .le. 0 ...
      print *,'error 3: cdf .le. 0'
#endif
!
9999  return
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
