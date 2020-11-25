!#include "f77_dcomplx.h"
      subroutine contur (jjj1, bc, jjj2, ba, max, axarr, imin, imax      &  
     & , istop, ayarr, jmin, jmax, jstop)
 
!**************************************************************************
!
!  contur
!
!  calls        gqcntn
!               gqnt
!               setcrt
!               q7plot
!
!**************************************************************************
 
!======================================================================
! 12/28/2001 ler pppl - change arguments to "call setcrt" to double
!                        precision
!=======================================================================
!
!c purpose: draw a contour plot
!           (it is assumed that ba(i,j) = fcn(axarr(i,j),ayarr(i,j)))
!
!c author: f. n. fritsch
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE tgcblk     
      USE tgchar     
      USE tgcolr     
      USE tgmap      
      USE tgpnts     
      USE tvgxx1
      USE q7quad
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER imin,imax,istop,jjj2,jmin,jmax,jstop,jjj1
      INTEGER ima,nmax,i,i1m1
      INTEGER ie,nt,j,istep,jstep,jma,nstep,nmin
      INTEGER itemz, max
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   wd,vp,czmax2,czmin2
      REAL   REAL
!============
!     common /tvgxx1/ ksent, kclip, krastr, ktyp, cxfact, cxcons
!    * , cyfact, cycons, cxmi, cxma, cymi, cyma, klx, kly, kpol
!    * , knamsb, knerr, crasmx, karrsz, cxarr(256)
!    * , cyarr(256), cxlast, cylast, kspace, kcfx, kcfy
!
!     common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
!    * , int, k01, k02, c01, cdf, czmin, czmax
!
!      dimension ba(2), bc(2), axarr(2), ayarr(2)
      dimension ba(*), axarr(*), ayarr(*), bc(*)
      REAL*8 ba, axarr, ayarr, bc
        dimension wd(4),vp(4)
!
!     equivalence (ima, nmax)
!     equivalence (i, i1m1)
 
!**************************************************************************
!
!  tgcommon  -  TV80 to GKS common blocks
!
!  description  This is used by routines in tv80gks to make a common
!               place for accessing commons.
!
!    National Energy Research Supercomputer Center
!    Lawrence Livermore National Laboratory
!    University of California
!    Livermore, California 94550, USA
!
!    (C) Copyright 1991 The Regents of the University of California.
!    All Rights Reserved.
!
!**************************************************************************
 
!       parameter (MAXPNT=200)
!       integer WRDSIZ
!       parameter (WRDSIZ=8)
!
!       logical dotrans
!       REAL   transx,transy
!       REAL   scalex,scaley
!       REAL   rotat
!       REAL   centrx,centry
!       REAL   trnmat(3,3)
!       logical matmade
!       character*8 tgname
!       logical doinit
!
!       common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
!    +    centrx,centry,trnmat,matmade,tgname,doinit
!
!       REAL   xvmin,xvmax,yvmin,yvmax
!       REAL   xwmin,xwmax,ywmin,ywmax
!       integer maptyp,iclip
!
!       common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
!    +    maptyp,iclip
!
!       REAL   chx,chy
!       integer ichcase
!       integer ichindx
!       REAL   chrot
!       integer ichangle
!       REAL   chparm(4,4)
!       REAL   chupx,chupy
!       integer ichfont
!       REAL   chaddx,chaddy
!       logical autofeed
!
!       common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
!    +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
!    +                  ichclip,minfont,maxfont,autofeed
!
!       REAL   red(16),green(16),blue(16)
!       integer icurclr
!
!       common /tgcolr/ icurclr,red,green,blue
!
!       integer numpnt
!       REAL   xpnt(MAXPNT),ypnt(MAXPNT)
!       integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
!       common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
!    +                  ilnindx,ilnspace,iptspace
 
 
 
! If called by user set routine top routine name
 
      if (tgname .eq. 'NONE') tgname = 'CONTUR'
 
 
 
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
! get contours for negative k01 option...
    1 k01 = -k01
      bc(k01) = bc(2)
      i1m1 = k01-1
      i = i1m1
      cdf = (bc(k01)-bc(1))/REAL(i1m1)
      do 2 j = 2, i1m1
         bc(j) = bc(j-1)+cdf
    2 continue
!
! set up addresses for q7plot..
    3 continue
!
! begin logic for sweeping through grid...
    4 continue
      istep = istop
      ima = imax-istep
      nmax = ima
      jstep = jstop
      jma = jmax-jstep
      nstep = jstep*max
      nmin = (jmin-1)*max
      nmax = nmin+ima
      ima = nmax
      nmin = nmin+imin
!
      do 40 itemz = nmin, nmax, istep
!
! element (i,jmim)
      cxc(1) = axarr(itemz)
      cyc(1) = ayarr(itemz)
      czc(1) = ba(itemz)
!
! element (i+istep,jmin)
      cxc(4) = axarr(itemz+istep)
      cyc(4) = ayarr(itemz+istep)
      czc(4) = ba(itemz+istep)
      if (czc(4) .ge. czc(1)) go to 22
      czmax = czc(1)
      czmin = czc(4)
      go to 23
!
   22 czmax = czc(4)
      czmin = czc(1)
   23 continue
      karrsz = itemz
!
      do 35 j = jmin, jma, jstep
      karrsz = karrsz+nstep
!
! element (i,j+jstep)
      cxc(2) = axarr(karrsz)
      cyc(2) = ayarr(karrsz)
      czc(2) = ba(karrsz)
!
! element (i+istep,j+jstep)
      cxc(3) = axarr(karrsz+istep)
      cyc(3) = ayarr(karrsz+istep)
      czc(3) = ba(karrsz+istep)
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
      cxc(1) = cxc(2)
      cyc(1) = cyc(2)
      czc(1) = czc(2)
      cxc(4) = cxc(3)
      cyc(4) = cyc(3)
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
!
9999  if (tgname .eq. 'CONTUR') tgname = 'NONE'
      return
 
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
