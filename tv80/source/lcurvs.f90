!#include "f77_dcomplx.h"
      subroutine lcurvs (jjj1, bc, jjj2, fcn, nx, ny)
 
!**************************************************************************
!
!  lcurvs
!
!  calls        gqcntn
!               gqnt
!               setcrt
!               q7plot
!
!**************************************************************************
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE tgcblk     
      USE tgchar     
      USE tgcolr     
      USE tgmap      
      USE tgpnts     
      USE tvgxx1
      USE tvgxx2
      USE q7quad
      USE tvmap1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER jjj2,nx,ny,jjj1
      INTEGER i,i1m1,it,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   dx,dy,czmacx,czmin2
      REAL   REAL
!============
      REAL*8 fcn
      external fcn
!=======================================================================
! 12/28/2001 ler pppl - explicitly define passed function 'fcn' as
!                        double precision
!                     - change first argument in first 'fcn' invocation
!                        from single to double precision
!=======================================================================
!
!c notes:
!
! 1) fourth argument (fcn) should be typed generic according
!    to the original documentation in ltss chapter 304.
!    however, for tv80cray, it is assumed that it is an ascii string
!
!c variable declarations:
!
!% None in tvgxx1 are used
!
!     common /tvgxx1/ ksent, kclip, krastr, ktyp, cxfact, cxcons
!    * , cyfact, cycons, cxmi, cxma, cymi, cyma, klx, kly, kpol
!    * , knamsb, knerr, crasmx, karrsz, cxarr(256)
!    * , cyarr(256), cxlast, cylast, kspace, kcfx, kcfy
!
!     common /tvgxx2/ kxarr(511), kyarr(511), kitcnt, klchar
!    * , kntenv, kdashv
!    * , knumch, kcoff, kcpt, kcprm(7,4), klindx, klindy
!    * , kchdx, kchdy, kindyo
!    * , klabel(4,7), kltv80(9), krasmx, kact
!    * , kdstat(8), kname(8), klenfi(8), klevd(8)
!    * , kdbuf(8), kdbufs(8), knfrm(8)
!    * , kepflg(8), kmxfil(8), kminbf(8), klerad(8)
!    * , kdfact(8), knumfi(8)
!
!     common /tvmap1/ kia, knarg, cx1, cx2, cy1, cy2
!
!     common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
!    * , int, k01, k02, c01, cdf, czmin, czmax
!
!      dimension bc(2)
      dimension bc(*)
      REAL*8 bc, x, y, xx, yy, dcx1
 
!        dimension wd(4),vp(4)
!
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
 
      if (tgname .eq. 'NONE') tgname = 'LCURVS'
 
 
 
!
! get min and max NDC window coordinates and user coordinates
!
      call tggetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,it)
!
! The user coordinate should not really be in log form.  But
! if they are, tggetmap will return exponents.  Convert the
! exponents into coordinates.
!
      if (it .eq. 2 .or. it .eq. 4) then
        cy1 = 10.0 **cy1
      cy2 = 10.0 **cy2
      endif
      if (it .eq. 3 .or. it .eq. 4) then
        cx1 = 10.0 **cx1
      cx2 = 10.0 **cx2
      endif
!
!c program statements:
!
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
      dx = (cx2-cx1)/nx
      dy = (cy2-cy1)/ny
      y = cy1
      do 40 i = 1, ny
      cyc(1) = y
      cyc(2) = y+dy
      cyc(3) = y+dy
      cyc(4) = y
      yy = y+dy
      dcx1 = cx1
      czc(1) = fcn(dcx1,y)
      czc(2) = fcn(dcx1,yy)
      if (czc(2) .ge. czc(1)) go to 22
      czmax = czc(1)
      czmin = czc(2)
      go to 23
!
   22 continue
      czmax = czc(2)
      czmin = czc(1)
   23 continue
      x = cx1
      do 35 j = 1, nx
      cxc(1) = x
      cxc(2) = x
      cxc(4) = x+dx
      cxc(3) = cxc(4)
      xx = cxc(4)
      czc(4) = fcn(xx,y)
      czc(3) = fcn(xx,yy)
      if (czc(3) .ge. czc(4)) go to 24
      czmacx = czc(4)
      czmin2 = czc(3)
      go to 25
!
   24 continue
      czmacx = czc(3)
      czmin2 = czc(4)
   25 continue
      if (czmax .ge. czmacx) go to 26
      czmax = czmacx
   26 if (czmin2 .ge. czmin) go to 27
      czmin = czmin2
   27 call q7plot (bc(1))
      if (kp .eq. 1) go to 9999
      czc(1) = czc(4)
      czc(2) = czc(3)
      czmax = czmacx
      czmin = czmin2
      x = xx
   35 continue
      y = yy
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
!     error exit if cdf .le. 0 ...
      print *,'error 3: cdf .le. 0'
!
9999  if (tgname .eq. 'LCURVS') tgname = 'NONE'
      return
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
