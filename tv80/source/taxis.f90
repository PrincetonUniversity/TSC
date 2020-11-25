      subroutine taxis (fxb, fyb, fww)
!**************************************************************************
!
!  taxis.F - axis drawing companion for gaxis
!
!  contents    taxis - draws axis
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
 
!=======================================================================
! 12/28/2001 ler pppl - dimension cliprect
!=======================================================================
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE tgcblk     
      USE tgchar     
      USE tgcolr     
      USE tgmap      
      USE tgpnts     
      USE tvaxis
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER itype,ierr,j
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL   cxb,cyb,cxe,cye,ctickx,cticky,coffx,coffy,cac,cbc,cpx
      REAL   cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2
!============
      REAL*8 fxb, fyb, fww, dcxb, dcyb, gxp, gxe, gyp, gye
!c created: may17.1976 by mick archuleta
!
!c revised: sep26.1980 by steven williams
!
! this routine can be called after a call to gaxis has been made.
! its purpose is to draw an axis line without tick marks. the
! information for the axis line is still available from the last call
! to gaxis. thus, fxb and fyb are used to specify the begin point
! of this tick line. the end point is calculated by taxis and it
! is in the same relative direction as the line drawn by gaxis. fww is
! a variable which tells taxis how long to make the tick
! marks and in which direction. if the same direction and length, then
! it should be 1, if the opposite direction and same length, then
! it should be -1.
!
!c variable declarations:
!
!     logical klr
!     common /tvaxis/ cxb, cyb, cxe, cye, ctickx, cticky
!    * , coffx, coffy, cac, cbc, kori, ksize, kn1, klr
!    * , klxsav, klysav, cpx(30), cpy(30), kdiv, kepclp
!
! xdiff,ydiff   - difference between this and last point
! fxbl,fybl     - linear version of fxb,fyb
! cxbl,cybl     - linear version of cxb,cyb
!
      REAL   xdiff,ydiff
      REAL   fxbl,fybl,cxbl,cybl
 
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
 
      REAL   cliprect(4)
 
 
 
! Set library entry routine name
 
        if (tgname .eq. 'NONE') tgname = 'taxis'
 
!
! Convert the starting point from taxis and gaxis to linear and find the
! difference in the x and y directions.
!
      call culxy (fxb,fyb,fxbl,fybl)
      dcxb = cxb
      dcyb = cyb
      call culxy (dcxb,dcyb,cxbl,cybl)
      xdiff = fxbl - cxbl
      ydiff = fybl - cybl
!
! Get world and NDC coordinates as well as map type, then set the
! map type to be linear.
!
      call tggetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,1)
!
! Turn clipping off so tic marks outside of clipping rect will be drawn
!
      call gqclip (ierr,kepclp,cliprect)
      call gsclip (0)
!
! Loop through the array cpx and cpy to place tic marks at each location.
!
      do 80 j = 1, kdiv
         gxp = cpx(j) + xdiff
         gyp = cpy(j) + ydiff
         call line (gxp, gyp, gxp+cac*ctickx*fww, gyp+cbc*cticky*fww)
   80 continue
!
! Reset the map and draw the axis
!
      call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call line (fxb, fyb, fxb+cxe-cxb, fyb+cye-cyb)
!
! Return Clipping back to its original setting
!
      call gsclip (kepclp)
!
      if (tgname .eq. 'taxis') tgname = 'NONE'
      return
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
