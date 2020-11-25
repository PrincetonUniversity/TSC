        subroutine mapp (rmax,xmin,xmax,ymin)
 
!*************************************************************************
!
!  mapp  -  Define an x,y axes as polar
!
!  synopsis     call mapp (rmax)
!               call mapp (rmax,xmin)
!               call mapp (rmax,xmin,xmax)
!               call mapp (rmax,xmin,xmax,ymin)
!
!               float rmax                      Max radius
!               float xmin,xmax,ymin            Window coordinates
!
!  description  rmax is the maximum radial distance.  If positive, a grid
!               is drawn.  xmin and xmax define the range in NDC where the
!               plot will be placed.  ymin is the minimum y in NDC.  ymax
!               is ymin+(xmax-xmin) since the plot must be square.  The
!               default size for xmin,xmax,ymin,ymax is .001 to .999
!
!*************************************************************************
 
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE tgcblk     
      USE tgchar     
      USE tgcolr     
      USE tgmap      
      USE tgpnts     
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     INTEGER maxpnt,ichclip,minfont,maxfont
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   xvl,yvb,xvr,xwr,ywt,xwl,ywb,yvt
!============
        REAL*8 rmax
        REAL*8 xmin,xmax,ymin
        integer tgchkmap
 
! include the standard common block
 
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
 
 
 
! Set library entry routine name
 
        if (tgname .eq. 'NONE') tgname = 'MAPP'
 
 
 
        xvl = 0.001 
        yvb = 0.001 
        xvr = 0.999 
 
        xwr = abs(rmax)
        ywt = xwr
        xwl = -xwr
        ywb = -xwr
 
        if (ymin .ne. -1.0 ) yvb = ymin
        if (xmax .ne. -1.0 ) xvr = xmax
      if (xmin .ne. -1.0 ) xvl = xmin
 
        yvt = yvb+(xvr-xvl)
 
! Although this is a polar mapping, check it as if it were linear
 
        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,1) .ne. 0) then
          goto 9999
        endif
 
        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,5)
 
        if (rmax .gt. 0) call maplab (0. ,xwr,0,2,1)
 
9999    if (tgname .eq. 'MAPP') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
