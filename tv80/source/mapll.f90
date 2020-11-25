        subroutine mapll (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)    
 
!*************************************************************************
!
!  mapll  -  Define an x and y axes in log-log coordinates
!
!  synopsis     call mapll (xleft,xright, ybot,ytop)
!               call mapll (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)
!
!               float xleft,xright,ybot,ytop    World coordinates
!               float xmin,xmax,ymin,ymax       Window coordinates
!
!  description  The world coordinates are mapped into the range specified
!               by the Window.  The window x and y range should be between
!               0.0 and 1.0.  0.0 and 1.0 is the default value.
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
      REAL   xvl,yvb,xvr,yvt,xwl,ywb,xwr,ywt
!============
        REAL*8 xleft,xright,ybot,ytop
        REAL*8 xmin,xmax,ymin,ymax
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
 
        if (tgname .eq. 'NONE') tgname = 'MAPLL'
 
 
 
      xvl = 0.0 
      yvb = 0.0 
      xvr = 1.0 
      yvt = 1.0 
      if (xmin .ne. -1.0 ) xvl = xmin
      if (ymin .ne. -1.0 ) yvb = ymin
      if (xmax .ne. -1.0 ) xvr = xmax
      if (ymax .ne. -1.0 ) yvt = ymax
 
        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop
 
        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,4) .ne. 0) then
          goto 9999
        endif
 
        call maplog (xwl,xwr)
        call maplog (ywb,ywt)
 
        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,4)
 
9999    if (tgname .eq. 'MAPLL') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
