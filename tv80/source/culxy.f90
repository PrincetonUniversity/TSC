        subroutine culxy (x,y,xnew,ynew)
 
!*************************************************************************
!
!  culxy  -  Convert user coordinates to linear-linear
!
!  synopsis     call culxy (x,y,xnew,ynew)
!               real x,y        Users x,y
!               real xnew,ynew  Transformed point
!
!  description  someday
!
!*************************************************************************
 
!       real x,y,xnew,ynew
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
      REAL   ynew,xnew
!============
      REAL*8 x,y
        REAL   xtemp,ytemp
 
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
 
 
 
        if (maptyp .eq. 1) then
 
          xtemp = x
          ytemp = y
 
        else if (maptyp .eq. 2) then
 
          xtemp = x
          if (y .le. 0.0 ) then
            print *,'Error in CULXY: Invalid log coordinate'
            ytemp = ywmin
          else
!           ytemp = alog10(y)
           ytemp = dlog10(y)
          endif
 
        else if (maptyp .eq. 3) then
 
          if (x .le. 0.0 ) then
            print *,'Error in CULXY: Invalid log coordinate'
            xtemp = xwmin
          else
!           xtemp = alog10(x)
           xtemp = dlog10(x)
          endif
          ytemp = y
 
        else if (maptyp .eq. 4) then
 
          if (x .le. 0.0 .or. y .le. 0.0 ) then
            print *,'Error in CULXY: Invalid log coordinate'
            xtemp = xwmin
            ytemp = ywmin
          else
!           xtemp = alog10(x)
           xtemp = dlog10(x)
!            ytemp = alog10(y)
            ytemp = dlog10(y)
          endif
 
        else if (maptyp .eq. 5) then
 
          xtemp = x * cos(y)
          ytemp = x * sin(y)
 
        endif
 
        xnew = xtemp
        ynew = ytemp
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
