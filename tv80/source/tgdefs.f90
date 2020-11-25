        block data tgdefs
 
!*************************************************************************
!
!  tgdefs - block data for tv80 to gks shell
!
!  description  Sets the common blocks to something meaningful
!
!*************************************************************************
 
 
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
!============
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
 
 
 
! include the site dependent info
 
!***************************************************************************
!
!  tgsite.inc - site dependent information
!
!***************************************************************************
 
! The minimum and maximum font numbers
 
!     integer minfnt,maxfnt
!     parameter (minfnt = 1, maxfnt = 20)
 
! The default font
 
      integer deffnt
      parameter (deffnt = 1)
 
! The character height and expansion factor for fonts 1, 2, 3 and 4
 
!     parameter (chrht1 = .009 , chexp1 = 0.87 )
!     parameter (chrht2 = .013 , chexp2 = 0.90 )
!     parameter (chrht3 = .016 , chexp3 = 0.98 )
!     parameter (chrht4 = .018 , chexp4 = 1.28 )
!
 
!       data doinit /.TRUE./
 
!       data tgname/"NONE    "/
!       data minfont/minfnt/
!       data maxfont/maxfnt/
 
!       data chparm/chrht1,.0156250 ,chexp1,.00781250 ,                  &  
!    &              chrht2,.0235294 ,chexp2,.0117647 ,                   &  
!    &              chrht3,.0312500 ,chexp3,.0156250 ,                   &  
!    &              chrht4,.0476190 ,chexp4,.0238095 /
!
!       data red   /0. ,1. ,1. ,0. ,0. ,0. ,1. ,1. ,                     &  
!    & .0 ,.5 ,.5 ,.0 ,.0 ,.0 ,.5 ,.5 /
!       data green /0. ,1. ,0. ,1. ,0. ,1. ,0. ,1. ,                     &  
!    & .0 ,.5 ,.0 ,.5 ,.0 ,.5 ,.0 ,.5 /
!       data blue  /0. ,1. ,0. ,0. ,1. ,1. ,1. ,0. ,                     &  
!    & .0 ,.5 ,.0 ,.0 ,.5 ,.5 ,.5 ,.0 /
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
