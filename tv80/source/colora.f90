        subroutine colora (colorstr)
 
!*************************************************************************
!
!  colora  -  Set the current color by name
!
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
      character*(*) colorstr
 
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
 
        if (tgname .eq. 'NONE') tgname = 'COLORA'
 
 
 
! search for the string that represents the color.
 
        if (colorstr(1:5) .eq. 'black' .or.                              &  
     &      colorstr(1:5) .eq. 'BLACK') then
          icurclr = 0
        else if (colorstr(1:5) .eq. 'white' .or.                         &  
     &      colorstr(1:5) .eq. 'WHITE') then
          icurclr = 1
        else if (colorstr(1:3) .eq. 'red' .or.                           &  
     &      colorstr(1:3) .eq. 'RED') then
          icurclr = 2
        else if (colorstr(1:5) .eq. 'green' .or.                         &  
     &      colorstr(1:5) .eq. 'GREEN') then
          icurclr = 3
        else if (colorstr(1:4) .eq. 'blue' .or.                          &  
     &      colorstr(1:4) .eq. 'BLUE') then
          icurclr = 4
        else if (colorstr(1:4) .eq. 'cyan' .or.                          &  
     &      colorstr(1:4) .eq. 'CYAN') then
          icurclr = 5
        else if (colorstr(1:7) .eq. 'magenta' .or.                       &  
     &      colorstr(1:7) .eq. 'MAGENTA') then
          icurclr = 6
        else if (colorstr(1:6) .eq. 'yellow' .or.                        &  
     &      colorstr(1:6) .eq. 'YELLOW') then
          icurclr = 7
        else
          call tgerror (1)
          if (tgname .eq. 'COLORA') tgname = 'NONE'
          return
        endif
 
        call tgflush
 
#ifndef NCAR_DUMMY
        call gsplci (icurclr)
        call gspmci (icurclr)
        call gstxci (icurclr)
        call gsfaci (icurclr)
#endif
 
        if (tgname .eq. 'COLORA') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
