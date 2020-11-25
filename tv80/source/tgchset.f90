        subroutine tgchset (ifont,index,iorient)
 
 
!**************************************************************************
!
!  tggks.f  -  Routines to initialize the gks drivers
!
!  contents    The routines in this file open GKS if it is not already
!        opened and then open a specific workstation given a
!        workstation identifier and a connection id.  The names
!        and the workstation that is opened will change depending
!        on what workstations are available with your implementation
!        of GKS.  After initializing a workstation, tginit should
!        be called to put the workstation in a state that is
!        consistent with what the library thinks it should be.
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
 
!**************************************************************************
!
!  cgmbin -  CGM binary
!
!  description    Open GKS if necessary, get the identity of the workstation
!        set the output file name, open and then activate the
!        workstation.
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
      INTEGER iorient,ifont
      INTEGER index
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   pi
!============
        REAL   size,upx,upy
        parameter (PI = 3.1415926535897932 )
 
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
 
#ifndef NCAR_DUMMY 
 
        ichindx = index
        ichfont = ifont
 
! font with character precision
 
        call gstxfp (ichfont,2)
 
! Set the rotation based on the angle in iorient
 
        if (iorient .eq. 0) then
          chrot = 0. 
          chupx = cos (chrot + PI/2. )
          chupy = sin (chrot + PI/2. )
        if (dotrans) then
          upx = cos (chrot + rotat + PI/2. )
          upy = sin (chrot + rotat + PI/2. )
        else
          upx = chupx
          upy = chupy
        endif
          call gschup (upx,upy)
          size = chparm(1,ichindx+1)
          if (dotrans) size = size * scaley
          call gschh (size)
          chaddx = - chparm(2,ichindx+1) * chupx
          chaddy = - chparm(2,ichindx+1) * chupy
          size = chparm(3,ichindx+1)
          if (dotrans) size = size * scalex/scaley
          call gschxp (size)
        else if (iorient .eq. 1) then
          chrot = PI/2. 
          chupx = cos (chrot + PI/2. )
          chupy = sin (chrot + PI/2. )
        if (dotrans) then
          upx = cos (chrot + rotat + PI/2. )
          upy = sin (chrot + rotat + PI/2. )
        else
          upx = chupx
          upy = chupy
        endif
          call gschup (upx,upy)
          size = chparm(1,ichindx+1)
          if (dotrans) size = size * scalex
          call gschh (size)
          chaddx = - chparm(2,ichindx+1) * chupx
          chaddy = - chparm(2,ichindx+1) * chupy
          size = chparm(3,ichindx+1)
          if (dotrans) size = size * scaley/scalex
          call gschxp (size)
        endif
#endif
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
