        subroutine tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,itype)
 
!*************************************************************************
!
!  tgmap.f  -  tv80gks mapping routines
!
!  contents    tggetmap  -  Get the linear map used by the plotting library
!        tgsetmap  -  Set the linear map used by the plotting library
!        tgchkmap - check to see if the specified mapping is legal
!        map  -  Define an x and y axes in linear coordinates
!        mapll  -  Define an x and y axes in log-log coordinates
!        mapsl  -  Define an x and y in linear-log coordinates
!        mapls  -  Define an x and y axes in log-linear coordinates
!        maps  -  Define an x and y axes in linear coordinates and
!           draw ticks
!        mapsll  -  Define an x and y axes in log coordinates and draw
!           ticks
!        mapssl  -  Define an x,y axes in linear-log coordinates and
!           draw ticks
!        mapsls  -  Define an x,y axes in log-linear coordinates and
!           draw ticks
!        mapg  -  Define an x and y axes in linear coordinates and draw
!           grid
!        mapgll  -  Define an x and y axes in log coordinates and
!           draw grid
!        mapgsl  -  Define an x,y axes in linear-log coordinates and
!           draw grid
!        mapgls  -  Define an x,y axes in log-linear coordinates and
!           draw grid
!        mapx  -  Select mapping routine at run time
!        mapp  -  Define an x,y axes as polar
!        maplog  -  Convert log coordinates of a map
!        dders  -  Set clipping on or off
!
!  description  The mapping routines let you define x and y axes and the
!               window in the picture plane into which you wish to draw.
!
!    National Energy Research Supercomputer Center
!    Lawrence Livermore National Laboratory
!    University of California
!    Livermore, California 94550, USA
!
!    (C) Copyright 1991 The Regents of the University of California.
!    All Rights Reserved.
!
!*************************************************************************
 
!*************************************************************************
!
!  tggetmap  -  Get the linear map used by the plotting library
!
!  synopsis     call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,itype)
!               real vl,vr      Viewport x range
!               real vb,vt      Viewport y range
!               real wl,wr      Window x range
!               real wb,wt      Window y range
!               integer itype   Type of map
!
!  description  Allows a single routine to get the map that will be used
!               convert the users coordinates into window the workstation
!               window linear units.
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
      INTEGER itype
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   vr,vb,vt,wl,wr,wb,wt,vl
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
 
 
 
        vl = xvmin
        vr = xvmax
        vb = yvmin
        vt = yvmax
        wl = xwmin
        wr = xwmax
        wb = ywmin
        wt = ywmax
        itype = maptyp
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
