        subroutine mapx (imap,xleft,xright,ybot,ytop,                    &  
     &                   xmin,xmax,ymin,ymax)
 
!*************************************************************************
!
!  mapx  -  Select mapping routine at run time
!
!  synopsis     call mapx (imap,xleft,xright, ybot,ytop)
!               call mapx (imap,xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)
!
!               integer imap                    Mappimg routine
!               float xleft,xright,ybot,ytop    World coordinates
!               float xmin,xmax,ymin,ymax       Window coordinates
!
!  description  The world coordinates are mapped into the range specified
!               by the Window.  The window x and y range should be between
!               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
!               routine is very similar to SPPS SET.  The difference is
!               the last argument to set.  A 1 will specify linear x and
!               y mapping.
!
!               imap Routine  imap Routine  imap Routine
!               ------------  ------------  ------------
!               1    map      5    mapg     9    maps
!               2    mapll    6    mapgll   10   mapsll
!               3    mapsl    7    mapgsl   11   mapssl
!               4    mapls    8    mapgls   12   mapsls
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
        integer imap
        REAL*8 xleft,xright,ybot,ytop
        REAL*8 xmin,xmax,ymin,ymax
        REAL*8 wl,wr,wb,wt,vl,vr,vb,vt
 
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
 
        if (tgname .eq. 'NONE') tgname = 'MAPX'
 
 
 
        if (imap .lt. 1 .or. imap .gt. 12) then
          call tgerror (11)
          goto 9999
        endif
 
        if (imap .gt. 4) then
          vl = 0.11328 
          vb = 0.11328 
          vr = 1.0 
          vt = 1.0 
        else
          vl = 0.0 
          vb = 0.0 
          vr = 1.0 
          vt = 1.0 
        endif
      if (xmin .ne. -1.0 ) vl = xmin
      if (ymin .ne. -1.0 ) vb = ymin
      if (xmax .ne. -1.0 ) vr = xmax
      if (ymax .ne. -1.0 ) vt = ymax
 
        wl = xleft
        wb = ybot
        wr = xright
        wt = ytop
 
        goto (1,2,3,4,5,6,7,8,9,10,11,12),imap
 
1       call map (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
2       call mapll (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
3       call mapsl (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
4       call mapls (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
5       call mapg (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
6       call mapgll (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
7       call mapgsl (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
8       call mapgls (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
9       call maps (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
10      call mapsll (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
11      call mapssl (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999
 
12      call mapsls (wl,wr,wb,wt,vl,vr,vb,vt)
 
9999    if (tgname .eq. 'MAPX') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
