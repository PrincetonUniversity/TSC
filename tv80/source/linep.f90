!#include "f77_dcomplx.h"
        subroutine linep (x1,y1, x2,y2, iexponent)
 
!*************************************************************************
!
!  linep  -  Draw a line with points
!
!  synopsis     call linep (x1,y1, x2,y2)
!               call linep (x1,y1, x2,y2, iexponent)
!
!               real x1,y1              First point
!               real x2,y2              Second
!               integer iexponent       Exponent to use for screen res
!
!  description  Draws a line from x1,y1 to x2,y2 using a series of points.
!               The points are assumed to be a distance of 2**iexponent
!               apart.  Where the screen is assumed to be 2**10 x 2**10
!               (1024 x 1024).  Initially iexponent is assumed to be 2.
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
      INTEGER iexponent,k,numpts
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   deltax,deltay
      REAL   REAL
!============
      REAL*8 x1, y1, x2, y2
        REAL   x1p,y1p,x2p,y2p
      REAL*8 dx1p, dy1p, dx2p, dy2p, ddeltax, ddeltay
        REAL   vl,vr,vb,vt,wl,wr,wb,wt
        REAL   distx,disty,dist
        integer it
 
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
 
        if (tgname .eq. 'NONE') tgname = 'LINEP'
 
 
 
! Convert to plotter coordinates
 
        call culxy (x1,y1,x1p,y1p)
        call clnxy (x1p,y1p,x1p,y1p)
        call culxy (x2,y2,x2p,y2p)
        call clnxy (x2p,y2p,x2p,y2p)
 
! Save the old mapping and reset the mapping to normalized coordinates
 
        call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
        call tgsetmap (vl,vr,vb,vt,vl,vr,vb,vt,1)
 
! find the distance between successive points by calculating distance
! between endpoints and dividing by 2**k.
 
        distx = x2p - x1p
        disty = y2p - y1p
        dist = sqrt (distx*distx + disty*disty)
      if (iexponent .gt. 0) iptspace = iexponent
        k = iptspace
        numpts = nint (dist / (2. **k / 1024. ))
        if (numpts .gt. 0) then
          deltax = distx / REAL(numpts)
          deltay = disty / REAL(numpts)
        endif
 
! Draw the line of points.  To avoid any round off error, the last point
! is explicitly plotted rather than calculated from a series of additions
! by the points routine.
 
        call point4 (x2p,y2p)
        if (numpts .gt. 0)                                               &  
     &    call points4 (x1p,y1p,numpts,0,0,deltax,deltay)
 
! Reset the mapping to the users coordinates
 
        call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
 
        if (tgname .eq. 'LINEP') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
