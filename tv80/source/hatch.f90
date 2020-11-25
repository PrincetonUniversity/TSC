!#include "f77_dcomplx.h"
      subroutine hatch( xvert, yvert, npoints, phi, ngrid, mode )
!
 
!**************************************************************************
!
!  hatch.F - area shading subroutine
!
!  contents    hatch - draws lines in specified area
!        hatchmap
!        hatchunm
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
 
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     H A T C H
!
!     Provide shading for a general polygonal region.  There is absolutely no
!     assumption made about convexity.  A polygon is specified by its vertices,
!     given in either a clockwise or counter-clockwise order.  The density of
!     the shading lines (or points) and the angle for the shading lines are
!     both determined by the parameters passed to the subroutine.
!
!     The calling sequence is:
!
!        call hatch( xvert, yvert, npoints, phi, ngrid, mode )
!
!     The input parameters are interpreted as follows:
!
!        xvert    -  An array or x coordinates for the polygon vertices
!
!        xvert    -  An array or x coordinates for the polygon vertices
!
!        npoints  -  The number of vertices in the polygon
!
!        phi      -  The angle for the shading, measured counter-clockwise
!                    in radians from the positive x-axis
!
!        ngrid    -  A parameter  determining the shading density:
!                      -n => every n-th shading line
!                        ....
!                      -2 => every other shading line
!                      -1 => every shading line
!                       0 => every shading line
!                      +1 => every other shading line
!                      +2 => every fourth shading line
!                        ....
!                      +n => every 2**n-th shading line
!
!        mode     -  A parameter determining the shading mode:
!                      -1 =>    no shading - boundary drawn
!                       0 =>  line-shading - boundary not drawn
!                       1 => point-shading - boundary not drawn
!                       2 =>  line-shading - boundary drawn
!                       3 => point-shading - boundary drawn
!
!     All coordinates are assumed to be in the user coordinate system as
!     specified in the most recent call to the tv80lib mapping subroutines.
!     Either Cartesian or polar coordinates are acceptable, although weird
!     (but aesthetic) results are produced using log-log or semi-log plots.
!     Direct tv80lib raster units (integer) may also be used.
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE tgcblk     
      USE tgchar     
      USE tgcolr     
      USE tgmap      
      USE tgpnts     
      USE hacomm
      USE module_hatch
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER mode,nvert,it,i
      INTEGER ivert,icount,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   vl,vr,vb,vt,wl
      REAL   wr,wb,wt,step,xmin,xmax,ymin,ymax,y,yhead,ytail,delx
      REAL   dely,delta,xkey,xtemp
      REAL   REAL
!============
      integer npoints, ngrid
!     real xvert(npoints), yvert(npoints), phi
      logical points,bound
 
!     common /hacomm/ sinphi,cosphi
 
!
!     This subroutine has to maintain an internal array of the transformed
!     coordinates.  This establishes a storage limitation.  The parameter
!     "nsize" is the maximum number of vertices allowed.  If the user wants
!     to increase the size of the arrays, make sure that the three common
!     blocks - hatchcom1, hatchcom2, hatchcom3 - are declared before this
!     routine and that the variable "ndimen" which is in common block
!     - hatchcom0 - is initialized at run time (not data-loaded) to the
!     correct dimension.  The maximum number of vertices is one less than
!     the dimension of the work arrays.
!
!      parameter ( nsize = 101 )
!
!     common / hatchcm0 / ndimen
!     common / hatchcm1 / cxvert(101)
!     common / hatchcm2 / cyvert(101)
!     common / hatchcm3 / xintercept(101)
 
! include the standard common block
 
      dimension xvert(*), yvert(*)
      REAL*8 xvert, yvert, phi, x1, x2, y1, y2
 
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
 
 
 
!
 
! Set library entry routine name
 
      if (tgname .eq. 'NONE') tgname = 'HATCH'
 
 
 
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     Check for valid number of vertices.
!
      if (npoints .lt. 2 .or. npoints .gt. ndimen-1) then
        goto 9999
      endif
!
!     See if point plotting was specified and whether to draw the
!     outline of the polygons.  Set flags for later.
!
      points = (mode .eq. 1) .or. (mode .eq. 3)
      bound  = (mode .eq. 2) .or. (mode .eq. 3) .or. (mode .eq. -1)
!
! If the last vertex is not the same as the first, set nvert to npoints+1
! to signify that we need to add a vertex.
!
      if (xvert(npoints) .ne. xvert(1) .or.                              &  
     &    yvert(npoints) .ne. yvert(1)) then
        nvert = npoints + 1
      else
        nvert = npoints
      endif
!
! Draw the border if required.  Return if no other processing is required.
!
      if (bound) then
        call trace (xvert,yvert,npoints,-1,-1,0.0 ,0.0 )
        if (nvert .ne. npoints)                                          &  
     &    call line (xvert(1),yvert(1),xvert(npoints),yvert(npoints))
      endif
      if (mode .eq. -1) goto 9999
!
! Save the old mapping and set values for rotation.
!
      call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
      sinphi = sin(phi)
      cosphi = cos(phi)
!
! Convert the users coordinate first into linear and then into
! normalized coordinates.  Various problems arise when the users
! input values have a very high or very low magnitude.  The easiest
! way to solve these problems is to convert the coordinates into
! normalized coordinats.
!
 
      do 5 i = 1,npoints
        call culxy (xvert(i),yvert(i),cxvert(i),cyvert(i))
        call clnxy (cxvert(i),cyvert(i),cxvert(i),cyvert(i))
5     continue
 
!
! Now that our input coordinate are all normalized, change the mapping
! into a normalized map.
!
 
      call tgsetmap (vl,vr,vb,vt,vl,vr,vb,vt,1)
 
!
! Compute the spacing from the parameters.  The finest spacing is 1024
! points across the screen.  This corresponds to either an input of
! ngrid = -1 or ngrid = 0.
!
      if (ngrid .ge. 0) then
        step = REAL(2 ** (min(ngrid,10))) / 1023. 
      else
        step = REAL(abs (ngrid)) / 1023. 
      endif
!
! Find the maximum and minimum x and y coordinates after rotation.  The
! user specifies an angle of rotation for the hatch lines.  It is easier
! to rotate the vertices and then compute horizontal lines in some sort
! of scanline algorithm.  The resulting 'scan' lines are then unrotated
! and drawn to produce the desired shading lines.  hatchmap will rotate
! the vertex. hatchunm will unrotate.
!
      call hatchmap (cxvert(1),cyvert(1))
      xmin = cxvert(1)
      xmax = xmin
      ymin = cyvert(1)
      ymax = ymin
      do 10 i = 2,npoints
        call hatchmap (cxvert(i),cyvert(i))
        xmax = max (cxvert(i),xmax)
        xmin = min (cxvert(i),xmin)
        ymax = max (cyvert(i),ymax)
        ymin = min (cyvert(i),ymin)
10    continue
      if (nvert .ne. npoints) then
        cxvert(nvert) = cxvert(1)
        cyvert(nvert) = cyvert(1)
      endif
!
! Do the scan line algorithm on the rotated vertices stored in cxvert
! and cyvert.
!
      do 20 y = ymin,ymax,step
        ivert = 0
        icount = 0
        do 30 i = 1,nvert-1
          yhead = y - cyvert(i+1)
          ytail = y - cyvert(i)
          if (sign(1. ,yhead) .ne. sign(1. ,ytail)) then
            icount = icount + 1
            delx = cxvert(i+1) - cxvert(i)
            dely = cyvert(i+1) - cyvert(i)
            delta = delx/dely * yhead
            xintercept (icount) = delta + cxvert(i+1)
          endif
30      continue
 
!
!       Sort the x intercept values.  Use a bubblesort because there aren't
!       very many of them (usually only two).
!
 
        do 40 i = 1,icount
          xkey = xintercept(i)
          do 50 j = 1,i-1
            if (xintercept(j) .gt. xkey) then
              xtemp = xkey
              xkey = xintercept(j)
              xintercept(j) = xtemp
            endif
50        continue
          xintercept(i) = xkey
40      continue
 
!
! All of the x coordinates for the shading segments along the current
! shading line are now known and are in sorted order.  All that remains
! is to draw them.  Process the x coordinates two at a time.
!
        do 60 i = 1, icount, 2
!
          x1 = xintercept(i)
          x2 = xintercept(i+1)
          y1 = y
          y2 = y
!
! Rotate back to original direction.
!
          call hatchunm (x1, y1)
          call hatchunm (x2, y2)
!
! See if plotting lines or points.
!
          if (points) then
            call linep (x1,y1, x2,y2, ngrid)
          else
            call line (x1, y1, x2, y2)
          endif
60      continue
20    continue
 
!
! Reset the map to what it was before hatch screwed it up. All of
! the 'goto 9999' statements in this subroutine are performed before
! the change to the mapping. So 9999 is places after the map reset.
!
 
      call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
 
9999  if (tgname .eq. 'HATCH') tgname = 'NONE'
      return
 
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
