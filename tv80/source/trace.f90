!#include "f77_dcomplx.h"
        subroutine trace (x,y,num,incx,incy,delx,dely)
 
!*************************************************************************
!
!  trace  -  Draw a series of lines
!
!  synopsis     call trace (x,y,num)
!               call trace (x,y,num,incx)
!               call trace (x,y,num,incx,incy)
!               call trace (x,y,num,incx,incy,  0.0)
!               call trace (x,y,num,incx,incy,  0.0, 0.0)
!               call trace (x,y,num,0   ,incy, delx)
!               call trace (x,y,num,0   ,0   , delx,dely)
!               call trace (x,y,num,incx,0   ,  0.0,dely)
!
!               real x,y                Endpoint of the lines
!               integer num             Number of endpoints
!               integer incx,incy       Increment between successive x
!                                       and y storage locations
!               real delx,dely          The amount to be added to the
!                                       first element in the array when
!                                       the increment is zero, in order
!                                       to generate succesive locations
!
!  description  Draws a series of vectors.  The default value for incx and
!               incy is 1.  If incx or incy is 0, then we generate the
!               succesive endpoint by adding delx or dely to the previous
!               endpoint.
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
      INTEGER incx,incy,num
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   delx,dely
      REAL   REAL
!============
        integer index,iaddx,iaddy
!        real deltx,delty,x(1),y(1)
!        real xpos,ypos
      REAL*8 deltx,delty,x(1),y(1)
      REAL*8 xpos,ypos
 
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
 
        if (tgname .eq. 'NONE') tgname = 'TRACE'
 
 
 
! Set the value of iaddx and deltx from arguments
 
        if (incx .eq. -1) then
        iaddx = 1
        deltx = 0.0 
      else if (incx .eq. 0) then
        iaddx = 0
        deltx = delx
      else
        iaddx = incx
        deltx = 0.0 
      endif
 
! Set the value of iaddy and delty from arguments
 
        if (incy .eq. -1) then
        iaddy = 1
        delty = 0.0 
      else if (incy .eq. 0) then
        iaddy = 0
        delty = dely
      else
        iaddy = incy
        delty = 0.0 
      endif
 
! If 3 argument version used, we  loop through the points in the array
 
        if (iaddx .eq. 1 .and. iaddy .eq. 1) then
 
          call setcrt (x(1),y(1))
          do 10 index = 2,num
            call vector (x(index),y(index))
10        continue
 
        else
 
! If more than 3 arguments were passed, life gets a bit more complicated, the
! user has control of the amount to add to the index and the difference to
! add to the points to generate the new points.
 
 
          call setcrt (x(1),y(1))
          do 20 index = 1,num-1
            xpos = x(1+index*iaddx)+deltx*REAL(index)
            ypos = y(1+index*iaddy)+delty*REAL(index)
            call vector (xpos,ypos)
20        continue
 
        endif
 
9999    if (tgname .eq. 'TRACE') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
