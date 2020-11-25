!#include "f77_dcomplx.h"
        subroutine tracec (ch,x,y,num,incx,incy,delx,dely)
 
!*************************************************************************
!
!  tracec  -  Draw a series of lines
!
!  synopsis     call tracec (ch,x,y,num)
!               call tracec (ch,x,y,num,incx)
!               call tracec (ch,x,y,num,incx,incy)
!               call tracec (ch,x,y,num,incx,incy,  0.0)
!               call tracec (ch,x,y,num,incx,incy,  0.0, 0.0)
!               call tracec (ch,x,y,num,0   ,incy, delx)
!               call tracec (ch,x,y,num,0   ,0   , delx,dely)
!               call tracec (ch,x,y,num,incx,0   ,  0.0,dely)
!
!               character ch            Character to plot at each point
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
!               endpoint.  The character ch will be place at each point
!               if it is more than ilnspace points farther away from the
!               last point where the character was plotted.
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
      INTEGER incx,incy,num,ierr
      INTEGER ioldhz,ioldvt,isqrsp
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   pxmin,pxmax,pymin,pymax,xmult,xadd,ymult,yadd,deltx
      REAL   delty
      REAL   REAL
!============
        character*1 ch
        integer index,iaddx,iaddy
!        real deltx,delty,x(1),y(1)
!        real xpos,ypos,xlin,ylin
      REAL*8 delx,dely,x(1),y(1)
      REAL*8 xpos,ypos
      REAL   xlin,ylin
        REAL   xlastc,ylastc,xnextc,ynextc,xdiffc,ydiffc
 
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
 
#ifndef NCAR_DUMMY 
 
! Set library entry routine name
 
        if (tgname .eq. 'NONE') tgname = 'TRACEC'
 
 
 
! allignment is centered in horizontal and vertical directions, character
! clipping is turned on if clipping was turned on.
 
        call gqtxal (ierr,ioldhz,ioldvt)
        call gstxal (2,3)
        if (iclip .eq. 1) call tgchclip (1)
        autofeed = .FALSE.
 
! calculate constants for conversion of linear user coordinates to 1024 x 1024
 
        pxmin = xvmin * 1024. 
        pxmax = xvmax * 1024. 
        pymin = yvmin * 1024. 
        pymax = yvmax * 1024. 
        xmult = (pxmax-pxmin)/(xwmax-xwmin)
        xadd = pxmin - xwmin * xmult
        ymult = (pymax-pymin)/(ywmax-ywmin)
        yadd = pymin - ywmin * ymult
        isqrsp = ilnspace*ilnspace
 
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
 
! plot the first point.  A character is plotted at this point if possible.
 
        call setcrt (x(1),y(1))
        call culxy (x(1),y(1),xlin,ylin)
        xlastc = xlin*xmult + xadd
        ylastc = ylin*ymult + yadd
        if (xlastc .ge. 0. .and. xlastc .lt. 1024. .and.                 &  
     &      ylastc .ge. 0. .and. ylastc .lt. 1024. ) then
          call setlch (x(1),y(1),ilncase,ilnindx,0,ilnfont)
          call gtext (ch,1,0)
        endif
 
        do 10 index = 1,num-1
          xpos = x(1+index*iaddx)+deltx*REAL(index)
          ypos = y(1+index*iaddy)+delty*REAL(index)
          call vector (xpos,ypos)
          call culxy (xpos,ypos,xlin,ylin)
          xnextc = xlin*xmult + xadd
          ynextc = ylin*ymult + yadd
          xdiffc = xnextc - xlastc
          ydiffc = ynextc - ylastc
          if ((xdiffc*xdiffc + ydiffc*ydiffc) .ge. isqrsp) then
            if (xnextc .ge. 0. .and. xnextc .lt. 1024. .and.             &  
     &          ynextc .ge. 0. .and. ynextc .lt. 1024. ) then
              call setlch (xpos,ypos,ilncase,ilnindx,0,ilnfont)
              call gtext (ch,1,0)
            endif
            xlastc = xnextc
            ylastc = ynextc
          endif
10      continue
 
! Reset the text horizontal, vertical alignment and character clipping.
 
        call gstxal (ioldhz,ioldvt)
        call tgchclip (0)
        autofeed = .TRUE.
 
        if (tgname .eq. 'TRACEC') tgname = 'NONE'
#endif

        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
