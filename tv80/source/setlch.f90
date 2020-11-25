        subroutine setlch (x,y,icase,isize,iorient,ifont)
 
!*************************************************************************
!
!  setlch  - set position, intensity and size of characters
!
!  synopsis     call setlch (x,y)
!               call setlch (x,y,icase)
!               call setlch (x,y,icase,isize)
!               call setlch (x,y,icase,isize,iorient)
!               call setlch (x,y,icase,isize,iorient,ifont)
!
!               real x,y        Position
!               integer icase   Case for future calls
!               integer isize   Size of characters
!               integer iorient Orientation of strings
!               integer itype   Font to use
!
!  description  x,y will be the position in user coordinates for the
!               next string to be plotted.
!
!               intens is the intensity of the characters.
!                 0 - (default) low
!                 1 - high
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
      INTEGER isize,iorient,ifont,icase
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   xtemp,ytemp,xmax,xmin,ymax,ymin
!============
      REAL*8 x, y
 
        integer itsize,itorient,itfont
        REAL   x2,x3,x4,y2,y3,y4
 
        logical tgisfont
 
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
 
        if (tgname .eq. 'NONE') tgname = 'SETLCH'
 
 
 
! Set temporary variables to the old values of the character attributes
 
        itsize = ichindx
        itorient = ichangle
        itfont = ichfont
 
! If a font was passed, check to see if the font number passed to this
! routine is valid.  If it is not, then itfont will not be changed.
 
        if (ifont .ne. -1) then
          if (tgisfont (ifont)) itfont = ifont
        endif
 
! angle to rotate
 
        if (iorient .ne. -1) itorient = iorient
 
! specify sizes from 0 to 3
 
        if (isize .ne. -1) itsize = isize
 
! specify which case to use
 
        if (icase .ne. -1) ichcase = icase
 
! convert x,y to window linear coordinates and then to NDC coordinates
 
10      call culxy (x,y,chx,chy)
        call clnxy (chx,chy,chx,chy)
        call tgchset (itfont,itsize,itorient)
 
! If autofeed is not on, then don't reposition character. Go to end.
 
        if (.not. autofeed) goto 9999
 
! Calculate the approximate corners of the rectangle that bounds the initial
! character of the string in clockwise order.  tgchset will have set chupx to
! the cosine of the angle of the up vector and chupy to the sine of that angle.
! The up vector is at an angle of chrot+PI/2.
 
        xtemp = cos(chrot)
        ytemp = sin(chrot)
        x2 = (chparm(2,itsize+1)*chupx) + chx
        y2 = (chparm(2,itsize+1)*chupy) + chy
        x3 = (chparm(4,itsize+1)*xtemp) + x2
        y3 = (chparm(4,itsize+1)*ytemp) + y2
        x4 = (chparm(4,itsize+1)*xtemp) + chx
        y4 = (chparm(4,itsize+1)*ytemp) + chy
 
! If the rectangle is outside the NDC area, move it into this area.
 
        xmax = max (chx,max(x2,max(x3,x4)))
        xmin = min (chx,min(x2,min(x3,x4)))
        ymax = max (chy,max(y2,max(y3,y4)))
        ymin = min (chy,min(y2,min(y3,y4)))
 
        if (xmin .lt. 0. ) chx = 0. 
        if (xmax .gt. 1. ) chx = 1. - (xmax-xmin)
        if (ymin .lt. 0. ) chy = 0. 
        if (ymax .gt. 1. ) chy = 1. - (ymax-ymin)
 
9999    if (tgname .eq. 'SETLCH') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
