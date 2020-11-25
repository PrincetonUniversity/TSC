!#include "f77_dcomplx.h"
!#define TRANSP 1
#define TRANSP 
#ifdef TRANSP
      SUBROUTINE EZinit
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      return
      END
      SUBROUTINE EZfini(FinFlag1,FinFlag2)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER FinFlag1,FinFlag2
      return
      END
      SUBROUTINE EZdra ( x, y, i )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 x, y
      INTEGER i
      return
      END
      SUBROUTINE EZcros( x, y, n )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8  x, y
      INTEGER n
      return
      END
      SUBROUTINE EZsets(i,j,k,l,a,b,c,d,kind)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 a,b,c,d
      INTEGER i,j,k,l,kind
      return
      END
      SUBROUTINE EZaxes ( incx , incxi , incy , incyi )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER incx,incxi,incy,incyi
      return
      END
      SUBROUTINE EZscal ( incx , incxi , incy , incyi )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER incx,incxi,incy,incyi
      return
      END
      SUBROUTINE EZsymb(x,y,n,ascii)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*1 ascii
      INTEGER i, n
      REAL*8    x(n), y(n)
      return
      END
      SUBROUTINE EZwrit(i,j,string,lcen,lor)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*(*) string
      INTEGER i,j,lcen,lor
      return
      END
      SUBROUTINE EZcurv (x,y,n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 x,y
      INTEGER n
      return
      END
      SUBROUTINE EZpnts(x,y,n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 x,y
      INTEGER n
      return
      END
      SUBROUTINE EZpnt2(x,y,n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 x,y
      INTEGER n
      return
      END
      SUBROUTINE EZbars(x,y,n,ch)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*1 ch
      REAL*8 x,y
      INTEGER n
      return
      END
      SUBROUTINE EZeras (i,ii)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,ii
      return
      END
      SUBROUTINE SysSetUp
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      return
      END
      SUBROUTINE SysClose
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      return
      END
      SUBROUTINE EZrCon(i,j,k,l,a,b,c,d,                                 &  
     &  kclev1,clevelin,kclev2,                                          &  
     &  AryContr,                                                        &  
     &  i1stdim,                                                         &  
     &  xAry, ixmin, ixmax, ixstep,                                      &  
     &  yAry, jymin, jymax, jystep )
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,j,k,l,kclev1,kclev2,i1stdim,ixmin,ixmax,ixstep
      INTEGER jymin,jymax,jystep
!============
      REAL*8    a,b,c,d
      REAL*8    clevelin, AryContr(i1stdim,jymax)
      REAL*8    xAry(ixmax), yAry(jymax)
      return
      END
      SUBROUTINE EZcurvmd (x,x1stdim,y,y1stdim,n)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n
!============
      INTEGER x1stdim, y1stdim
      REAL*8    x(x1stdim,n), y(y1stdim,n)
      return
      END
!
!
#else
!
!
!     EZ SG     begins                      ---------------------------|
!                                                                      |
!                                                                      |
!     EZ ---> SGlib by D. W. Ignat and F. C. Jobes,
!     Princeton Plasma Physics Laboratory
!     Princeton, NJ 08543
!     1990, 1991, 1992, ... , 2000
!     General idea borrowed from Dori Miller Barnes, PPPL
!
!----------------------------------------------------------------------
      SUBROUTINE EZinit
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL binitt, initt
      include "sgIbk.inc"
!     initt(0) does not call ERASE
!     initt(>0) does call ERASE
      call initt(240)
      call binitt
      level = 1
      return
      END
!----------------------------------------------------------------------
      SUBROUTINE EZfini(FinFlag1,FinFlag2)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL upause, finitt, tsend
!     upause(0) beep, wait for input, do not erase
!     upause(1) beep, wait for input, erase
!     finitt(i,j) puts cursor at i,j
!     FinFlag 1,2 not used in ezsg, but may be used in ez80 for cray
      INTEGER FinFlag1, FinFlag2
      include "sgIbk.inc"
      call tsend
      if ( level .eq. 3 ) then
           call upause(1)
           level=1
      else
           call finitt(0,0)
           level=0
      endif
      return
      END
!----------------------------------------------------------------------
      SUBROUTINE EZdra ( x, y, i )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL movabs, drwabs
      INTEGER i
      REAL*8    x, y
      include "sgIbk.inc"
      include "sgRbk.inc"
      INTEGER ixplot, iyplot
 
      if( itype .eq. 1 .or. itype .eq. 3 ) then
        ixplot = int( (x-x1)*xpix  ) + ix1
      else
        ixplot = int( log10(x/x1)*xpix ) + ix1
      endif
 
      if ( itype .eq. 1 .or. itype .eq. 2 ) then
        iyplot = int( (y-y1)*ypiy  ) + iy1
      else
        iyplot = int( log10(y/y1)*ypiy ) + iy1
      endif
      if (ixplot .le. ix2 .and. ixplot .ge. ix1 .and.                    &  
     &    iyplot .le. iy2 .and. iyplot .ge. iy1      ) then
 
            if(i.eq.0) call movabs(ixplot,iyplot-1)
            if(i.eq.1) call drwabs(ixplot,iyplot)
      else
            ixplot = min(ixplot,ix2)
            ixplot = max(ixplot,ix1)
            iyplot = min(iyplot,iy2)
            iyplot = max(iyplot,iy1)
            if(i.eq.0) call movabs(ixplot,iyplot)
            if(i.eq.1) call drwabs(ixplot,iyplot)
      endif
      return
      END
!----------------------------------------------------------------------
      SUBROUTINE EZcros ( xs, ys, n )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL movabs, drwabs
      INTEGER i, n
      REAL*8    xs(n), ys(n), x, y
      include "sgIbk.inc"
      include "sgRbk.inc"
      INTEGER ixplot, iyplot, size
      DATA    size / 1 /
 
      do 10 i=1,n
      x = xs(i)
      y = ys(i)
      if( itype .eq. 1 .or. itype .eq. 3 ) then
        ixplot = int( (x-x1)*xpix  ) + ix1
      else
        ixplot = int( log10(x/x1)*xpix ) + ix1
      endif
 
      if ( itype .eq. 1 .or. itype .eq. 2 ) then
        iyplot = int( (y-y1)*ypiy  ) + iy1
      else
        iyplot = int( log10(y/y1)*ypiy ) + iy1
      endif
 
      if (ixplot .le. ix2 .and. ixplot .ge. ix1 .and.                    &  
     &    iyplot .le. iy2 .and. iyplot .ge. iy1      ) then
 
            call movabs(ixplot-size,iyplot-size)
            call drwabs(ixplot+size,iyplot+size)
            call movabs(ixplot+size,iyplot-size)
            call drwabs(ixplot-size,iyplot+size)
      endif
 10   continue
 
      call movabs(ixplot     ,iyplot     )
 
      return
      END
!----------------------------------------------------------------------
      SUBROUTINE EZsets(i,j,k,l,a,b,c,d,kind)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL slimx, slimy, dlimx, dlimy, setwin, xtype, ytype
      include "sgIbk.inc"
      include "sgRbk.inc"
      INTEGER i,j,k,l,kind
      REAL*8    a,b,c,d
 
        call xneat(0)           !frame data exactly
        call yneat(0)           !try 15Jun94 to keep limits chosen
 
      call slimx(i,j)
      call slimy(k,l)
      call dlimx(a,b)
      call dlimy(c,d)
      call setwin
      level = 2
 
      if     ( kind .eq. 1)  then
         call xtype(1)
         call ytype(1)
         xpix = AREAL(j-i)/(b-a)
         ypiy = AREAL(l-k)/(d-c)
      else if ( kind .eq. 2) then
         call xtype(2)
         call ytype(1)
         xpix = AREAL(j-i)/log10(b/a)
         ypiy = AREAL(l-k)/(d-c)
      else if ( kind .eq. 3) then
         call xtype(1)
         call ytype(2)
         xpix = AREAL(j-i)/(b-a)
         ypiy = AREAL(l-k)/log10(d/c)
      else if ( kind .eq. 4) then
         call xtype(2)
         call ytype(2)
         xpix = AREAL(j-i)/log10(b/a)
         ypiy = AREAL(l-k)/log10(d/c)
      else
         call LSCstop( ' kind .ne. 1,2,3,4 in ezsets ')
      endif
 
      x1=a
      x2=b
      y1=c
      y2=d
      ix1=i
      ix2=j
      iy1=k
      iy2=l
      itype=kind
      level=3
      return
      END
!----------------------------------------------------------------------
      SUBROUTINE EZrCon (i,j,k,l,a,b,c,d,                                &  
     &  kclev1,clevelin,kclev2,                                          &  
     &  AryContr,                                                        &  
     &  i1stdim,                                                         &  
     &  xAry, ixmin, ixmax, ixstep,                                      &  
     &  yAry, jymin, jymax, jystep )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      include "sgIbk.inc"
      include "sgRbk.inc"
 
      INTEGER i,j,k,l,kind,                                              &  
     &        kclev1,kclev2,i1stdim,                                     &  
     &        ixmin,ixmax,ixstep,                                        &  
     &        jymin,jymax,jystep
      REAL*8    a,b,c,d
      REAL*8    clevelin, AryContr(i1stdim,jymax)
      REAL*8    xAry(ixmax), yAry(jymax)
      kind = 1
!
!     call EZsets  (i,j,k,l,a,b,c,d,kind)
      call agraphdi(i,j,k,l,a,b,c,d, ' ',' ')
!
      level=2
      x1=a
      x2=b
      y1=c
      y2=d
      ix1=i
      ix2=j
      iy1=k
      iy2=l
 
      call rcontr (                                                      &  
     &  kclev1,clevelin,kclev2,                                          &  
     &  AryContr,                                                        &  
     &  i1stdim,                                                         &  
     &  xAry, ixmin, ixmax, ixstep,                                      &  
     &  yAry, jymin, jymax, jystep )
 
      level=3
      return
 
      END
!
!----------------------------------------------------------------------
!
!  draw graph and map virtual space for contour plot
!  modified Mar 93 by di to try to keep limits from being changed
!  as advised by MLT
!
      SUBROUTINE agraphdi(minx,maxx,miny,maxy,xmin,xmax,ymin,ymax,       &  
     &	  xlab,ylab)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
        CHARACTER xlab,  ylab
!       CHARACTER xlab,  ylab ! added Jan 2000 to keep g77 happy; dont know if this is exactly correct
	DIMENSION xlab(1), ylab(1)
        INTEGER   minx,maxx,miny,maxy
        REAL*8      xmin,xmax,ymin,ymax
        REAL*8      x(2), y(2)
 
        call binitt
!       call xfrm(2)            ! short ticks
	call xfrm(2)
	call yfrm(2)
        call xneat(0)           !frame data exactly
        call yneat(0)
	call slimx(minx,maxx)
	call slimy(miny,maxy)
	call dlimx(xmin,xmax)
	call dlimy(ymin,ymax)
	call npts(2)
	x(1)=xmin
	x(2)=xmax
	y(1)=ymin
	y(2)=ymax
	call check(x,y)
	call agdspl(x,y)        ! draw left/bottom axes, label ticks
	call agtchs(xlab,0)     ! put user labels
	call agtcvs(ylab,0)
	call xloctp(0)
	call ylocrt(0)
	call grid
	return
	end
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZaxes ( incx , incxi , incy , incyi )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL movabs, drwabs, xneat, yneat, xtics, ytics,               &  
     &                         xmtcs, ymtcs, xfrm , yfrm ,               &  
     &                         xmfrm, ymfrm, xlab , ylab ,               &  
     &           npts,  check, agdspl
      INTEGER             incx , incxi , incy , incyi
      REAL*8 x(2),y(2)
      include "sgIbk.inc"
      include "sgRbk.inc"
!     Draw a frame, right side and top.
      call movabs(ix2,iy1)
      call drwabs(ix2,iy2)
      call drwabs(ix1,iy2)
!
      call xneat (2)
!     0 endpoints at raw data limits
!     1 endpoints at next minor tick
!     2 endpoints at next major tick
      call yneat (2)
      call xtics ( incx )
      call ytics ( incy )
      call xmtcs ( incxi)
      call ymtcs ( incyi)
      call xfrm(3)
      call yfrm(3)
      call xmfrm(2)
      call ymfrm(2)
 
      call xlab  ( 1 )
!     0 no label
!     1 label major ticks
!     2 label end ticks
      call ylab  ( 2 )
 
      x(1)=x1
      x(2)=x2
      y(1)=y1
      y(2)=y2
       call npts(2)
       call check(x,y)
      call agdspl(x,y)
      return
      END
 
!        set form of major ticks, default=5, 0=no axis
!all xfrm(3)            ! 1=no ticks,   2=short out,   3=short inside
!all yfrm(3)            ! 4=short thru, 5=grid thru,   6=grid inside
!         set form of minor ticks, default=2
!all xmfrm(2)
!all ymfrm(2)
!         set tick label frequency, default=1
!all xlab(2)            ! 0=no label,   1=major ticks, 2=end ticks on
!all ylab(2)
!         set tick length, default=height of 1 character
!all xlen(8)            ! length in screen units
!all ylen(8)            ! minor ticks are half as long
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZscal ( incx , incxi , incy , incyi )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL movabs, drwabs, xneat, yneat, xtics, ytics,               &  
     &                         xmtcs, ymtcs, xfrm , yfrm ,               &  
     &                         xmfrm, ymfrm, xlab , ylab ,               &  
     &           npts,  check, agdspl
      INTEGER             incx , incxi , incy , incyi
      REAL*8 x(2),y(2)
      include "sgIbk.inc"
      include "sgRbk.inc"
!
!     Like ezaxes, but do not draw the frame around the graph
!
      call xneat (2)
!     0 endpoints at raw data limits
!     1 endpoints at next minor tick
!     2 endpoints at next major tick
      call yneat (2)
      call xtics ( incx )
      call ytics ( incy )
      call xmtcs ( incxi)
      call ymtcs ( incyi)
      call xfrm(3)
      call yfrm(3)
      call xmfrm(2)
      call ymfrm(2)
 
      call xlab  ( 1 )
!     0 no label
!     1 label major ticks
!     2 label end ticks
      call ylab  ( 2 )
 
      x(1)=x1
      x(2)=x2
      y(1)=y1
      y(2)=y2
       call npts(2)
       call check(x,y)
      call agdspl(x,y)
      return
      END
 
!        set form of major ticks, default=5, 0=no axis
!all xfrm(3)            ! 1=no ticks,   2=short out,   3=short inside
!all yfrm(3)            ! 4=short thru, 5=grid thru,   6=grid inside
!         set form of minor ticks, default=2
!all xmfrm(2)
!all ymfrm(2)
!         set tick label frequency, default=1
!all xlab(2)            ! 0=no label,   1=major ticks, 2=end ticks on
!all ylab(2)
!         set tick length, default=height of 1 character
!all xlen(8)            ! length in screen units
!all ylen(8)            ! minor ticks are half as long
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZsymb(x,y,n,ascii)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ezdra, agths
!     Probably needs some work adjusting the location to get the
!     character centered at the desired place.  Later, man.
      CHARACTER*1 ascii
      INTEGER i, n
      REAL*8    x(n), y(n)
      do 10 i=1,n
           call EZdra ( x(i), y(i), 0 )
           call agths(ascii,1)
 10   continue
      return
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZwrit(i,j,string,lcen,lor)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL agtmhs, agtmvs
      CHARACTER string*(*)
      CHARACTER*1 dollar
      INTEGER     i,j,lcen,lor, indx, kx, ky, MAX, n
      REAL*8        CHARX, CHARY
      DATA        dollar,   MAX,    CHARX,  CHARY  /                     &  
     &            '$'   ,   75 ,    7.0_R8,  11.0_R8/
      kx=i
      ky=j
      n=0
 
      do 10 indx=1,MAX
 
        if(string(indx:indx) .eq. dollar) go to 30
        n=n+1
10    continue
 
30    continue
      if(n.eq.0) return
      if(lcen.eq.0) go to 100
 
      kx= i - int( AREAL(n)*CHARX +0.5_R8)
      ky= j + int( AREAL(n)*CHARY +0.5_R8)
 
      if(lor.eq.0) ky=j
      if(lor.ne.0) kx=i
 
100   continue
 
      if(lor.eq.0) call agtmhs(kx,ky,string,n)
      if(lor.ne.0) call agtmvs(kx,ky,string,n)
      return
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZcurv (x,y,n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ezdra
      INTEGER i, n
      REAL*8    x(n), y(n)
      call EZdra ( x(1) , y(1) , 0 )
      do 10 i=1,n
      call EZdra ( x(i) , y(i) , 1 )
10    continue
      return
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZcurvmd (x,x1stdim,y,y1stdim,n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ezdra
      INTEGER i, n
      INTEGER x1stdim, y1stdim
      REAL*8    x(x1stdim,n), y(y1stdim,n)
        call EZdra ( x(1, 1) , y(1, 1) , 0 )
      do 10 i=1,n
        call EZdra ( x(1, i) , y(1, i) , 1 )
10    continue
      return
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZpnts(x,y,n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ezdra
      INTEGER i, n
      REAL*8    x(n), y(n)
!     You can plot a dot by calling pntabs(ix,iy) or pointa(x,y),
!     which does a move followed by a draw to the same location internally.
      do 10 i=1,n
      call EZdra (x(i), y(i), 0)
      call EZdra (x(i), y(i), 1)
!     call pointa(x(i), y(i))
 10   continue
      return
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZpnt2(x,y,n)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ezdra
      INTEGER i, n
      REAL*8    x(n), y(n)
      do 10 i=1,n,    2
      call EZdra (x(i), y(i), 0)
      call EZdra (x(i), y(i), 1)
!     call pointa(x(i), y(i))
 10   continue
      return
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZbars(x,y,n,ch)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ezdra
      CHARACTER*1 ch, chx, chy
      INTEGER i, n
      REAL*8    x(n), y(n), zero
      DATA    zero, chx, chy / 0.00_R8, 'x', 'y' /
      if      (ch .eq. chy) then
        do 10 i=1,n
        call EZdra (x(i), zero, 0)
        call EZdra (x(i), y(i), 1)
 10     continue
      else if (ch .eq. chx) then
        do 20 i=1,n
        call EZdra (zero, y(i), 0)
        call EZdra (x(i), y(i), 1)
 20     continue
      else
        do 30 i=1,n
        call EZdra (x(i), zero, 0)
        call EZdra (x(i), y(i), 1)
 30   continue
      endif
      return
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE EZeras (i,ii)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL initt
      INTEGER i, ii
      call initt (240)
      return
      END
!
!     ------------------------------------------------------------------
      SUBROUTINE SysSetUp
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      return
      END
      SUBROUTINE SysClose
      USE sgIbk
      USE sgRbk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      return
      END
!                                                                      |
!                                                                      |
!     EZ SG   ends                          ---------------------------|
!
!
#endif
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
