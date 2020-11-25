      subroutine plotv (xx1, yy1, xx2, yy2, fact)
 
!*************************************************************************
!
!  plotv  -  Draw an arrow
!
!  synopsis     call plotv (x1,y1,x2,y2)
!               call plotv (x1,y1,x2,y2,fact)
!               real x1,y1      Start of arrow
!               real x2,y2      Head of arrow
!               real fact       Some sort of factor to indicate head size
!
!  description  Creates an arrow pointing from x1,y1 to x2,y2.  fact
!               specifies how large the arrowhead will be.  fact should
!               be between 0. and 1.  This routine was extracted from
!               tv80lib and rewritten.
!
!*************************************************************************
 
!=======================================================================
! 12/28/2001 ler pppl - change 'call line' and 'call vector' arguments
!                        from single precision to double precision
!=======================================================================
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE tgcblk     
      USE tgchar     
      USE tgcolr     
      USE tgmap      
      USE tgpnts     
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER it
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   con,gx1,gy1,gx2,gy2,vl,vr,vb,vt,wl,wr,wb,wt,dx,dy,pro
      REAL   ax,by,ay,bx
!============
      REAL*8 xx1, yy1, xx2, yy2, fact, dgx1, dgx2, dgy1, dgy2
!
!c purpose: plot a vector (pat crowleys method - may, 1968)
!           (argument h, if given, is fixed length of arrowhead -
!            otherwise range)
!           (length of arrowhead proportional to length of v)
!
!  Default legth of arrowhead
 
        parameter (con=0.12 )
 
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
 
        if (tgname .eq. 'NONE') tgname = 'PLOTV'
 
 
 
! Convert the coordinates into normalized coordinates.  Then convert the map
! into a normalized linear map.
 
        call culxy (xx1,yy1,gx1,gy1)
        call clnxy (gx1,gy1,gx1,gy1)
        call culxy (xx2,yy2,gx2,gy2)
        call clnxy (gx2,gy2,gx2,gy2)
        call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
        call tgsetmap (vl,vr,vb,vt,vl,vr,vb,vt,1)
!
        dx = gx2-gx1
        dy = gy2-gy1
        pro = con
        if (fact .ne. -1.0 ) then
          pro = dx*dx+dy*dy
          if (pro .gt. 0) pro = fact/sqrt(pro+pro)
        endif
!
        dgx1 = gx1
        dgy1 = gy1
        dgx2 = gx2
        dgy2 = gy2
        call line (dgx1, dgy1, dgx2, dgy2)
        ax = (dy-dx)*pro
        by = -(dy-dx)*pro
        ay = (-dy-dx)*pro
        bx = (-dy-dx)*pro
        call line (dgx2+ax, dgy2+ay, dgx2, dgy2)
        call vector (dgx2+bx, dgy2+by)
 
! Reset the map to the original map
 
        call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
 
        if (tgname .eq. 'PLOTV') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
