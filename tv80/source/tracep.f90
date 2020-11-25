        subroutine tracep (xarray,yarray,num,k,incx,incy)
 
!*************************************************************************
!
!  tracep  -  Draw a series of dotted lines.
!
!  synopsis     call tracep (xarray,yarray,num)
!               call tracep (xarray,yarray,num,k)
!               call tracep (xarray,yarray,num,k,incx)
!               call tracep (xarray,yarray,num,k,incx,incy)
!
!               real xarray(num)        X coordinates
!               real yarray(num)        Y coordinates
!               integer num             Number of coordinates
!               integer k               Exponent to pass to linep
!               integer incx            Increment between x coordinates
!               integer incy            Increment between y coordinates
!
!  description  Plots a series of points between the points given in
!               xarray and yarray.  The points are drawn by making a
!               call to linep.
!
!*************************************************************************
 
!       real xarray(1)
!        real yarray(1)
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE tgcblk     
      USE tgchar     
      USE tgcolr     
      USE tgmap      
      USE tgpnts     
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ixindx,iyindx,i
!============
      REAL*8 xarray(1)
      REAL*8 yarray(1)
        integer num
        integer k
        integer incx
        integer incy
!
! The value of k needs to be stored between successive calls to tracep.
! The initial value of k is 2.
!
        integer ix,iy
 
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
 
        if (tgname .eq. 'NONE') tgname = 'TRACEP'
 
 
 
        ix = 1
        iy = 1
 
        if (incx .ne. -1) ix = incx
      if (incy .ne. -1) iy = incy
      if (k .ne. -1) iptspace = k
 
        ixindx = 1
        iyindx = 1
        do 10 i = 1,num-1
          call linep (xarray(ixindx),yarray(iyindx),                     &  
     &      xarray(ixindx+ix),yarray(iyindx+iy),iptspace)
          ixindx = ixindx + ix
          iyindx = iyindx + iy
10      continue
 
        if (tgname .eq. 'TRACEP') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
