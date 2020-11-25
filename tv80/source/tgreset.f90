        subroutine tgreset ()
 
!*************************************************************************
!
!  tgreset  -  Reinitialize the interface
!
!  synopsis     call tgreset ()
!
!  description  This routine allows the user to muck with gks and then
!        call this routine to reset gks to what the library thinks
!        it should be.
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
      INTEGER iopstate
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
 
 
 
! Set library entry routine name
 
        if (tgname .eq. 'NONE') tgname = 'TGRESET'
 
 
 
! Make sure that gks was opened (operation state .gt. 0)
 
        call gqops (iopstate)
        if (iopstate .eq. 0) then
          call tgerror (20)
          call exit (1)
        endif
 
! Make sure the user has called tginit to initialize the library
 
        if (doinit) then
          call tgerror (21)
          call exit (1)
        endif
 
! Set the viewport, window and current transform
 
        call gsvp (1,xvmin,xvmax,yvmin,yvmax)
        call gswn (1,xvmin,xvmax,yvmin,yvmax)
        call gselnt (1)
 
! Set the clipping
 
        call gsclip (iclip)
 
! Set font number, character precision to stroke, up vector, character height,
! and character expansion factor.
 
        call gstxfp (ichfont,2)
        call gschup (chupx,chupy)
        call gschh (chparm(1,ichindx+1))
        call gschxp (chparm(3,ichindx+1))
 
! Set text path to right and allignment to left bottom
 
        call gstxp (0)
        call gstxal (1,5)
 
! Set color
 
        call gsplci (icurclr)
        call gspmci (icurclr)
        call gstxci (icurclr)
        call gsfaci (icurclr)
 
        if (tgname .eq. 'TGRESET') tgname = 'NONE'
        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
