        subroutine catmat
 
!**********************************************************************
!
!  catmat  -  Concatenate matrices
!
!  synopsis     call catmat ()
!
!  description  Concatenates the matrices to create one transformation
!               matrix.  The effect will be to move the center to the
!               origin, scale, rotate, move back out to the original
!               position, and then move the amount specified in the
!               translation.
!
!**********************************************************************
 
!=======================================================================
! 12/28/2001 ler pppl - Change first two argments in 'call culxy' from
!                        single precision to double precision
!=======================================================================
 
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
        REAL   trxnew,trynew
        REAL   ctxnew,ctynew
 
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
 
 
 
        REAL   matrix (3,3)
        REAL*8 dcentrx, dcentry
 
 
! convert translation and center into the proper units
 
        call clnxy (xwmin+transx,ywmin+transy,trxnew,trynew)
        trxnew = trxnew - xvmin
        trynew = trynew - yvmin
 
        dcentrx = centrx
        dcentry = centry
        call culxy (dcentrx,dcentry,ctxnew,ctynew)
        call clnxy (ctxnew,ctynew,ctxnew,ctynew)
 
! erase old trasformation matrix
 
        call makidmat (trnmat)
 
! move to the origin
 
        call makidmat (matrix)
        matrix (3,1) = -ctxnew
        matrix (3,2) = -ctynew
        call multmat (matrix)
 
! scale
 
        call makidmat (matrix)
        matrix (1,1) = scalex
        matrix (2,2) = scaley
        call multmat (matrix)
 
! rotate
 
        call makidmat (matrix)
        matrix (1,1) = cos (rotat)
        matrix (1,2) = sin (rotat)
        matrix (2,1) = - matrix(1,2)
        matrix (2,2) = matrix(1,1)
        call multmat (matrix)
 
! move center back where it was, add in the translation
 
        call makidmat (matrix)
        matrix (3,1) = ctxnew + trxnew
        matrix (3,2) = ctynew + trynew
        call multmat (matrix)
 
        matmade = .TRUE.
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
