        subroutine cntxy (x,y,xnew,ynew)
 
!**********************************************************************
!
!  tgtran.F  -  transformation routines
!
!  contents    cntxy  -  convert normalized point to a transformed point
!        init2d  -  Initialize transformation matrix
!        tran2d  -  Set translation
!        scal2d  -  Set scale
!        rot2d  -  Set rotation
!        center  -  Set the center of scalings and rotations
!        catmat  -  Concatenate matrices
!        makidmat  -  Make an identity matrix
!        multmat  -  Multiply the transformation matrix by another
!           matrix
!
!  description  All routines in tglib which do transformations will
!               need to be fed through the transformation routines
!               to obtain the correct coordinates.
!
!    National Energy Research Supercomputer Center
!    Lawrence Livermore National Laboratory
!    University of California
!    Livermore, California 94550, USA
!
!    (C) Copyright 1991 The Regents of the University of California.
!    All Rights Reserved.
!
!**********************************************************************
 
!**********************************************************************
!
!  cntxy  -  convert normalized point to a transformed point
!
!  synopsis     call cntxy (x,y,xnew,ynew)
!               real xnew,ynew       transformed x,y
!               real x,y             Point to transform
!
!**********************************************************************
 
 
! The temporary variables have two purposes.  Obviously they save the
! converted value of x and y so that the conversion routines don't
! have to be called twice.  By saving the values in temporary variables
! the input x,y can be the same as the output x,y.  For example:
!   call cntxy (x,y,x,y)
 
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
! idecl:  explicitize implicit REAL declarations:
      REAL   y,xnew,ynew,x
!============
        REAL   xtemp,ytemp
 
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
 
 
 
        if (.not. dotrans) then
          xnew = x
          ynew = y
          return
        endif
 
        if (.not. matmade) call catmat
 
        xtemp = x
        ytemp = y
        xnew = xtemp*trnmat(1,1) + ytemp*trnmat(2,1) + trnmat(3,1)
        ynew = xtemp*trnmat(1,2) + ytemp*trnmat(2,2) + trnmat(3,2)
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
