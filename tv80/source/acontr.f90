      subroutine acontr (jjj1, bc, jjj2, ba, max, imin, imax, istop      &  
     & , jmin, jmax, jstop)
 
!**************************************************************************
!
!  contur.F - contouring routines
!
!  contents    acontr
!        rcontr
!        contur
!        xcontr
!        lcurvs
!        q7plot
!        q8plot
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
 
!
!....added autofeed to tgchar common in gaxisi
!....added block data statement to define ndimen
!...changed inta to int on lines 3168 and 3172
 
!**************************************************************************
!
!  acontr
!
!  calls        xcontr
!
!**************************************************************************
 
!=======================================================================
! 12/28/2001 ler pppl - change last two arguments in call to xcontr to
!                        double precision
!=======================================================================
!
!c variable declarations:
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
      INTEGER jjj2,imin,imax,istop,jmin,jmax,jstop,jjj1
      INTEGER i
      INTEGER max
!============
      dimension bc(*), ba(*)
      REAL*8 bc, ba
!     dimension ba(1), bc(1)
 
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
 
 
 
! If called by user set routine top routine name
 
      if (tgname .eq. 'NONE') tgname = 'ACONTR'
 
 
 
!
!c program statements:
!
!% The following call to xcontur used to only have 12 arguments.  The
!% last 2 arguments were set to -1 to make the call consistent with
!% the definition.  The last 2 arguments should not be used by xcontr.
!
      i = 1
      call xcontr (i, jjj1, bc(1), jjj2, ba, max, imin, imax, istop      &  
     & , jmin, jmax, jstop, -1.0d0 , -1.0d0 )
!
9999  if (tgname .eq. 'ACONTR') tgname = 'NONE'
      return
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
