!#include "f77_dcomplx.h"
        subroutine tginit (iwkid)
 
!*************************************************************************
!
!  tginit  -  Initialize the interface
!
!  synopsis     call tginit (iwkid)
!        integer iwkid     workstation identifier
!
!  description  Initialize the tv80 to gks shell so that gks will be in
!               a state similar to the default state of tv80lib.  This
!               must be called after a every workstation is opened.
!        iwkid is the workstation identifier associated with the
!        workstation.
!
!*************************************************************************
 
 
!=======================================================================
! 12/28/2001 pppl ler - change 1 to variable i_one in call to gqdsp
!                       ... this change done earlier than 12/28/2001
!=======================================================================
! The tgdefs BLOCK DATA subprogram must be mentioned in a called routine
! in order for it to be included into the program.  tginit must be
! called for the shell to work properly, so it is placed here.
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE tgcblk     
      USE tgchar     
      USE tgcolr     
      USE tgmap      
      USE tgpnts     
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iwkid,ierr,itype,i_one
      INTEGER lx,ly,numclr,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   dummy,rx,ry,diff
      REAL   REAL
!============
        external tgdefs
 
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
 
 
 
! include the site dependent info
 
!***************************************************************************
!
!  tgsite.inc - site dependent information
!
!***************************************************************************
 
! The minimum and maximum font numbers
 
!     integer minfnt,maxfnt
!     parameter (minfnt = 1, maxfnt = 20)
 
! The default font
 
      integer deffnt
      parameter (deffnt = 1)
 
! The character height and expansion factor for fonts 1, 2, 3 and 4
 
!     parameter (chrht1 = .009 , chexp1 = 0.87 )
!     parameter (chrht2 = .013 , chexp2 = 0.90 )
!     parameter (chrht3 = .016 , chexp3 = 0.98 )
!     parameter (chrht4 = .018 , chexp4 = 1.28 )
 
#ifndef NCAR_DUMMY

! Set library entry routine name
 
        if (tgname .eq. 'NONE') tgname = 'TGINIT'
 
 
 
! If the this is the first time the initialization routine was called,
! set up the defaults of the library.
 
        if (doinit) then
 
          doinit = .FALSE.
          dotrans = .FALSE.
 
! set global variables
 
          ichfont = deffnt
          numpnt = 0
          autofeed = .TRUE.
 
! set the default mapping
 
          call map (0.d0 ,1.d0 ,0.d0 ,1.d0 ,0.d0 ,1.d0 ,0.d0 ,1.d0 )
 
! tv80lib defaults to no clipping, turn it off
 
          call dders (1)
 
! set line attributes
 
          call setpch (0,1,ichfont,100)
 
! set the colormap and current color for the device
 
          call colora ('WHITE')
 
! set default font, tv80 font size index 0, rotation 0, clipping off
 
          call setch (0.0d0 ,0.0d0 ,0,0,0,ichfont)
          call tgchclip (0)
 
! set text path to left to right and allignment to bottom left
 
          call gstxp (0)
          call gstxal (1,5)
 
! set transformation to none
 
        call init2d ()
 
        endif
 
! Alter the workstation to conform to tv80 by centering the window
 
        call gqwkc (iwkid,ierr,dummy,itype)
        i_one = 1
        call gqdsp (itype,ierr,i_one,rx,ry,lx,ly)
 
        if (rx .gt. ry) then
          diff = REAL(rx - ry)/2. 
          call gswkvp (iwkid,diff,REAL(ry)+diff,0. ,REAL(ry))
        else if (rx .lt. ry) then
          diff = REAL(ry - rx)/2. 
          call gswkvp (iwkid,0. ,REAL(rx),diff,REAL(rx)+diff)
        endif
 
! set up to 16 colors
 
        call gqeci (iwkid,0,ierr,numclr,dummy)
        if (numclr .eq. 0) numclr = 16
        do 10 i = 1,min(numclr,16)
          call gscr (iwkid,i-1,red(i),green(i),blue(i))
10      continue
 
        if (tgname .eq. 'TGINIT') tgname = 'NONE'
#endif

        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
