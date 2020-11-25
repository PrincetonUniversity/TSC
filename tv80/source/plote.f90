        subroutine plote
!*************************************************************************
!
!  tginit.F  -  initialization and necessary routines in TV80GKS
!
!  contents     plote -  end plotting
!        plotea -  flush buffers
!               endpl -  flush buffer and reset mapping
!               frame -  Flush buffers and frame advance
!               colori -  Set the current color by index
!               colora -  Set the current color by name
!               tginit -  Initialize the interface
!        tgreset -  Reinitialize the interface
!               tgdefs - block data for tv80 to gks shell
!
!    National Energy Research Supercomputer Center
!    Lawrence Livermore National Laboratory
!    University of California
!    Livermore, California 94550, USA
!
!    (C) Copyright 1991 The Regents of the University of California.
!    All Rights Reserved.
!
!*************************************************************************
 
!*************************************************************************
!
!  plote  -  end plotting
!
!  synopsis     call plote
!
!  description  Flush the line buffer, deactivate all workstations that
!               are currently active, close all that are open and then
!               close gks.
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
      INTEGER istate,ierr,numws,iws
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
 
#ifndef NCAR_DUMMY 
 
! Set library entry routine name
 
        if (tgname .eq. 'NONE') tgname = 'PLOTE'
 
 
 
        call tgflush
 
! While GKS is in the state of workstation active (3), get the workstation
! identifier and then deactivate the workstation.  When all workstations
! are deactivated, the state will be workstation open (2).
 
10    call gqops (istate)
      if (istate .eq. 3) then
        call gqacwk (1,ierr,numws,iws)
        call gdawk (iws)
        goto 10
      endif
 
! While GKS is in the state of workstation open (2), get the workstation
! identifier and then close the workstation.  When all workstations
! are closed, the state will be workstation gks open (1).
 
20    call gqops (istate)
      if (istate .eq. 2) then
        call gqopwk (1,ierr,numws,iws)
        call gclwk (iws)
        goto 20
      endif
 
        call gclks
 
        doinit = .TRUE.
 
        if (tgname .eq. 'PLOTE') tgname = 'NONE'
        return
 
! Some sort of disaster happened. Do an emergency close and return.
 
999     call geclks
 
! Reset init variable so if the user wants to use the library without
! restarting his program, the library will initialize again.
 
        doinit = .TRUE.
 
        if (tgname .eq. 'PLOTE') tgname = 'NONE'
#endif

        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
