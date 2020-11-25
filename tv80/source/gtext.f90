        subroutine gtext (string,num,ioffset)
 
!*************************************************************************
!
!  gtext  -  plot a string
!
!  synopsis     call gtext (string)
!               call gtext (string,num)
!               call gtext (string,num,ioffset)
!
!  description  Plots string of num characters starting at ioffset.
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
      INTEGER ierr,kepclp
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   cliprect
!============
        character*(*) string
        integer num,ioffset
 
        integer inum,ioff
        character*256 strbuf
        dimension cliprect(4)
      REAL   x,y
 
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
 
#ifndef NCAR_DUMMY 
 
! Set library entry routine name
 
        if (tgname .eq. 'NONE') tgname = 'GTEXT'
 
 
 
        if (num .gt. 256) then
          call tgerror (7)
          goto 9999
        endif
 
! Get the parameters into local variables.  The defaults are the length
! of the string and an offset of 0
 
        if (num .eq. -1) then
        inum = len (string)
      else
        inum = num
      endif
 
      if (ioffset .eq. -1) then
        ioff = 1
      else
        ioff = ioffset + 1
      endif
 
! Convert the string to the appropriate case and save in strbuf
 
        if (ichcase .eq. 0) then
          call tgtoup (string(ioff:ioff+inum-1),strbuf,inum)
        elseif (ichcase .eq. 1) then
          call tgtolow (string(ioff:ioff+inum-1),strbuf,inum)
        else
          strbuf = string(ioff:ioff+inum-1)
        endif
 
! Reposition the currrent x,y in normalized coordinates if needed
 
        if (autofeed .and. (chx .gt. 1.0 .or. chx .lt. 0.0 .or.          &  
     &      chy .gt. 1.0 .or. chy .lt. 0.0 )) then
          if (chx .gt. 1.0 ) then
            chx = abs (chaddx)
          else if (chx .lt. 0.0 ) then
            chx = 1.0 - abs (chaddx)
          endif
          if (chy .gt. 1.0 ) then
            chy = abs (chaddy)
          else if (chy .lt. 0.0 ) then
            chy = 1.0 - abs (chaddy)
          endif
          call frame (-1)
        endif
 
! Turn off clipping so the labels will appear if character clipping is off
 
        if (ichclip .eq. 0) then
          call gqclip (ierr,kepclp,cliprect)
          call gsclip (0)
        endif
 
! Write string and update the current x,y in normalized coordinates
 
        if (dotrans) then
          call cntxy (chx,chy,x,y)
          call gtx (x,y,strbuf(1:inum))
        else
          call gtx (chx,chy,strbuf(1:inum))
        endif
 
! Return clipping to its original state if character clipping is off
 
        if (ichclip .eq. 0) call gsclip (kepclp)
 
        chx = chx + chaddx
        chy = chy + chaddy
 
9999    if (tgname .eq. 'GTEXT') tgname = 'NONE'
#endif

        return
 
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
