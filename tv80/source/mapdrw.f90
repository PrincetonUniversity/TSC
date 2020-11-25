      subroutine mapdrw (string,numchr,ioffset,axis,chrpos,isx,itic)
!**************************************************************************
!
!  mapdrw.f - drawing routine for gaxis
!
!  contents    mapdrw - draws text and tick marks for gaxis
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
 
 
!**************************************************************************
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER it
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,sdata
      REAL   ticlen
!============
      character*(*) string
      integer numchr,ioffset
      REAL   axis,chrpos
      integer isx,itic
      REAL*8 x1, y1, x2, y2
!
! numchr - number of characters in the string
!
! ioffset - offset in the string to the first character
!
! axis - value on the axis to print on
!
! chrpos - The character position at axis
!
! isx - Axis number
!       = 0 (y axis)
!       = 1 (x axis)
!       = 2 (polar axis)
!
! itic - Print tic or grid
!       = 0 (print grid line)
!       = 1 (print tic mark)
!       = 2 (omit marks)
!
 
! Get the current mapping, the type (it) has to be linear
 
      call tggetmap (xvmin,xvmax,yvmin,yvmax,                            &  
     &    xwmin,xwmax,ywmin,ywmax,it)
 
! Check to see if we are plotting on the x axis
 
      if (isx .eq. 1) then
        sdata = (axis-xwmin) * (xvmax-xvmin)/(xwmax-xwmin) + xvmin
        sdata = (sdata * 42.5 ) + 1.5 
        x1 = sdata
        y1 = chrpos
        call setch (x1,y1,0,1,1,-1)
!ccccc          call setch (sdata,chrpos,0,1,1,-1)
        call gtext (string,numchr,ioffset)
        if (itic .eq. 0) then
          x1 = axis
          x2 = axis
          y1 = ywmin
          y2 = ywmax
          call line (x1, y1, x2, y2)
!ccccc            call line (axis,ywmin,axis,ywmax)
        else if (itic .eq. 1) then
          ticlen = .01 * (ywmax-ywmin)
          x1 = axis
          x2 = axis
          y1 = ywmin
          y2 = ywmin+ticlen
          call line (x1, y1, x2, y2)
!ccccc            call line (axis,ywmin,axis,ywmin+ticlen)
          y1 = ywmax
          y2 = ywmax-ticlen
          call line (x1, y1, x2, y2)
!ccccc            call line (axis,ywmax,axis,ywmax-ticlen)
        endif
 
! Check to see if we are plotting on the y axis
 
      else if (isx .eq. 0) then
        sdata = (axis-ywmin) * (yvmax-yvmin)/(ywmax-ywmin) + yvmin
        sdata = (sdata * 42.5 ) + .5 
        x1 = chrpos
        y1 = sdata
        call setch (x1,y1,0,1,0,-1)
!ccccc          call setch (chrpos,sdata,0,1,0,-1)
        call gtext (string,numchr,ioffset)
        if (itic .eq. 0) then
          x1 = xwmin
          y1 = axis
          x2 = xwmax
          y2 = axis
          call line (x1, y1, x2, y2)
!ccccc            call line (xwmin,axis,xwmax,axis)
        else if (itic .eq. 1) then
          ticlen = .01 * (xwmax-xwmin)
          x1 = xwmin
          y1 = axis
          x2 = xwmin+ticlen
          y2 = axis
          call line (x1, y1, x2, y2)
!ccccc            call line (xwmin,axis,xwmin+ticlen,axis)
          x1 = xwmax
          x2 = xwmax-ticlen
          call line (x1, y1, x2, y2)
!ccccc            call line (xwmax,axis,xwmax-ticlen,axis)
        endif
 
! We are not on the x or y axis, we are plotting polar, put ticks in x direction
 
      else
        sdata = (axis-xwmin) * (xvmax-xvmin)/(xwmax-xwmin) + xvmin
        sdata = (sdata * 42.5 ) + 1.5 
        x1 = sdata
        y1 = chrpos
        call setch (x1,y1,0,1,1,-1)
!ccccc          call setch (sdata,chrpos,0,1,1,-1)
        call gtext (string,numchr,ioffset)
        if (itic .eq. 1) then
          ticlen = .01 * (ywmax-ywmin)
          x1 = axis
          y1 = 0. 
          x2 = axis
          y2 = ticlen
          call line (x1, y1, x2, y2)
!ccccc            call line (axis,0.,axis,ticlen)
        endif
      endif
 
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
