        subroutine tgerror (ierror)
 
!*************************************************************************
!
!  tgerror  -  Print an error message
!
!  synopsis     call tgerror (ierror)
!               integer ierror   Error number
!
!  description  Prints the name stored in tgname, followed by the
!        error message.  Negative error numbers indicate error
!        messages that are somehow dependent on LRLTRAN.  Positive
!        numbers indicate more realistic errors.
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
!     INTEGER maxpnt,ichclip,minfont,maxfont
!============
        integer ierror
 
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
 
 
 
        if (ierror .eq. -1) then
          print *,tgname,': Invalid number of arguments'
 
        else if (ierror .eq. -2) then
          print *,tgname,': Invalid camera id'
 
        else if (ierror .eq. -3) then
          print *,tgname,': Unused workstation identifier not found'
 
        else if (ierror .eq. 1) then
          print *,tgname,': Invalid color'
 
        else if (ierror .eq. 2) then
          print *,tgname,': Bad format'
 
        else if (ierror .eq. 3) then
          print *,tgname,': no more than 30 labels'
 
        else if (ierror .eq. 4) then
          print *,tgname,': no less than 2 labels'
 
        else if (ierror .eq. 5) then
          print *,tgname,': bdiv(2) .le. bdiv(1)'
 
        else if (ierror .eq. 6) then
          print *,tgname,': axis line is a point'
 
        else if (ierror .eq. 7) then
          print *,tgname,': string exceeded maximum of 256'
 
        else if (ierror .eq. 8) then
          print *,tgname,                                                &  
     &      ': Negative coordinates reset to 1e-5 to 1e+5'
 
        else if (ierror .eq. 9) then
          print *,tgname,': Minimum log map value reset'
 
        else if (ierror .eq. 10) then
          print *,tgname,': Invalid value for argument'
 
        else if (ierror .eq. 11) then
          print *,tgname,': Invalid map number'
 
        else if (ierror .eq. 12) then
          print *,tgname,': Bad mapping. left .eq. right'
 
        else if (ierror .eq. 13) then
          print *,tgname,': Bad mapping. bottom .eq. top'
 
        else if (ierror .eq. 14) then
          print *,tgname,': Bad mapping. left .gt. right'
 
        else if (ierror .eq. 15) then
          print *,tgname,': Bad mapping. bottom .gt. top'
 
        else if (ierror .eq. 16) then
          print *,tgname,                                                &  
     &      ': Bad mapping. Log of negative x is undefined'
 
        else if (ierror .eq. 17) then
          print *,tgname,                                                &  
     &      ': Bad mapping. Log of negative y is undefined'
 
        else if (ierror .eq. 18) then
          print *,tgname,                                                &  
     &      ': Bad mapping. Virtual x should be 0 <= x <= 1'
 
        else if (ierror .eq. 19) then
          print *,tgname,                                                &  
     &      ': Bad mapping. Virtual y should be 0 <= y <= 1'
 
        else if (ierror .eq. 20) then
          print *,tgname,': GKS not opened'
 
        else if (ierror .eq. 21) then
          print *,tgname,': Library not initialized'
 
        else if (ierror .eq. 22) then
          print *,tgname,': Insufficient significant bits in x range'
 
        else if (ierror .eq. 23) then
          print *,tgname,': Insufficient significant bits in y range'
 
        endif
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
