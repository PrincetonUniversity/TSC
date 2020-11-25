      subroutine gaxisa (fxb, fyb, fxe, fye, iori, isize, iside          &  
     &  , form, idiv, cdivin)
 
 
!**************************************************************************
!
!  gaxisa.F - axis labeling routine
!
!  contents    gaxisa - axis parameter calculations
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
      USE tvaxis
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER isize,iside,idiv,iori
      INTEGER jside,ierr
      INTEGER ioldhz,ioldvt,j,iabs
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL   cxb,cyb,cxe,cye,ctickx,cticky,coffx,coffy,cac,cbc,cpx
!     REAL   cpy,
      REAL   cliprect,bsl,side,ss,gxs,gys
!============
      REAL*8 fxb, fyb, fxe, fye, dcxb, dcyb, dcxe, dcye
 
      character*(*) form
      character*(*) cdivin(*)
 
      integer lenfor
      character*20 frmstr
!
!c author: mike archuleta dec05.1975
!
!c modifier: steven williams sep26.1980
!
! this subroutine enables the sophisticated user to label
! plots or graphs in just about any method desired. this
! routine is intimately tied to tv80cray and baselib,
! since it assumes the existence of a common block from tv80cray
! (tvgxx1), subroutines setch, crtbcd, and line from tv80cray,
! and subroutines zmovechr, zcetoa, zcftoa, zcitoa and zcotoa
! from baselib.
! assumed inline functions are alog10, sign and sqrt
!
! fxb    the beginning x coordinate of the label
! fyb    the beginning y coordinate of the label
! fxe    the end x coordinate of the label
! fye    the end y coordinate of the label
! iori    the orientation of the text string to be plotted
!         a 0 implies horizontal and a 1 is vertical
! isize   the size of the characters (0-3)
! iside   which side of the axis the label are to be drawn.
!         a 0 implies the left and a 1 implies the right.
! iform   the format of the plotted string. typical uses
!         would be 3ha10, 4hf5.2, 5he20.8, 3hi10, 2ho4
! idiv    the number of labeled marks
! idivin  an array containing information to be plotted
!
! idiv is the key which specifies the type of labeling to do.
! if idiv is negative, then the array adivin contains iabs(idiv)
! values. if the iform is 'a', then the adivin array
! contains the alphanumeric text to be plotted. if the
! iform is not 'a', then the adivin array values are the
! numbers to be used as label at their location.
!
! if idiv is zero, then nice numbered labels are generated
! between adivin(1) and adivin(2).
!
! if idiv is positive, then idiv labels will be plotted
! between adivin(1) and adivin(2).
!
!c variable declarations:
!
!     logical klr
!     common /tvaxis/ cxb, cyb, cxe, cye, ctickx, cticky
!    * , coffx, coffy, cac, cbc, kori, ksize, kn1, klr
!    * , klxsav, klysav, cpx(30), cpy(30), kdiv, kepclp
!
      REAL   cx1,cx2,cy1,cy2,cxmi,cxma,cymi,cyma
      integer itype
      dimension cliprect(4)
      dimension bsl(4)
!
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
 
 
 
!
      save bsl
      data bsl(1) /128. /
      data bsl(2) /85. /
      data bsl(3) /64. /
      data bsl(4) /42. /
!
! Set library entry routine name
 
      if (tgname .eq. 'NONE') tgname = 'GAXISA'
 
 
 
!
!c program statements:
!
! store the arguments into local variables
      cxb = fxb
      cyb = fyb
      cxe = fxe
      cye = fye
      kori = iori
      ksize = isize
      kdiv = idiv
      side = -1. 
      jside = iside
      if (jside .eq. 1) side = 1. 
 
      lenfor = len (form)
!
! Get world and NDC coordinates as well as map type, then set the
! map type to be linear.  Set the allignment to left, centered.
! turn the clipping off so text outside the viewport will show up.
!
      call tggetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,1)
      call gqtxal (ierr,ioldhz,ioldvt)
      call gstxal (1,3)
      call gqclip (ierr,kepclp,cliprect)
      call gsclip (0)
!
! If coordinates in x/y were in log, set klxsav/klysav to 1 else 0
!
      klxsav = 0
      klysav = 0
      if (itype .eq. 3 .or. itype .eq. 4) klxsav = 1
      if (itype .eq. 2 .or. itype .eq. 4) klysav = 1
!
! change the argument which specifies the side of the
! axis the label is to be drawn on if the writing is going
! top to bottom or right to left
!
! The following code uses operands that don't exist in fortran.
! They are rewritten to something that should be equivalent.
!
!      if (kori .eq. 0 .and. cye .lt. cyb) jside = .not. jside
!      if (kori .eq. 1 .and. cxe .lt. cxb) jside = .not. jside
!      klr = .not. (kori .xor. jside)
!
      if (kori .eq. 0 .and. cye .lt. cyb) then
        if (jside .eq. 1) then
        jside = 0
      else
        jside = 1
      endif
      endif
      if (kori .eq. 1 .and. cxe .lt. cxb) then
        if (jside .eq. 1) then
        jside = 0
      else
        jside = 1
      endif
      endif
      klr = (kori .eq. jside)
! calculate the character offset
      coffx = (cx2-cx1)/((cxma-cxmi)*bsl(ksize+1))
      coffy = (cy2-cy1)/((cyma-cymi)*bsl(ksize+1))
! figure out how long the tick marks should be
      ctickx = (cx2 - cx1) * .01 
      cticky = (cy2 - cy1) * .01 
! calculate the normal to the axis
      cac = cyb-cye
      cbc = cxe-cxb
      ss = sqrt(cac*cac+cbc*cbc)
      if (ss .eq. 0) go to 64
! normalize the normal and point it in the right direction
      cac = side*cac/ss
      cbc = side*cbc/ss
!
! This version of gaxis requires integer input for the labels
!
! Check the input format and get the length of the format.  If 'I10' is
! passed, then the length is 10.
!
      if (form(1:1) .ne. 'a' .and. form(1:1) .ne. 'A') then
        print *,'GAXISA: invalid format'
      goto 999
      endif
 
      kn1 = 0
      do 3 j = 2,lenfor
        if (form(j:j) .lt. '0' .or. form(j:j) .gt. '9') then
        print *,'GAXISA: invalid format'
        goto 999
      endif
      kn1 = kn1 * 10 + ichar (form(j:j)) - 48
 3    continue
      frmstr = '(' // form(1:lenfor) // ')'
!
! jump if the format was bad
!
      if (kn1 .eq. 0) go to 60
      if (kn1 .gt. 80) go to 60
!
! use the array bdiv and its information to plot the label
!
   10 kdiv = iabs(kdiv)
      if (kdiv .gt. 30) go to 61
      if (kdiv .lt. 2) go to 62
! space out the ascii information in equal increments
   16 gxs = (cxe-cxb)/(kdiv-1)
      gys = (cye-cyb)/(kdiv-1)
      cpx(1) = cxb
      cpy(1) = cyb
! fill the value array (bpp), position x array (cpx)
! and the position y array (cpy)
      do 17 j = 2, kdiv
         cpx(j) = cpx(j-1)+gxs
         cpy(j) = cpy(j-1)+gys
   17 continue
!
! Print out the labels
!
   40 continue
 
      do 45 j = 1, kdiv
       call gaxdrw (cpx(j), cpy(j), cdivin(j))
   45 continue
!
! draw the axis line
   50 call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call line (cxb*1.0d0 , cyb*1.0d0 , cxe*1.0d0 , cye*1.0d0 )
      goto 9999
!
! these are the error messages
!
   60 call tgerror (2)
      go to 999
   61 call tgerror (3)
      go to 999
   62 call tgerror (4)
      go to 999
   63 call tgerror (5)
      go to 999
   64 call tgerror (6)
!
 999  call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
!
9999  call gsclip (kepclp)
      call gstxal (ioldhz,ioldvt)
      if (tgname .eq. 'GAXISA') tgname = 'NONE'
      return
 
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
