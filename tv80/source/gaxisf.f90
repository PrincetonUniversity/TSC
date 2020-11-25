      subroutine gaxisf (fxb, fyb, fxe, fye, iori, isize, iside          &  
     &  , form, idiv, adivin)
 
 
!**************************************************************************
!
!  gaxisf.F - axis labeling routine
!
!  contents    gaxisf - axis parameter calculations
!        gaxdrw - axis label drawing from calculations
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
      INTEGER ioldhz,ioldvt,j,iabs,j20,j10
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   cliprect,bdiv,bpp,bsk,bsl,side,ss,divs,gq,gdx,gux
      REAL   g1,start,gxs,gys,gds
      REAL   end
!============
      REAL*8 fxb, fyb, fxe, fye, dcxb, dcyb, dcxe, dcye
 
      character*(*) form
      REAL*8 adivin(*)
 
!      real adivin(*)
 
      integer lenfor
      character*80 outstr
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
      dimension bdiv(30)
      dimension bpp(3,35)
      dimension bsk(4)
      dimension bsl(4)
 
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
 
!
      save bsk
      data bsk(1) /17. /
      data bsk(2) /12. /
      data bsk(3) /9. /
      data bsk(4) /7. /
!
      save bsl
      data bsl(1) /128. /
      data bsl(2) /85. /
      data bsl(3) /64. /
      data bsl(4) /42. /
 
! Set library entry routine name
 
      if (tgname .eq. 'NONE') tgname = 'GAXISF'
 
 
 
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
! This version of gaxis requires floating point input for the labels
!
! Check the input format and get the length of the format.  If 'F10' is
! passed, then the length is 10.
!
      if (form(1:1) .ne. 'f' .and. form(1:1) .ne. 'F' .and.              &  
     &    form(1:1) .ne. 'e' .and. form(1:1) .ne. 'E') then
        print *,'GAXISF: invalid format'
      goto 999
      endif
 
      kn1 = 0
      do 3 j = 2,lenfor
        if (form(j:j) .eq. '.') goto 4
        if (form(j:j) .lt. '0' .or. form(j:j) .gt. '9') then
        print *,'GAXISF: invalid format'
        goto 999
      endif
      kn1 = kn1 * 10 + ichar (form(j:j)) - 48
 3    continue
 4    frmstr = '(' // form(1:lenfor) // ')'
!
! jump if the format was bad
!
      if (kn1 .eq. 0) go to 60
      if (kn1 .gt. 80) go to 60
!
! jump according to the type of labeling to be done
      if (kdiv) 10, 20, 30
!
! use the array bdiv and its information to plot the label
!
   10 kdiv = iabs(kdiv)
!
      do 11 j = 1, kdiv
        bdiv(j) = adivin(j)
   11 continue
!
! fill the value array (bpp), position x array (cpx)
! and the position y array (cpy)
      divs = bdiv(kdiv)-bdiv(1)
      if (divs .le. 0) go to 63
      do 15 j = 1, kdiv
         bpp(1,j) = bdiv(j)
         gq = (bdiv(j)-bdiv(1))/divs
         cpx(j) = (cxe-cxb)*gq+cxb
         cpy(j) = (cye-cyb)*gq+cyb
   15 continue
      go to 40
!
! generate nice numbered labels between bdiv(1) and bdiv(2)
   20 kdiv = max(cxma-cxmi,cyma-cymi)*bsk(ksize+1)
      bdiv(1) = adivin(1)
      bdiv(2) = adivin(2)
   21 gdx = (bdiv(2)-bdiv(1))/kdiv
! jump if bdiv(2) comes before bdiv(1)
      if (gdx .le. 0) go to 63
      gux = alog10(gdx)
      if (gux .lt. 0. ) gux = gux-1. 
! get the floor of the step size
      j20 = gux
! only use steps of .1, .2, or .5
      g1 = gdx*(10.**(-j20))
      j10 = g1+0.5 
      if (j10-5) 22, 25, 23
   22    if (j10-2) 25, 25, 24
   23 j10 = 5
      go to 25
   24 j10 = 2
   25 gdx = j10*(10.**j20)
! change gdx if it would create more than kdiv
! labels.
      if (gdx*kdiv .le. bdiv(2)-bdiv(1)) gdx = 2*gdx
! find the starting point. must be the first nice number
! either on or after the specified starting point
      start = int(bdiv(1)/gdx)*gdx
!     write(*,*) start
      if (start .lt. bdiv(1)) start = start + gdx
! see if we need more than kdiv labels
      kdiv = kdiv/2
   26 end = kdiv*gdx+start
      kdiv = kdiv+1
! the end point has to be the first nice number on or before
! the specified end point
      if (end - bdiv(2)) 26, 27, 265
  265 end = end - gdx
      kdiv = kdiv - 1
! fill the value array (bpp), position x array (cpx)
! and the position y array (cpy)
   27 do 28 j = 1, kdiv
         bpp(1,j) = start+(j-1)*gdx
         gq = (bpp(1,j)-bdiv(1))/(bdiv(2)-bdiv(1))
         cpx(j) = (cxe-cxb)*gq+cxb
         cpy(j) = (cye-cyb)*gq+cyb
   28 continue
      go to 40
!
! generate kdiv labels between bdiv(1) and bdiv(2)
! jump if bad number of divisions
   30 if (kdiv .gt. 30) go to 61
      if (kdiv .lt. 2) go to 62
! calcualte the x and y step sizes
      gxs = (cxe-cxb)/(kdiv-1)
      gys = (cye-cyb)/(kdiv-1)
      bdiv(1) = adivin(1)
      bdiv(2) = adivin(2)
   31 gds = (bdiv(2)-bdiv(1))/(kdiv-1)
! jump if bdiv(2) comes before bdiv(1)
      if (gds .le. 0) go to 63
      bpp(1,1) = bdiv(1)
      cpx(1) = cxb
      cpy(1) = cyb
! fill the value array (bpp), position x array (cpx)
! and the position y array (cpy)
      do 32 j = 2, kdiv
         bpp(1,j) = bpp(1,j-1)+gds
         cpx(j) = cpx(j-1)+gxs
         cpy(j) = cpy(j-1)+gys
   32 continue
!
! Print out the labels
!
   40 continue
 
! this is the floating point format
      do 45 j = 1, kdiv
       write (outstr,frmstr) bpp(1,j)
         call gaxdrw (cpx(j), cpy(j), outstr)
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
      if (tgname .eq. 'GAXISF') tgname = 'NONE'
#endif

      return
 
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
