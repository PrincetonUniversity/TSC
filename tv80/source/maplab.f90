!#include "f77_dcomplx.h"
      subroutine maplab (fmin, fmax, islog, isx, itic)
!=======================================================================
! 12/19/2001 ler pppl - fix rounding when shifting digits by dividing by
!                        10.0d0 instead of 10.0
!=======================================================================
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER islog,isx,itic,jlimch,jmaxch,ierr,kepclp,it,jdxe
      INTEGER i1,inta,icomp,jzeros,jstare,jzeroe,jende,jnice
      INTEGER jmaxe,jprint,kexp,n,i,jdatae,jdshl,jprinn,jnp,jnp1,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   fmax,fmin,cliprect,ccpl1,vl,vr,vb,vt,wl,wr,wb,wt
      REAL   fscale,wl8,wr8,wb8,wt8,gdx,gstart,gend,diff,g1,gexp1
      REAL   gsavsp,expsta,data1
      REAL   REAL, log10a
!============
      REAL*8 x1, y1, x2, y2
!
!c purpose: set number of divisions and get ascii value
!           that should be displayed at tick marks
!
!c revised: oct20.1980 by steven williams
!
!c arguments:
!
! fmin = minimum x or y axis value
!
! fmax = maximum x or y axis value
!
! islog = axis log indicator
!      = 1 (do log scaling)
!      = 0 (do linear scaling)
!
! isx = x or y axis
!      = 2 (polar x axis)
!      = 1 (x axis)
!      = 0 (y axis)
!
! itic = print tic marks
!      = 1 (print tics)
!      = 0 (omit tics)
!
!c enhancements:
!
! gdx = distance between successive tick marks
! gend = axis end point
! gsavsp = saved starting point coordinate
! gstart = axis start point
! jblaw = word filled with ascii blank characters
! jdshl = number of character positions to shift decimal point left
! jdxe = integer exponent of gdx
! jende = integer exponent of gend
! jlimch = "limit characters"
!        = upper limit on desirable number of characters to print
!          in tick mark labels
! jmaxch = "maximum characters"
!        = maximum number of characters that can be printed
!          in a tick mark label
! jmaxe = maximum of jstare and jende
! jnice = 0 (=> mean data - need an exponent factor at bottom of axis)
!       = 1 (=> nice data - dont need an exponent factor at bottom of axis)
! jprinn = number of trailing characters not to print in tick labels
! jprint = base for number of characters to print in tick labels
! jsavst = saved string
! jstare = integer exponent of gstart
! jzeroe = "zero end"
!        = -1 (=> gend .lt. 0.0)
!        = 0  (=> gend .eq. 0.0)
!        = 1  (=> gend .gt. 0.0)
! jzeros = "zero start"
!        = -1 (=> gstart .lt. 0.0)
!        = 0  (=> gstart .eq. 0.0)
!        = 1  (=> gstart .gt. 0.0)
! ndotz = "dot zeros"
!       = a 2-word character string starting with a dot
!         and followed by zeros
!
!c variable declarations:
!
!        integer kascii(3),kprint,kskip,knp
        integer kprint,kskip,knp
      character*20 kstring
        REAL   cdata,cspoin
!
!      dimension ndotz(2)
      character*20 ndotz
!      dimension ns(10)
      character*20 ns
      dimension cliprect(4)
!
      save
!      data jblaw /"        "/
      data jlimch /6/
! can only have 12 digits to right of decimal point with baselib
! subroutine zcetoa - need 7 more for sign, one digit to left of
! decimal point, decimal point and e+## trailer phrase => 19 characters
      data jmaxch /18/
      data ndotz /'.0000000000000000000'/
      data ccpl1 /42.5 /

#ifndef NCAR_DUMMY

!
! Turn off clipping so the labels will appear
!
        call gqclip (ierr,kepclp,cliprect)
        call gsclip (0)
!
! Do the initial setup and draw a border around the window
!
        call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
        call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,1)
 
        if (isx .eq. 2) then
          knp = (vr-vl) * 11. 
          cspoin = ccpl1*(vt+vb)-5.0 
          fscale = .5 * (vt-vb)/(wr+wr)
          call line4 (wl,(wt+wb)/2. ,wr,(wt+wb)/2. )
          call line4 ((wr+wl)/2. ,wb,(wr+wl)/2. ,wt)
        else
          if (isx .eq. 1) then
            knp = (vr-vl) * 13. 
            cspoin = 2*ccpl1*vb-4.0 
            fscale = (vr-vl)/(wr-wl)
          else
            knp = (vt-vb) * 13. 
            cspoin = 2*ccpl1*vl-4.0 
            fscale = (vb-vt)/(wb-wt)
          endif
          wl8 = wl
          wr8 = wr
          wb8 = wb
          wt8 = wt
          call setcrt4 (wl,wb)
          call vector4 (wr,wb)
          call vector4 (wr,wt)
          call vector4 (wl,wt)
          call vector4 (wl,wb)
        endif
 
        if (cspoin .lt. 1) cspoin = 1.0 
 
! jump to 5000 if log axis scaling is to be performed
 
      if (islog .eq. 1) go to 5000
!
      if (knp .le. 0) knp = 1
      gdx = (fmax-fmin)/knp
!
! get the floor of the step size
      jdxe = log10a(gdx)
!
! we will only use step sizes of 0.1, 0.2 and 0.5
      i1 = gdx*(10.0 **(-jdxe))+0.000001 
      if (i1-5) 1100, 1400, 1200
 1100 if (i1-2) 1400, 1400, 1300
 1200 i1 = 5
      go to 1400
!
 1300 i1 = 2
 1400 gdx = i1*(10.0 **jdxe)
!
! get the starting position of the text
      gstart = inta(fmin/gdx)*gdx
      if (icomp(gstart,fmin) .eq. -1) gstart = gstart+gdx
!
! get the new number of divisions
      knp = inta((fmax-gstart)/gdx)
      go to 1550
 1500 knp = knp+1
 1550 gend = knp*gdx+gstart
      if (knp .gt. 34) go to 1600
      if (icomp(gend,fmax)) 1500, 1700, 1600
 1600 knp = knp-1
!
! get exponent of start point
 1700 jzeros = icomp(gstart,0.0 )
      if (jzeros .eq. 0) jstare = 0
      if (jzeros .ne. 0) jstare = log10a(gstart)
!
! get exponent of end point
      jzeroe = icomp(gend,0.0 )
      if (jzeroe .eq. 0) jende = 0
      if (jzeroe .ne. 0) jende = log10a(gend)
!
! specify mean or nice numbers
      jnice = 1
      if ((jdxe .le. -(jlimch))                                          &  
     & .or. (jdxe .eq. -(jlimch-1) .and. gstart .lt. 0.0 )               &  
     & .or. (jdxe .eq. -(jlimch-1) .and. gend .lt. 0.0 )                 &  
     & .or. (jstare .eq. (jlimch) .and. gstart .lt. 0.0 )                &  
     & .or. (jende .eq. (jlimch) .and. gend .lt. 0.0 )                   &  
     & .or. ((jlimch) .lt. jstare)                                       &  
     & .or. ((jlimch) .lt. jende))                                       &  
     & jnice = 0
!
! get maximum number of characters to print
      jmaxe = max(jstare,jende)
      if (jnice .eq. 0 .and. jmaxe .eq. 0 .and. jzeros .eq. 0)           &  
     & jmaxe = jende
      if (jnice .eq. 0 .and. jmaxe .eq. 0 .and. jzeroe .eq. 0)           &  
     & jmaxe = jstare
      if (jnice .eq. 0) then
         jprint = jmaxe-jdxe+1
      elseif (jnice .gt. 0) then
         if (jmaxe .lt. 0) then
            jprint = 1-jdxe
         elseif (jdxe .gt. 0) then
            jprint = jmaxe+1
         else
            jprint = jmaxe-jdxe+2
         endif
      endif
!     add 1 character for leading minus sign
      if (gstart .lt. 0. .or. gend .lt. 0. ) jprint = jprint+1
!
! output warning message
      if (jprint .gt. jmaxch) then
         jprint = jmaxch
         print *,'warning: can only print ',jmaxch,                      &  
     &           ' characters in tick mark labels'
      endif
!
      if (jnice .ge. 1) go to 2000
!
! get step exponent value
!
! divide the original calculation '10.0**float(jdxe)' by 10 in order
! for the internal write to write out the same exponent as zcetoa.
!
      diff = 1.0 /(ccpl1*fscale)
      cdata = gstart-diff
      g1 = (10.0 **REAL(jdxe)) / 10.0 
!     g1 = 10.0d0**(jdxe-1)
      gexp1 = 1.00001 *g1
!
! convert floater [gexp1] to ascii e<jmaxch>.<jmaxch-7> format string
!      kascii(1) = jblaw
!      kascii(2) = jblaw
!      kascii(3) = jblaw
!      call zcetoa (kascii(1), 0, gexp1, jmaxch, jmaxch-7)
      kstring = ' '
      write (kstring,'(E18.11)') gexp1
!
! print step exponent value with no tic or grid marks
      kprint = 4
      kskip = jmaxch-4
      kexp = 1
      cspoin = cspoin-2.0 
      call mapdrw (kstring,kprint,kskip,cdata,cspoin,isx,2)
      cspoin = cspoin+2.0 
!
! now go through the number of divisions loop
 2000 n = knp+1
      do 4099 i = 1, n
!
! get the ascii value to be plotted
      g1 = -(i-1)*gdx
      cdata =(gstart-g1)*1.000001 
!
      if (icomp(gstart,g1) .eq. 0) then
!
!        get (dummy) exponent for zero
         jdatae = min(jdxe,0)
!
!        get character code for zero
!
! I think that the original code actually generated one more character
! than necessary.
!
!         kascii(1) = 8h 0.00000
!         kascii(2) = 8h0000000e
!         kascii(3) = 8h+00
         kstring = ' 0.00000000000e+00'
      else
!
!        get exponent for data
         jdatae = log10a(cdata)
!
!        get character code for data
!        convert floater [cdata] to ascii e<jmaxch>.<jmaxch-7> format string
!         kascii(1) = jblaw
!         kascii(2) = jblaw
!         kascii(3) = jblaw
!         call zcetoa (kascii(1), 0, cdata, jmaxch, jmaxch-7)
!
! This code isn't clear in some parts what it is doing.  Convert the
! number to the same format that zcetoa uses so that code changes
! aren't necessary.
!
         kstring = ' '
         write (kstring,'(E18.11)') cdata / 10.0 
         kstring(2:jmaxch-5) = kstring(3:jmaxch-4)
         kstring(2:2) = kstring(3:3)
         kstring(3:3) = '.'
      endif
!
! right justify character code for data
      if (jnice .ge. 1) go to 2600
!
! have mean data
!      call zmovechr (kascii, 2, kascii, 3, jmaxch-7)
!
! This appears to move the fractional portion over the decimal.
! so 1.23456789012e+12 becomes 1234567890122e+12
!
      kstring(3:jmaxch-5) = kstring(4:jmaxch-4)
      if (icomp(cdata,0. ) .eq. 0) then
         jprint = 1
      else
         jprint = jdatae-jdxe+1
      endif
      go to 3000
!
! have nice data
 2600 if (jdatae .le. -1) go to 2700
!
! nice data has non-negative exponent
!
! This appears to do the same as the above zmovechr except with 'nice'
! data.  So  1.23456789012e+12 may become 1234567890122e+12
!
!      call zmovechr (kascii, 2, kascii, 3, jdatae)
      kstring(3:2+jdatae) = kstring(4:3+jdatae)
!      call zmovechr (kascii, jdatae+2, ndotz, 0, 1)
      kstring(jdatae+3:jdatae+3) = ndotz(1:1)
      go to 2800
!
! nice data has negative exponent
 2700 jdshl = -jdatae-1
!      call zmovechr (ns, 0, kascii, 0, jmaxch)
      ns(1:jmaxch) = kstring(1:jmaxch)
!      call zmovechr (kascii, jdshl+3, ns, 3, jmaxch-7-jdshl)
      kstring(jdshl+4:jmaxch-7) = ns(4:jmaxch-jdshl-4)
!      call zmovechr (kascii, jdshl+2, ns, 1, 1)
      kstring(jdshl+3:jdshl+3) = ns(2:2)
!      call zmovechr (kascii, 1, ndotz, 0, jdshl+1)
      kstring(2:jdshl+2) = ndotz(1:jdshl+1)
!
! shift data right
 2800 jprinn = jprint-(max(jdatae,-1)+1)+min(jdxe,0)
      if (jdxe .le. -1) jprinn = jprinn-1
!      call zmovechr (ns, 0, kascii, 0, jmaxch)
      ns(1:jmaxch) = kstring(1:jmaxch)
!      call zmovechr (kascii, jprinn, ns, 0, jprint)
      kstring(jprinn+1:jprinn+jprint) = ns(1:jprint)
!      call zmovechr (kascii, 0, jblaw, 0, jprinn)
      kstring(1:jprinn) = ' '
!
! get number of significant digits to print for data
 3000 kprint = jprint
      if (cdata .lt. 0.0 ) kprint = kprint+1
!
! get number of characters to skip for data
      kskip = 1
      if (cdata .lt. 0.0 ) kskip = 0
!
! indicate no exponent
 3100 kexp = 0
!
! print data
      gsavsp = cspoin
      cspoin = cspoin-REAL(kprint)+4. 
      if (cspoin .lt. 1.0 ) cspoin = 1.0 
      call mapdrw (kstring,kprint,kskip,cdata,cspoin,isx,itic)
!
! restore cspoin value
      cspoin = gsavsp
!
 4099 continue
!
      goto 9999
!
! this is log scaling
 5000 jnp = (2.0 *knp)/(fmax-fmin)
      if (jnp .gt. 8) jnp = 8
      kexp = 0
!
! now go through the number of divisions loop
      i1 = fmax-fmin+1
      do 6099 i = 1, i1
!
! write out exponent with line
      expsta = fmin+i-1.0 
      cdata = expsta
      data1 = 10.0 **cdata
!
!     cft-dependent encode statement has been converted to baselib calls
!     encode (12, 111, kascii) data1
!
!     convert floater [data1] to ascii e12.5 format string
!      call zmovechr (kascii(1), 0, "                ", 0, 16)
!      call zcetoa (kascii(1), 0, data1, 12, 5)
!
! Only the exponent seems to be of interest here.  Divide by 10
! to get the same exponent returned from zcetoa.
!
      kstring = ' '
      write (kstring,'(E12.5)') data1 / 10.0 
!
      kprint = 4
      kskip = 8
      call mapdrw (kstring,kprint,kskip,cdata,cspoin,isx,itic)
!
      if (jnp .le. 0) go to 6099
!
! specify to print out rightmost character in kstring array with line
      kprint = 1
      kskip = 7
      cspoin = cspoin+3.0 
!
      if (i .eq. i1) goto 9999
!
! now go through parts loop
!
! The following call to mapdrw seems to only print one character at
! a time.  The previous internal write does not seem to influence
! this loop.
!
      jnp1 = jnp+1
      do 7099 j = 2, jnp1
      cdata = expsta+alog10(REAL(j))
!      kascii(1) = 48+j
      kstring(1:8) = '        '
      kstring(8:8) = char (48+j)
      call mapdrw  (kstring,kprint,kskip,cdata,cspoin,isx,itic)
 7099 continue
!
      cspoin = cspoin-3.0 
!
 6099 continue
!
! Reset the map to what it was before entering maplab
!
 9999   call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
        call gsclip (kepclp) 
#endif

        return
        END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
