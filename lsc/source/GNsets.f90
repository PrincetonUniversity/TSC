!#include "f77_dcomplx.h"
!
!----------------------------------------------------------------------
!
      SUBROUTINE GNsets(i,j,k,l,a,b,c,d,kind)
      USE gnuI
      USE gnuR
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,j,k,l,kind
      REAL*8    a,b,c,d
      REAL*8    xorigin, yorigin
      REAL*8    xsize  , ysize
      REAL*8 AREAL

      xorigin = AREAL(i)  /1024._R8
      yorigin = AREAL(k)  / 768._R8
      xsize   = AREAL(j-i)/1024._R8
      ysize   = AREAL(l-k)/ 768._R8
      if (level .eq. 1) then
         write (iUnit,'(a)') 'set multiplot # by GNsets'
         write (iUnit,'(a)') 'set nolabel   # by GNsets'
         level = 2
      endif
 
      write(iUnit,'(''set origin '' , f5.2,'','', f5.2)')                &  
     &                                             xorigin, yorigin
 
      if (j-i .eq. l-k) then
      write(iUnit,'(''se si sq ''   , f5.2,'','', f5.2)')                &  
     &                                               xsize,   ysize
      else
      write(iUnit,'(''se si nosq '' , f5.2,'','', f5.2)')                &  
     &                                               xsize,   ysize
      endif
 
      write(iUnit,'(''se xr ['',1pe9.2,'':'',1pe9.2,'']'')')             &  
     &                                                            a,b
 
      write(iUnit,'(''se yr ['',1pe9.2,'':'',1pe9.2,'']'')')             &  
     &                                                            c,d
      level=3
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
