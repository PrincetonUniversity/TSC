!#include "f77_dcomplx.h"
!
!----------------------------------------------------------------------
!
      SUBROUTINE GNlist(i,j,string,lcen,lor)
!                                 ,lcen,lor) ignored
      USE gnuI
      USE gnuR
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

!      xsize   = float(j-i)/1024.
!      ysize   = float(l-k)/ 768.
 
      CHARACTER string*(*)
      CHARACTER*1 dollar
      INTEGER     i,j,lcen,lor, indx, MAX, n
      REAL*8        xplace, yplace
      DATA        dollar,   MAX/                                         &  
     &            '$'   ,   75 /
      REAL*8 AREAL
 
      xplace = AREAL(i)/1024._R8
      yplace = AREAL(j)/ 768._R8
 
      n=0
      do 10 indx=1,MAX
        if(string(indx:indx) .eq. dollar) go to 30
        n=n+1
10    continue
 
30    continue
      if(n.eq.0) return
 
100   continue
 
      write(iUnit,'(''se la "'',                                         &  
     &                a,                                                 &  
     &              ''" at screen '',                                    &  
     &                f6.3,                                              &  
     &              '','',                                               &  
     &                f6.3,                                              &  
     &              '' font "Courier"''                                  &  
     &                )')               string(1:n), xplace, yplace
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
