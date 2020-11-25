!#include "f77_dcomplx.h"
!     ------------------------------------------------------------------
      SUBROUTINE setpwrl(pwrlevel, n, nflat)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i, nflat, nfl
      INTEGER newn, iALGORIT, iORIGINL, iGEOMETR, iLINEAR
      REAL*8    pwrlevel(n)
      REAL*8                                                             &  
     &        ONE, ONEpt05, ONEpt10, ONEpt15, ONEptEPS
      DATA    ONE, ONEpt05, ONEpt10, ONEpt15, ONEptEPS /                 &  
     &        1.0_R8,    1.05_R8,    1.10_R8,    1.15_R8,     1.15_R8/
! 1.05^100 = 130; 1.10^100 = 10^4; 1.15^100 = 10^6; 1.20^100 = 10^8
 
      DATA    iALGORIT, iORIGINL, iGEOMETR, iLINEAR /                    &  
     &               3,        1,        2,       3 /
      REAL*8 AREAL
 
      if(nflat .lt. 1)then
         nfl = 1
      else
         nfl = nflat
      endif
 
      if(nflat .ge. n-1)then
         nfl = 1
      else
         nfl = nflat
      endif
 
      if (iALGORIT .eq. iORIGINL) then
!
!     Original code.  Raise by 2 * each time
!     BEGIN
!     geometric series
      pwrlevel(1) = 1._R8/ (2._R8** (n - nfl))
      do 10 i = 2, n - nfl + 1
          pwrlevel(i) = 2._R8* pwrlevel(i - 1)
  10  continue
      do 20 i = n - nfl + 1, n
          pwrlevel(i) = pwrlevel(n - nfl + 1)
  20  continue
!     END
!     Original code.  Raise by 2 * each time
 
      else if (iALGORIT .eq. iGEOMETR) then
!     New Geometric Ramp Up
 
      newn = n
      do i = newn, newn-nfl, -1
         pwrlevel(i) = ONE
      enddo
 
      do i = newn-nfl-1, 1, -1
         pwrlevel(i) = pwrlevel(i + 1)/ONEptEPS
      enddo
 
      else if (iALGORIT .ge. iLINEAR) then
!     Linear Ramp Up
      newn = n
      do i = newn, newn-nfl, -1
         pwrlevel(i) = ONE
      enddo
 
      pwrlevel(newn) = ONE
      do i = 1, newn-nfl
         pwrlevel(i) = AREAL(i)/AREAL(newn-nfl)
      enddo
 
      endif
 
      return
      END
!     cycle.F ends                                                     |
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
