!
!c
!     LowOf1(ch); LSCPause; MaxOfAry; MinOfAry ------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE LowOf1(ch)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*1 ch
      INTEGER int
      int = ichar(ch)
      if (int .le. 90 .and. int .ge. 65) int = int + 32
      ch = char(int)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
