!     -----------------------------------------------------------------|
      SUBROUTINE DolAtEnd (SomeString)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*(*) SomeString
      CHARACTER*1 DOLLAR
      INTEGER i,lastch, length
      DATA DOLLAR /'$'/
      length = len(SomeString)
      do 10 i=length,1,-1
	lastch=i
        if (SomeString(i:i) .ne. ' ') go to 11
 10   continue
      SomeString(1:1) = DOLLAR
      return
!
!
 11   continue
      lastch=lastch+1
      if (lastch .gt. length) lastch=length
      SomeString(lastch:lastch) = DOLLAR
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
