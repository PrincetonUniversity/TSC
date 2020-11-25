!
!----------------------------------------------------------------------
!
      SUBROUTINE GNtitl(Title)
      USE gnuI
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*(*) Title
      CHARACTER*24  MyTitle
      INTEGER         l
      l = INDEX(Title,'$')
      if (l .gt. 1) then
         l = l-1
         write(MyTitle,'(a)') Title(1:l)
      else
         MyTitle = Title
      endif
      write (iUnit,'(''se tit "'', a, ''"'')') MyTitle
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
