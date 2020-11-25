!
!----------------------------------------------------------------------
!
      SUBROUTINE GNfini
      USE gnuI
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      if ( level .ge. 2 ) then
           write (iUnit,'(a)') 'set nomultiplot # by GNfini '
           write (iUnit,'(a)') 'set title       # by GNfini '
           level=1
      else
           write (iUnit,'(a)') '# end of file #'
           close (iUnit)
           level=0
      endif
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
