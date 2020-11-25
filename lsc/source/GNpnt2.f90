!
!----------------------------------------------------------------------
!
      SUBROUTINE GNpnt2(x,y,n)
      USE gnuI
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, n
      REAL*8    x(n), y(n)
 
      call GNzero(x,y,n)
 
      write (iUnit,'(a)') "plot '-' with dots"
      do i = 1, n, 2
         if (y(i) .eq. 0.0_R8) then
         write( iUnit, '(1x, 1pe9.2, 1x, ''0'' )' ) x(i)
         else
         write( iUnit, '(1x, 1pe9.2, 1x, 1pe9.2)' ) x(i), y(i)
         endif
      enddo
 
      write (iUnit,'(a)') 'e'
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
