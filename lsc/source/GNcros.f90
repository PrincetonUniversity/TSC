!
!----------------------------------------------------------------------
!
      SUBROUTINE GNcros ( xs, ys, n )
      USE gnuI
      USE gnuR
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, n
      REAL*8    xs(n), ys(n), x, y
      INTEGER ixplot, iyplot, size
      DATA    size / 1 /
 
      call GNzero(x,y,n)
 
      write (iUnit,'(a)') "plot '-' with points"
      do i = 1, n
         write( iUnit, '(1x, 1pe9.2, 1x, 1pe9.2)' ) xs(i), ys(i)
      enddo
      write (iUnit,'(a)') 'e'
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
