!
!----------------------------------------------------------------------
!
      SUBROUTINE GNcurv(x,y,n)
      USE gnuI
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, n
      REAL*8    x(n), y(n)
 
      call GNzero(x,y,n)
 
      write (iUnit,'(a)') "plot '-' with lines"
      i = 1
         write( iUnit, '(1x, 1pe9.2, 1x, 1pe9.2)' ) x(i), y(i)
 
      do i = 2, n-1
         if ( y(i)   .ne. 0.00_R8.or.                                    &  
     &        y(i-1) .ne. 0.00_R8.or.                                    &  
     &        y(i+1) .ne. 0.00_R8) then
         write( iUnit, '(1x, 1pe9.2, 1x, 1pe9.2)' ) x(i), y(i)
         endif
      enddo
 
      i = n
         write( iUnit, '(1x, 1pe9.2, 1x, 1pe9.2)' ) x(i), y(i)
 
      write (iUnit,'(a)') 'e'
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
