!     -----------------------------------------------------------------|
      SUBROUTINE viwrite(n, vec, name)
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,iout,nby5,i,nout,j
!============
      CHARACTER *(*) name
      INTEGER vec, vmin, vmax
      DIMENSION vec(n)
      call VecIMnMx(vec, n, vmin, vmax)
      write(nLSCcom2, 100)name, n, vmin, vmax
 100  format(/, a, /, ' n = ', i3 ,/, ' min = ', i5 ,                    &  
     &             /, ' max = ', i5)
      iout = 1
      if(n .lt. 5)then
              nby5 = 1
      else
              nby5 = n / 5
      endif
      do 20 i = 1, nby5
              nout = n - iout + 1
              if(nout .gt. 5)then
                      nout = 5
              endif
              write(nLSCcom2, 101)(vec(j), j = iout, iout + nout - 1)
 101          format(i5, 2x, i5, 2x, i5, 2x, i5, 2x, i5)
              iout = iout + nout
 20   continue
      write(nLSCcom2, 102)
 102  format(/)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
