      subroutine points4(xvec4, yvec4, num, incx, incy, delx4, dely4)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER incx,incy,num,i
!============
      REAL*4 xvec4(*), yvec4(*), delx4, dely4
      REAL*8 xvec8(1000), yvec8(1000), delx8, dely8
      do i = 1, num
         xvec8(i) = xvec4(i)
         yvec8(i) = yvec4(i)
      enddo
      delx8 = delx4
      dely8 = dely4
      call points(xvec8, yvec8, num, incx, incy, delx8, dely8)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
