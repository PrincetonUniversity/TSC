      subroutine tridiag (a,b,c,d,x,e,f,m)
 
!     tri-diagonal matrix solver where m is the number of rows,
!     k = m-1, and matrix equation is as illustrated.
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER m,n,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 b,c,d,x,e,f,a
!============
      dimension a(1),b(1),c(1),d(1),x(1),e(1),f(1)
!============      
 
      e(1) = b(1)
      f(1) = d(1)/e(1)
      do 2 n = 2,m
      e(n) = b(n)-a(n)*c(n-1)/e(n-1)
2     f(n) = (d(n)-a(n)*f(n-1))/e(n)
      x(m) = f(m)
      do 4 i = 2,m
      n = m-i+1
4     x(n) = f(n)-c(n)*x(n+1)/e(n)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
