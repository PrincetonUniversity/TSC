!#include "library_names.h"
      subroutine f03aae(a1,n1,n2,a2,a3,n3)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n1,n2,n3
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 a2,d
      REAL*8 a3(*),a1(n1,*)
      INTEGER ipiv(n1)
      INTEGER i
!============

!     using NAG routines
!     call       f03aaf(a1,n1,n2,a2,a3,n3)


!     using simple LU decomposition from Numerical Recipe
!     call ludcmp(a1,n1,n2,ipiv,d)
!     do i = 1, n1
!        d = d * a1(i,i)
!     enddo
!     a2 = d
!     n3 = 0


!     using LAPACK
      call dgetrf(n2,n2,a1,n1,ipiv,n3)
      d = 1.0_R8
      do i = 1, n2
         if(ipiv(i) .ne. i) d = -d
      enddo
      do i = 1, n2
         d = d * a1(i,i)
      enddo
      a2 = d


      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
