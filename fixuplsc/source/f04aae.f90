!#include "library_names.h"
      subroutine f04aae(a1,n1,a2,n2,n3,n4,a3,n5,a4,n6)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n1,n2,n3,n4,n5,n6
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 a1(n1,*),a2(n2,*),a3(n5,*),a4(*)
      INTEGER ipiv(n3)
      INTEGER i, j
      REAL*8 d
!============


!     using NAG routines 
!     call       f04aaf(a1,n1,a2,n2,n3,n4,a3,n5,a4,n6)

 
!!      do i = 1, n3
!!        do j = 1, n4
!!          a3(i, j) = a2 (i, j)
!!        enddo
!!      enddo


!     using ludcmp and lubksb
!     n6 = 0
!     call ludcmp(a1,n1,n3,ipiv,d)
!     do j = 1, n4
!        call lubksb(a1,n1,n3,ipiv,a3(1,j))
!     enddo


!     using LAPACK routines 
      call dgesv(n3,n4,a1,n1,ipiv,a2,n2,n6)
      a3(1:n3,1:n4) = a2(1:n3,1:n4)  


      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
