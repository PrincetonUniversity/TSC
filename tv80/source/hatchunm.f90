      subroutine hatchunm (xmap, ymap)
 
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     H A T C H U N M A P
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      USE hacomm
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     REAL   sinphi,cosphi
!============
      REAL*8 xmap, ymap, xunmap, yunmap
!     common /hacomm/ sinphi,cosphi
 
! Rotate through an angle of phi to original orientation.
 
      xunmap = cosphi*xmap - sinphi*ymap
      yunmap = sinphi*xmap + cosphi*ymap
      xmap = xunmap
      ymap = yunmap
      return
 
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
