      subroutine hatchmap (xuser, yuser)
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     H A T C H M A P
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
      REAL   xuser,yuser,xmap,ymap
!     common /hacomm/ sinphi,cosphi
!
! Rotate through an angle of -phi.
!
      xmap = cosphi*xuser + sinphi*yuser
      ymap = - sinphi*xuser + cosphi*yuser
      xuser = xmap
      yuser = ymap
      return
 
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
