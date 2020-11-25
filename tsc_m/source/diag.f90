      subroutine diag(array,nameofarray)
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 array
!============
      dimension array(pnx,pnz)
      character*10 nameofarray
!============      
!
!     write(nterm,6670) nameofarray
!6670 format(a20,/)
!     do j=17,15,-1
!
!     write(nterm,6666)j, (array(i,j),i=10,15)
!6666 format(" j=",i3,"i=10,15",1p6e12.4)
!     enddo
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
