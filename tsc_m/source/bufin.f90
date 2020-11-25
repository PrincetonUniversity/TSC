      subroutine bufin(lu, xbeg, xfin)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER lu,lnw,loc,i
      INTEGER len
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xbeg,xfin
!============
      dimension xbeg(*)
!============      
!
      lnw = loc(xbeg(2)) - loc(xbeg(1))
      len = (loc(xfin) - loc(xbeg))/lnw + 1
!ccccc      write(*,*) ' len=', len
      read(lu) (xbeg(i),i=1,len)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
