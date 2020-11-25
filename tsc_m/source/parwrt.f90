        subroutine  parwrt (nout, np1, nr1, ntest)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nr1,ntest,nout
!============
      character*8 ne1,np1
        ne1 = '        '
        if (nr1.gt.ntest)   ne1 = '<-- yes '
        write (nout,1602)   np1,nr1,ntest,ne1
 1602   format (1x,a8,i5,5x,i5,a8)
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
