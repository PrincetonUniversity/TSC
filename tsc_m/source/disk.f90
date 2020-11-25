      subroutine disk(i1,i2,big1,big2)
!.................................................................
!..rxw/end
!
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i2,i1,iadscr,irec
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 big1,big2
!============
      dimension big1(1),big2(1)
!============      
      iadscr=0
      do 450 irec=1,nrecord
      big1(irec) = pltsav(iadscr+i1)
      big2(irec) = pltsav(iadscr+i2)
      iadscr = iadscr + lenscr
  450 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
