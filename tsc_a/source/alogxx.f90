      function alogxx(a)
!
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 a,alogxx
!============
      if(a.lt.1.0E-99_R8) a = 1.0E-99_R8
      alogxx = log10(a)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
