      subroutine deneval(ps,denval)
!
!           6.94 surfplot
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER jj,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 denval,ps,fac1,fac2
!============
      do 100 jj=3,npsit
      j = jj
      if(ps.lt.xsv(j)) go to 101
  100 continue
      go to 102
  101 fac1 = (ps-xsv(j-1))/(xsv(j)-xsv(j-1))
      fac2 = (xsv(j) - ps)/(xsv(j)-xsv(j-1))
      denval = (fac1*densum(j)+fac2*densum(j-1))*udsd
      if(denval .lt. 0) go to 102
      go to 103
  102 denval = 0._R8
  103 return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
