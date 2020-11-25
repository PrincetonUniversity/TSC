      function qget(ps)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ps,qget
!============
      do 10 j=2,npsit
      if(xsv2(j).gt.ps) go to 20
   10 continue
      qget = qprof2(npsit)
      return
   20 qget = ((ps-xsv2(j-1))*qprof2(j) + (xsv2(j)-ps)*qprof2(j-1))       &  
     &     / (xsv2(j) - xsv2(j-1))
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
