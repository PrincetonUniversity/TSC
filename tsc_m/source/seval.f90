      subroutine seval(psv,sval,ig,jg)
!
      USE CLINAM
      USE SAPROP
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ig,jg,jv,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sval,psv,frac
!============
      jv = 2
      frac = 1.0_R8
      if(psv.gt.xsv(2)) then
      do 15 j=3,npsit+1
      frac = (psv-xsv(j-1))/(xsv(j)-xsv(j-1))
      jv = j
      if(frac.ge.0.0_R8.and. frac.le.1.0_R8) go to 21
   15 continue
      frac = 1.0_R8
                        endif
   21 sval = frac*sradion(jv)+(1._R8-frac)*sradion(jv-1)
      if(psv.gt.xsv(npsit+1)) sval = 0._R8
      if(iexvc(ig,jg+1).gt.0) sval = 0._R8
!
      return
   10 continue
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
