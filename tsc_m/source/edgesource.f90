      subroutine edgesource
!**********************************************************************
!......6.92 wplot
!
!.....defines edge source terms for density equation for isurf=1 and idens=1
!
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 arg
!============
      if(acoef(871).le.0) return
      do 10 j=2,npsit
      arg = (j - (npsit+0.5_R8))/(npsit - 1._R8)
      sraveedg(j) = acoef(871)*exp(acoef(872)*arg)*usdd/usdt
   10 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
