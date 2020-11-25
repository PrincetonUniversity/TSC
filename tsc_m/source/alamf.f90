      function alamf(ps)
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ps,alamf,alam,psilima,pmax,psiaa
!============
      alamf = 0._R8
      alam = hypermult*udst
      psilima = hyperfrac*psimin + (1._R8-hyperfrac)*psilim
      if(ps.gt.psilima) return
      pmax = psilima-psimin
      if(pmax.le.0) return
      psiaa=max((ps-psimin)/pmax,1.E-8_R8)
      psiaa = min(psiaa,1.0_R8)**acoef(66)
      alamf = min(1._R8-psiaa,1.0_R8)*alam
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
