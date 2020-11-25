      subroutine delpdef
!
      USE CLINAM
      USE SAPROP
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pplas
!============
      if(lrswtch.gt.0 .and. psimin.ge.psilim) psilim = psimin+1._R8
      pplas = psilim - psimin
      if(pplas.ne.0) delpsi = 1._R8/pplas
!
      if(kcycle.le.0 .or. pplash.eq.0) pplash=pplas
      phalo = psilim + whalos*pplash
      if(whalos .ge. 100._R8) phalo = 1.E6_R8
!
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
