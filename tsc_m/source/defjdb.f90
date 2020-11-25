      subroutine defjdb
!
!...subroutine to calculate <j.B>/<B**2> as a function of
!...psi for equilibrium calculations
!
      USE CLINAM
      USE SAPROP
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 arg,gval,gpval,gppval
!============
      if(ifunc .ne. 7) return
!
      do 100 j=2,npsit
      arg=abs((xsv2(j)-psimin)/(psilim-psimin))
      if(arg .gt. 1.0_R8) arg=0.999999_R8
      ajaveq(j)=acoef(991)*(1.0_R8-arg**acoef(992))**acoef(993) +        &  
     &(1.0_R8-acoef(991))*(acoef(994)**2*arg**acoef(995)*                &  
     &(1.0_R8-arg)**acoef(996))/(acoef(994)**2+(acoef(997)-arg)**2)
      call geval(xsv2(j),2,gval,gpval,gppval,imag,jmag)
      ajaveq(j)=ajaveq(j)*gval*xmja2(j)/(vp2(j)*bsqar(j))
 100  continue
      ajaveq(1)=ajaveq(2)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
