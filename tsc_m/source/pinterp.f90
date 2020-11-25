      function pinterp(xx,zz,ii,jj)
!......7.20 pinterp
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!....cubic interpolation version
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ii,jj
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xx,zz,pinterp,xpass,zpass,graddum,dpxdum,dpzdum,gsdum
      REAL*8 psdum,pxzdum,pxxdum,pzzdum
!============
      xpass = xx
      zpass = zz
      if(zz.lt.0 .and. isym.eq.1) zpass = -zz
      call grap(1,zpass,xpass,graddum,dpxdum,dpzdum,gsdum,psdum,         &  
     &          pxzdum,pxxdum,pzzdum,1)
      pinterp = psdum
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
