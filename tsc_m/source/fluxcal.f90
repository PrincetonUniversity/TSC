      subroutine fluxcal(xcpass,zcpass,ipass,psipass,rmajpass,aminpas)
!-----------------------------------------------------------------------
!......7.10 ff
!
!
!.....special subroutine used in calculating dimensionless
!     gains for rfp option, ie=21,22,23,24
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ipass,n,ii,jj
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zcpass,psipass,rmajpass,aminpas,xcpass,sum1,theta
      REAL*8 delth,xpt,zpt,xx,zz,ppt
!============
      sum1 = 0._R8
      theta = 0._R8
      nn = 100
      delth = 2._R8*pi/nn
!
      do 100 n=1,nn
      theta = theta + delth
      xpt = rmajpass + aminpas*cos(theta)
      zpt =        aminpas*sin(theta)
!
      xx = xpt
      zz = zpt
      if(zz.lt.0 .and. isym.eq.1) zz = -zpt
      ii = (xx-ccon)/deex + 2
      jj = (zz-zzero)/deez + 2
!
      call gf(ineg,nmult,xx,zz,xcpass,zcpass,ppt)
!
      if(ipass.le.4) sum1 = sum1 + cos(theta*ipass)*ppt
      if(ipass.gt.4) sum1 = sum1 + sin(theta*(ipass-4))*ppt
  100 continue
!
      psipass = sum1/nn
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
