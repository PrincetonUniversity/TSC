      subroutine fluxmod(psi1,psi2,psi3,psi4,psisn1,rmajpass,aminpas)
!
!......for use in rfp feedback with ie=21,22,23,24
!
!   ---> note that if this subroutine is changed, then
!        subroutine fluxcal must be changed in a consistent way
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,ii,jj
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 psi2,psi3,psi4,psisn1,rmajpass,aminpas,psi1,sum1,sum2
      REAL*8 sum3,sum4,sum5,theta,delth,xpt,zpt,xx,zz,ppt,pinterp
!============
      sum1 = 0._R8
      sum2 = 0._R8
      sum3=0._R8
      sum4=0._R8
      sum5 = 0._R8
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
      if(isym.eq.1 .and. zz.lt.0) zz = -zpt
      ii = (xx-ccon)/deex + 2
      jj = (zz-zzero)/deez + 2
      ppt = tpi*pinterp(xx,zz,ii,jj)
!
      sum1 = sum1 + cos(theta)*ppt
      sum2 = sum2 + cos(2._R8*theta)*ppt
      sum3 = sum3 + cos(3._R8*theta)*ppt
      sum4 = sum4 + cos(4._R8*theta)*ppt
      sum5 = sum5 + sin(theta)*ppt
  100 continue
!
      psi1 = sum1/nn
      psi2 = sum2/nn
      psi3 = sum3/nn
      psi4 = sum4/nn
      psisn1 = sum5/nn
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
