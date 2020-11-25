      subroutine appvolt
!......6.90 cplot
!
      USE CLINAM
      USE SAPROP
      USE SCR3
      USE FEED

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ii,iw,jw,ig,iabs
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 viabturn,viabturno
!============
!     common/feed/ viabturn(18),viabturno(18)
!
      do 100 ii=1,nwire
      iw = iwire(ii)
      jw = jwire(ii)
      etay(iw,jw) = rswire(ii)*dxdz/(tpi*xary(iw))
      ptsave(iw,jw) = -etay(iw,jw)*xary(iw)*cwire0(ii)/dxdz
      ig = iabs(igroupw(ii))
      if(iseries(ig).ne.0) cwire0(ii) = 0._R8
      if(ig.le.18) ptsave(iw,jw) = -viabturn(ig)
      if(ig.eq.6 .or. ig.eq.7  .or. ig.eq.15 .or. ig.eq.16) then
      ptsave(iw,jw) = ptsave(iw,jw) + reboun
                                                            endif
  100 continue
      do 101 ig=1,18
      viabturn(ig) = 0.5_R8*(viabturn(ig)+viabturno(ig))
      viabturno(ig) = viabturn(ig)
  101 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
