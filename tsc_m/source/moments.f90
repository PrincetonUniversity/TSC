      subroutine moments(bth,brho)
!.....3.95 moments
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,iflip
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 brho,bth,sum1,sum2,theta,delth,zcord,rcord,gradsq
      REAL*8 dpsidx,dpsidz,gsval,psval,psixz,psixx,psizz,bx,bz,br,bt
!============
      sum1 = 0._R8
      sum2 = 0._R8
      theta = 0._R8
      nn = 100
      delth = 2._R8*pi/nn
!
      do 100 n=1,nn
      theta = theta + delth
      zcord =        1.196_R8*sin(theta)
      rcord = 2.65_R8+ 1.196_R8*cos(theta)
      iflip=0
      if(zcord.lt.0) then
                     iflip=1
                     zcord = -zcord
                      endif
      call grap(3,zcord,rcord,gradsq,dpsidx,dpsidz,gsval,psval,          &  
     &          psixz,psixx,psizz,1)
      if(iflip.eq.1) dpsidz = -dpsidz
      bx =  dpsidz/rcord
      bz = -dpsidx/rcord
      br = cos(theta)*bx + sin(theta)*bz
      bt =-sin(theta)*bx + cos(theta)*bz
      sum1 = sum1 + sin(theta)*br
      sum2 = sum2 + cos(theta)*bt
!
  100 continue
!
      bth = sum2/nn
      brho = sum1/nn
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
