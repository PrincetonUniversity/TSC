      subroutine gaps3(g1,g2,g3,g4,g5,g6)
!
!...this subroutine calculates the gaps to be used as input
!...for the shape control algorithm from the GA/LLNL group

      USE CLINAM
      USE ITER1
      
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER imode,isw,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 g2,g3,g4,g5,g6,dg1,dg2,dg3,dg4,dg5,dg6,g1,xp1,zp1,grsq
      REAL*8 dpsidx,dpsidz,gsval,psval,psixz,psixx,psizz,xp2,zp2
      REAL*8 xp3,zp3,xp4,zp4,xp5,zp5,xp6,zp6
!============
      REAL*8, DIMENSION(ntpts) :: xcon,zcon
      
      pi = 3.14159265358979_R8
      imode = 2
      isw =1
      do 1000 i=1,ntpts
      if(tpro(i) .le. time .and. tpro(i+1) .gt. time) then
      do 1001 j=1,ncnt
      xcon(j) = xcon0(i,j) + ((xcon0(i+1,j)-xcon0(i,j))/         &
     & (tpro(i+1)-tpro(i)))*(tpro(i+1)-time)
      zcon(j) = zcon0(i,j) + ((zcon0(i+1,j)-zcon0(i,j))/         &
     & (tpro(i+1)-tpro(i)))*(tpro(i+1)-time)
 1001 continue
      go to 1002
      endif
 1000 continue
 1002 continue
!
!...g1 = the gap between the plasma separatrix and the inner divertor
!...strike point    wall(5.325,-4.735)   plasma(4.991,-5.282)
!
      xp1 = xcon(1)
      zp1 = zcon(1)
!     xp1 = 4.200
!     zp1 = -3.800
      call grap(imode,zp1,xp1,grsq,dpsidx,dpsidz,gsval,psval, &
     &      psixz,psixx,psizz,isw)
      g1 = (psval-psisep)/max(sqrt(grsq),1.0E-8_R8)
!
!...g2 = the gap between the plasma separtrix and the outer divertor
!...strike point    wall(8.144,-5.285)   plasma(8.287,-6.009)
! 
      xp2 = xcon(2)
      zp2 = zcon(2)
!     xp2 = 5.600
!     zp2 = -4.500
      call grap(imode,zp2,xp2,grsq,dpsidx,dpsidz,gsval,psval, &
     &      psixz,psixx,psizz,isw)
      g2 = (psval-psisep)/max(sqrt(grsq),1.0E-8_R8)
! 
!...g3 = the gap between the separatrix and the ICRF antenna
!...   wall(11.245,1.922)   plasma(10.565,-0.374)
! 
      xp3 = xcon(3)
      zp3 = zcon(3)
!     xp3 = 8.200
!     zp3 = 0.420
      call grap(imode,zp3,xp3,grsq,dpsidx,dpsidz,gsval,psval, &
     &      psixz,psixx,psizz,isw)
      g3 = (psval-psisep)/max(sqrt(grsq),1.0E-8_R8)
! 
!...g4 = the gap between the scrape-off layer (SOL) and the FW in
!...the high ripple region   wall(10.134,4.559)   plasma(10.070,4.305)
! 
      xp4 = xcon(4)
      zp4 = zcon(4)
!     xp4 = 7.170
!     zp4 = 3.050
      call grap(imode,zp4,xp4,grsq,dpsidx,dpsidz,gsval,psval, &
     &      psixz,psixx,psizz,isw)
      g4 = (psval-psisep)/max(sqrt(grsq),1.0E-8_R8)
! 
!...g5 = the gap between the SOL and FW at the top inner corner of
!...the FW   wall(7.618,6.085)   plasma(6.930,5.889)
! 
      xp5 = xcon(5)
      zp5 = zcon(5)
!     xp5 = 5.420
!     zp5 = 4.050
      call grap(imode,zp5,xp5,grsq,dpsidx,dpsidz,gsval,psval, &
     &      psixz,psixx,psizz,isw)
      g5 = (psval-psisep)/max(sqrt(grsq),1.0E-8_R8)
! 
!...g6 = the gap between the SOL and FW at the inner-most point of
!...the plasma boundary   wall(4.927,1.825)   plasma(5.166,1.973)
! 
      xp6 = xcon(6)
      zp6 = zcon(6)
!     xp6 = 4.200
!     zp6 = 1.000
      call grap(imode,zp6,xp6,grsq,dpsidx,dpsidz,gsval,psval, &
     &      psixz,psixx,psizz,isw)
      g6 = (psval-psisep)/max(sqrt(grsq),1.0E-8_R8)
! 
!.... calculate the time derivative of the gaps
! 
      dg1 = (g1 - g1old) / dts
      dg2 = (g2 - g2old) / dts
      dg3 = (g3 - g3old) / dts
      dg4 = (g4 - g4old) / dts
      dg5 = (g5 - g5old) / dts
      dg6 = (g6 - g6old) / dts
      g1old = g1
      g2old = g2
      g3old = g3
      g4old = g4
      g5old = g5
      g6old = g6
      if (kcycle.lt.2) then
      dg1 = 0.
      dg2 = 0.
      dg3 = 0.
      dg4 = 0.
      dg5 = 0.
      dg6 = 0.
      endif
! 
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
