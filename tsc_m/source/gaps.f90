      subroutine gaps(g1,g2,g3,g4,g5,g6,dg1,dg2,dg3,dg4,dg5,dg6)
!
!C        subroutine gaps(g1,g2,g3,g4,g5,g6)
!
!..   This subroutine calculates the gaps to be used as input
!      for the shape control algorithm from the CREATE group.
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER lgapcnt,imode,isw
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 g2,g3,g4,g5,g6,dg1,dg2,dg3,dg4,dg5,dg6,g1,xg1,zg1,vx
      REAL*8 vz,grsq,dpsidx,dpsidz,gsval,psval,psixz,psixx,psizz
      REAL*8 grpsi,xg2,zg2,xg3,zg3,xg4,zg4,xg5,zg5,xg6,zg6,g1m2
      REAL*8 g2m2,g3m2,g4m2,g5m2,g6m2,g1m1,g2m1,g3m1,g4m1,g5m1,g6m1
      REAL*8 zzip
!============
      if (kcycle .eq. 1) then
!C         open(61,file='gapsout',status='unknown')
         lgapcnt = 49
      end if
!
      pi = 3.14159265358979_R8
!
      imode = 2
      isw =1
 
!.. g1 = the gap between the plasma separtrix and the inner divertor
!    strike point.  @ (5.00,-4.89)
!
      xg1 = 5.00_R8
      zg1 = -4.89_R8
!
!CC      xd1 = xg1 - xsep(1)
!CC      zd1 = zg1 - zsep(1)
!
!CC      g1 = sqrt(xd1*xd1 + zd1*zd1)
!
      vx = 0.526_R8
      vz = -.851_R8
!
      call grap(imode,zg1,xg1,grsq,dpsidx,dpsidz,gsval,psval,            &  
     &      psixz,psixx,psizz,isw)
!
      grpsi = dpsidx * vx + dpsidz * vz
!
      g1 = 1.0_R8*((psilim - psval) / grpsi)
!
!.. g2 = the gap between the plasma separtrix and the outer divertor
!    strike point.  @ (8.57, -5.94)
!
      xg2 = 8.57_R8
      zg2 = -5.94_R8
!
!CC      xd2 = xg2 - xsep(1)
!CC      zd2 = zg2 - zsep(1)
!
!CC      g2 = sqrt(xd2*xd2 + zd2*zd2)
!
      vx = -.857_R8
      vz = -.516_R8
!
      call grap(imode,zg2,xg2,grsq,dpsidx,dpsidz,gsval,psval,            &  
     &      psixz,psixx,psizz,isw)
!
      grpsi = dpsidx * vx + dpsidz * vz
!
      g2 = 1.0_R8*((psilim - psval) / grpsi)
!
!.. g3 = the gap between the separatrix and the ICRF antenna
!    @ (11.25, 1.26).  This gap is measured on the horizontal
!    plane, so it is calculated by (psilim - psi(x,z)) / dpsidx(x,z)
!
      vx = -.993_R8
      vz = 0.121_R8
!
      xg3 = 11.25_R8
      zg3 = 1.26_R8
!
      call grap(imode,zg3,xg3,grsq,dpsidx,dpsidz,gsval,psval,            &  
     &      psixz,psixx,psizz,isw)
!
      grpsi = dpsidx * vx + dpsidz * vz
!
      g3 = 1.0_R8*((psilim - psval) / grpsi)
!
!.. g4 = the gap between the scrape-off layer (SOL) and the FW in
!       the high ripple region  @ (10.61, 3.99)
!
      vx = -.845_R8
      vz = -.535_R8
!
      xg4 = 10.61_R8
      zg4 = 3.99_R8
!
      call grap(imode,zg4,xg4,grsq,dpsidx,dpsidz,gsval,psval,            &  
     &      psixz,psixx,psizz,isw)
!
      grpsi = dpsidx * vx + dpsidz * vz
!
      g4 = 1.0_R8*((psilim - psval) / grpsi)
!
!.. g5 = the gap between the SOL and FW at the top inner corner of
!       the FW  @ (6.69,6.18)
!
      vx = 0.085_R8
      vz = -.996_R8
!
      xg5 = 6.69_R8
      zg5 = 6.18_R8
!
      call grap(imode,zg5,xg5,grsq,dpsidx,dpsidz,gsval,psval,            &  
     &      psixz,psixx,psizz,isw)
!
      grpsi = dpsidx * vx + dpsidz * vz
!
      g5 = 1.0_R8*((psilim - psval) / grpsi)
!
!.. g6 = the gap between the SOL and FW at the inner-most point of
!       the plasma boundary @ (4.91,1.70)
!
      vx = 1.00_R8
      vz = 0.00_R8
!
      xg6 = 4.91_R8
      zg6 = 1.70_R8
!
      call grap(imode,zg6,xg6,grsq,dpsidx,dpsidz,gsval,psval,            &  
     &      psixz,psixx,psizz,isw)
!
      grpsi = dpsidx * vx + dpsidz * vz
!
      g6 = 1.0_R8*((psilim - psval) / grpsi)
!
!.... Calculate the time derivative of the gaps....
!
!*******Make sure to use correct value of time step dt = dts
!
      dg1 = (g1 - g1m2) / (2._R8*dts)
      dg2 = (g2 - g2m2) / (2._R8*dts)
      dg3 = (g3 - g3m2) / (2._R8*dts)
      dg4 = (g4 - g4m2) / (2._R8*dts)
      dg5 = (g5 - g5m2) / (2._R8*dts)
      dg6 = (g6 - g6m2) / (2._R8*dts)
!
!.... move kcycle-1 value of gap to kcycle-2 value
!
      g1m2 = g1m1
      g2m2 = g2m1
      g3m2 = g3m1
      g4m2 = g4m1
      g5m2 = g5m1
      g6m2 = g6m1
!
!.... move present value of gaps to kcycle-1 value for next time step
!
      g1m1 = g1
      g2m1 = g2
      g3m1 = g3
      g4m1 = g4
      g5m1 = g5
      g6m1 = g6
!
      if (kcycle.lt.3) then
         dg1 = 0._R8
         dg2 = 0._R8
         dg3 = 0._R8
         dg4 = 0._R8
         dg5 = 0._R8
         dg6 = 0._R8
       end if
!
!
!.... plasma current in MA = tpi*tcurdtp*udsi*1.e-6
!
        zzip =  tpi * tcurdtp * udsi * 1.E-6_R8
!
       lgapcnt = lgapcnt + 1
!
       if (lgapcnt .eq. 50) then
          write(61,1001) times, g1, g2, g3, g4, g5, g6, zzip
          lgapcnt = 0
       end if
!
 1001  format(1x,e15.8,1x,e15.8,1x,e15.8,1x,e15.8,1x,e15.8,1x,e15.8,     &  
     &      1x,e15.8,1x,e15.8)
!
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
