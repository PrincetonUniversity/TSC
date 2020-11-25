      subroutine volt(v1,v2,v3,v4,v5,v6,v7)
!
!...  This routine calculates the voltages on the control coils
!...  based on the input of the gaps, time derivative of the
!...   gaps and the plasma current using the CREATE gain matrix.
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER lgapcnt,i,j,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 v2,v3,v4,v5,v6,v7,v1,vdwgain,gapdw,dgapdw,brhs,gap0
      REAL*8 zaip0
!============
      dimension vdwgain(7,13), gapdw(6), dgapdw(6), brhs(13)
      dimension gap0(6)
!============      
!
!... skip definition of gain matrix after code initiation
!
      if (kcycle .eq. 1) then
      if( numargs .lt. 1 ) then
         filename =  'gapsout'
      else
         filename =  'gapsout' // '.' // trim(suffix)
      end if
         open(61,file=trim(filename),status='unknown')
         write(61,1002)
         lgapcnt = 49
      end if
!
 1002  format(1x,'*** Volt ****')
!
!C... gaps nominal value definition....
!
!C       gap0(1) = 0.227
!C       gap0(2) = 0.227
!C       gap0(3) = 0.213
!C       gap0(4) = 0.203
!C       gap0(5) = 0.522
!C       gap0(6) = 0.155
!C       zaip0 = 24.00
!
!C..  test nominal values from TSC equil...
        gap0(1) = 0.240_R8
        gap0(2) = 0.180_R8
        gap0(3) = 0.215_R8
        gap0(4) = 0.200_R8
        gap0(5) = 0.586_R8
        gap0(6) = 0.190_R8
        zaip0 = 24.00_R8
!
!CC   if (kcycle .gt. 100) goto 51
!CC      if (kcycle .eq. 0)  lgapcnt = 49
!
      vdwgain(1,1) = 4.8755_R8
      vdwgain(1,2) = -1.2944_R8
      vdwgain(1,3) = -0.7191_R8
      vdwgain(1,4) = -0.9203_R8
      vdwgain(1,5) = -1.4753_R8
      vdwgain(1,6) = 2.6513_R8
      vdwgain(1,7) = 0.0222_R8
      vdwgain(1,8) = 0.1094_R8
      vdwgain(1,9) = -0.0348_R8
      vdwgain(1,10) = -0.0875_R8
      vdwgain(1,11) = -0.1353_R8
      vdwgain(1,12) = -0.0348_R8
      vdwgain(1,13) = 0.0843_R8
!
      vdwgain(2,1) = 0.4483_R8
      vdwgain(2,2) = 0.3270_R8
      vdwgain(2,3) = -2.5008_R8
      vdwgain(2,4) = 2.1002_R8
      vdwgain(2,5) = 3.4840_R8
      vdwgain(2,6) = 0.5160_R8
      vdwgain(2,7) = -.0306_R8
      vdwgain(2,8) = -0.0107_R8
      vdwgain(2,9) = -0.0104_R8
      vdwgain(2,10) = -0.2863_R8
      vdwgain(2,11) = 0.1314_R8
      vdwgain(2,12) = 0.1990_R8
      vdwgain(2,13) = 0.0208_R8
!
      vdwgain(3,1) = 1.3958_R8
      vdwgain(3,2) = 0.3502_R8
      vdwgain(3,3) = -1.6612_R8
      vdwgain(3,4) = 3.1105_R8
      vdwgain(3,5) = 0.8397_R8
      vdwgain(3,6) = -0.2545_R8
      vdwgain(3,7) = -0.0278_R8
      vdwgain(3,8) = 0.0027_R8
      vdwgain(3,9) = 0.0038_R8
      vdwgain(3,10) = -0.1405_R8
      vdwgain(3,11) = 0.3387_R8
      vdwgain(3,12) = 0.0419_R8
      vdwgain(3,13) = -0.0294_R8
!
      vdwgain(4,1) = 0.7163_R8
      vdwgain(4,2) = 0.3301_R8
      vdwgain(4,3) = 0.6813_R8
      vdwgain(4,4) = 0.4142_R8
      vdwgain(4,5) = 0.4425_R8
      vdwgain(4,6) = -0.1330_R8
      vdwgain(4,7) = 0.0152_R8
      vdwgain(4,8) = -0.0015_R8
      vdwgain(4,9) = 0.0087_R8
      vdwgain(4,10) = 0.1213_R8
      vdwgain(4,11) = 0.0334_R8
      vdwgain(4,12) = 0.0031_R8
      vdwgain(4,13) = -0.0079_R8
!
      vdwgain(5,1) = 0.7472_R8
      vdwgain(5,2) = 1.5808_R8
      vdwgain(5,3) = 0.2244_R8
      vdwgain(5,4) = 0.2786_R8
      vdwgain(5,5) = 0.2500_R8
      vdwgain(5,6) = -0.4218_R8
      vdwgain(5,7) = 0.0295_R8
      vdwgain(5,8) = 0.0006_R8
      vdwgain(5,9) = 0.0456_R8
      vdwgain(5,10) = 0.0628_R8
      vdwgain(5,11) = 0.0279_R8
      vdwgain(5,12) = -0.0015_R8
      vdwgain(5,13) = -0.0279_R8
!
      vdwgain(6,1) = 0.2254_R8
      vdwgain(6,2) = 3.4867_R8
      vdwgain(6,3) = -1.1662_R8
      vdwgain(6,4) = -0.7877_R8
      vdwgain(6,5) = -0.8265_R8
      vdwgain(6,6) = 1.3806_R8
      vdwgain(6,7) = 0.0074_R8
      vdwgain(6,8) = 0.0123_R8
      vdwgain(6,9) = 0.0941_R8
      vdwgain(6,10) = -0.1809_R8
      vdwgain(6,11) = -0.0970_R8
      vdwgain(6,12) = -0.0068_R8
      vdwgain(6,13) = 0.0437_R8
!
      vdwgain(7,1) = 0.8504_R8
      vdwgain(7,2) = 1.4542_R8
      vdwgain(7,3) = 1.4316_R8
      vdwgain(7,4) = -0.7587_R8
      vdwgain(7,5) = -2.6042_R8
      vdwgain(7,6) = -1.1999_R8
      vdwgain(7,7) = 0.0225_R8
      vdwgain(7,8) = 0.0311_R8
      vdwgain(7,9) = 0.0605_R8
      vdwgain(7,10) = 0.1753_R8
      vdwgain(7,11) = 0.0291_R8
      vdwgain(7,12) = -0.1442_R8
      vdwgain(7,13) = -0.0930_R8
!
!...  multiply all gains by 10^3 for correct normalization...
!
      do 11 i=1,7
         do 12 j = 1,13
           vdwgain(i,j) = 1000._R8* vdwgain(i,j)
!C           vdwgain(i,j) = 100. * vdwgain(i,j)
 12      continue
 11   continue
!
!.... start here after the gain matrix has been defined
 51   continue
!
      call gaps(gapdw(1),gapdw(2),gapdw(3),gapdw(4),gapdw(5),            &  
     & gapdw(6),dgapdw(1),dgapdw(2),dgapdw(3),dgapdw(4),dgapdw(5),       &  
     &  dgapdw(6))
!
      do 100 i = 1,6
!
         ii = i + 7
!
         brhs(i) = (gapdw(i) - gap0(i))
         brhs(ii) = dgapdw(i)
!
 100  continue
!
!.... plasma current in MA = tpi*tcurdtp*udsi*1.e-6
!
      brhs(7) = (tpi * tcurdtp * udsi * 1.E-6_R8) - zaip0
!
      v1 = 0._R8
      v2 = 0._R8
      v3 = 0._R8
      v4 = 0._R8
      v5 = 0._R8
      v6 = 0._R8
      v7 = 0._R8
!
      do 200 i = 1,13
!
         v1 = v1 + vdwgain(1,i) * brhs(i)
         v2 = v2 + vdwgain(2,i) * brhs(i)
         v3 = v3 + vdwgain(3,i) * brhs(i)
         v4 = v4 + vdwgain(4,i) * brhs(i)
         v5 = v5 + vdwgain(5,i) * brhs(i)
         v6 = v6 + vdwgain(6,i) * brhs(i)
         v7 = v7 + vdwgain(7,i) * brhs(i)
!
 200  continue
!
      if (v1 .gt. 2.0_R8) v1 = 2.0_R8
      if (v1 .lt. -2.0_R8) v1 = -2.0_R8
      if (v2 .gt. 10.0_R8) v2 = 10.0_R8
      if (v2 .lt. -10.0_R8) v2 = -10.0_R8
      if (v3 .gt. 10.0_R8) v3 = 10.0_R8
      if (v3 .lt. -10.0_R8) v3 = -10.0_R8
      if (v4 .gt. 10.0_R8) v4 = 10.0_R8
      if (v4 .lt. -10.0_R8) v4 = -10.0_R8
      if (v5 .gt. 10.0_R8) v5 = 10.0_R8
      if (v5 .lt. -10.0_R8) v5 = -10.0_R8
      if (v6 .gt. 10.0_R8) v6 = 10.0_R8
      if (v6 .lt. -10.0_R8) v6 = -10.0_R8
      if (v7 .gt. 10.0_R8) v7 = 10.0_R8
      if (v7 .lt. -10.0_R8) v7 = -10.0_R8
!
!C1turn      v1 = v1 * 29.5716
!
       lgapcnt = lgapcnt + 1
!
       if (lgapcnt .eq. 50) then
          write(61,1001) times, v1, v2, v3, v4, v5, v6, v7
          lgapcnt = 0
       end if
!
!C          print*, 'V =', times, v1, v2, v3, v4, v5, v6, v7
!
 1001  format('***',1x,e15.8,1x,e15.8,1x,e15.8,1x,e15.8,1x,e15.8,1x,     &  
     &      e15.8,1x,e15.8,1x,e15.8)
!
!
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
