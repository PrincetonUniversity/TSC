!     ------------------------------------------------------------------
      SUBROUTINE VolCalc
      USE params
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ugrid
      INTEGER ips
      REAL*8 xlookup, yreturn
      ips = Npsi
      iVlAry(ips) = iVlVec(NpsiJ)
!
!     compute integral of volume by interpolation
      do 10 ips = 1, Npsi-1
        xlookup = MidAry(ips)
        call linr1d(NpsiJ, MidVec , iVlVec, xlookup, yreturn)
        iVlAry(ips) = yreturn
 10   continue
!
!     compute dVol which is centered
        ips =1
        dVol(ips) = iVlAry(ips)
      do 30 ips = 2, Npsi
        dVol(ips) = iVlAry(ips) - iVlAry(ips-1)
 30   continue
 
!        call LSCpause
!        write(6,'('' psi index, centered volume elements,'',
!     ^            '' integrated vol:'')')
!        do 50 ips=1,npsi
!          write(6,'(1x,i4, 1x, 1pe10.3, 4x,1pe10.3)')
!     ^         ips, dVol(ips), iVlAry(ips)
! 50     continue
!        call LSCpause
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
