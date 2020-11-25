!
!     -----------------------------------------------------------------
      SUBROUTINE RfDamp
      USE DqlBins
      USE FeBins
      USE params
      USE power_mod
      USE ProfBody
      USE RayBins
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER izn, iry, iv, ips
!     uses decrement in ray energy along each ray to compute local
!     deposition and to deposit same into appropriate velocity
!     bin, differential volume element
!     assume power(izone, iray) is the power arriving at psi surface,
!     velocity zone izind(izone, iray), ivind(izone, iray) and use
!     differences to compute local deposition
!
      call BrodCast(NPSIDIM * NVELDIM, PRay,    0._R8)
      call BrodCast(NPSIDIM,           PRaytot, 0._R8)
      do 10 iry = 1, nrays
        do 20 izn = 1, nzones - 1
 
          iv = ivind(izn, iry)
!                                       If  ivind=0, the ray was stopped before
!                                       reaching this izone (izn).
!                                       No calculation is appropriate.
          if (iv .eq. 0) go to 20
          ips = izind(izn, iry)
          PRay(iv, ips) = PRay(iv, ips) + ( power(izn  ,iry) -           &  
     &                                      power(izn+1,iry)  )
 20     continue
 10   continue
!
      call Smooth(Pray, nv, NVELDIM, npsi, nsmoo, qlsm)
!
      do 30 ips = 1, NPSIDIM
        do 31 iv  = 1, NVELDIM
          PRayTot(ips) = PRayTot(ips) + PRay(iv,ips)
 31     continue
 30   continue
!
      PRaySum = 0._R8
      do 40 ips = 1, npsi
        PRaySum = PRaySum + PRaytot(ips)
 40   continue
!
      PPwrSum = 0._R8
      do 50 iry = 1, nrays
        PPwrSum = PPwrSum + ( power(1,iry) - power(nzones,iry) )
 50   continue
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
