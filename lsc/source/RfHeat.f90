!     -----------------------------------------------------------------
      SUBROUTINE RfHeat
      USE DqlBins
      USE FeBins
      USE params
      USE power_mod
      USE ProfBody
      USE RayBins
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ips, iv, iSMOi
!                                       compute power flow resulting
!                                       from quasilinear diffusion
      call BrodCast(NPSIDIM,           PQltot, 0._R8)
      call BrodCast(NVELDIM * NPSIDIM, Pql,    0._R8)
      iSMOi = mod(iSMO,2) + 1
!
!     The quasilinear power deposited in each velocity bin in each
!     psi bin is:
!     Pql = m v Dql (-df/dv) dv dVol/dpsi dpsi
      do 20 ips = 1, npsi
              do 10 iv = 1, nv
                 Pql(iv, ips) =  - Vpar(iv) *                            &  
     &                       Dql(iv, ips,iSMO) * dfdv(iv, ips,iSMOi) *   &  
     &                       dvsym(iv) * dVol(ips) * PwrNorm
                 Pqltot(ips) = Pqltot(ips) + Pql(iv, ips)
 10           continue
 20   continue
!
      call Smooth(Pql , nv, NVELDIM, npsi, nsmoo, qlsm)
      PqlSum = 0._R8
      do 30 ips = 1, npsi
         PqlSum = PqlSum + Pqltot(ips)
 30   continue
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
