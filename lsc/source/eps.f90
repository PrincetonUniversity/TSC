!
!     eps -- dielectric tensor              ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE eps ( r, z, kpar2, kper2 )
!============
! idecl:  explicitize implicit REAL declarations:
      USE dielec
      USE MKSetc
      USE params
      USE PIetc
      USE RayBins
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
      EXTERNAL plasma2d, bsi0, bsi1
      REAL*8                                                             &  
     &    r,z,kpar2,kper2,                                               &  
     &    HUGE,veow2,elam,elamc,emli0,emli1,                             &  
     &    exy0,bsi0,bsi1,                                                &  
     &                        psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
 
      DATA HUGE / 1.0E+30_R8/
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      If ( psi .ge. HUGE ) return
      ecyc  = omc/fghz
      ecyc2 = ecyc*ecyc
      epsq  = pe2/fghz**2
      ipsq  = pi2/fghz**2
      veow2 = 0.0445E-04_R8*tee/fghz**2
      elamc = veow2/ecyc2
      elam  = elamc * Kper2
!     This is the electron Lambda - - (Kper RhoE)^2
!     Now form the quantities exp(-elam) I0,1 (elam)
      emli0 = bsi0(elam)
      emli1 = bsi1(elam)
      epar  = 1._R8- epsq
      aion = aio/fghz**4
      aelc = ael
      eper = 1._R8+ epsq/ecyc2  - ipsq
      eper = eper - kper2 * ( aion + aelc )
      exy0  = epsq/ecyc
      exy   = exy0
      return
 
      ENTRY epsdr  ( r, z, Kpar2, Kper2 )
!     derivatives of the dielectric tensor elements
!     by (k-perp)**2 ; (k-par)**2 ; (omega)**2
!     this last is multiplied by (omega**2)
 
      d11er = -(aion + aelc)
      d11ar = 0._R8
      d11w0 = ipsq + 2._R8*aion*kper2
 
      d12er = 0._R8
      d12ar = 0._R8
      d12w0 = -0.5_R8*exy
 
      d33er = 0._R8
      d33ar = 0._R8
      d33w0 = epsq
      return
 
      END
!                                                                      |
!                                                                      |
!     eps ends                              ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
