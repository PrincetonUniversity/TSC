!
!     ------------------------------------------------------------------
!
      SUBROUTINE JsplCalc
      USE Jrf
      USE params
      USE ProfBody
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip,iFillJray, iGotRuna
      REAL*8    jd
      DATA    iFillJray /0/
      REAL*8    ZERO
      DATA    ZERO/                                                      &  
     &         0.0_R8/
!     jd      current driven with the E field incremented by dEdcAmnt
!     compute rf driven current given power dissipation vs vpar, psi
!
 
      call GetEdc(dEdcAmnt)
!     .                                 Fill array EdcAry with dEdcAmnt over
!     .                                 the actual Edc supplied from TSC
!     .                                 dEdcAmnt is assigned in block data
      do 10 ip = 1, npsi
        call jnorm(ip)
!     .                                 set up normalization,
!     .                                 tabulate u lookup values;
!     .                                 These depend on Ez
        call mkj(ip,jd, iGotRuna,iFillJray )
!     .                                 compute jrf driven, jd
        jsp(ip) = jd
 10   continue
      call GetEdc(ZERO)
!     call GetEdc(0.0)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
