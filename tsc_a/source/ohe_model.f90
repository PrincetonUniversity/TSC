      subroutine ohe_model(                                              &  
     &   chibohm,     zq,         zeps,       zrhostar,   zrlt           &  
     & , zrltcrit,    zra,        chii,       chie,       dh             &  
     & , dimp,        zkappa,     zkapexp,    zbound,     ierr           &  
     & , zcoef1,      zcoef2,     zcoef3,     zcoef4,     switch)
!
! Ottaviani-Horton-Erba Transport Model
!
! Inputs:
!    chibohm:   Bohm diffusivity, defined as T_e/eB, in MKS units.
!         *** all other inputs are dimensionless ***
!    zq:        The local value of q, the safety factor
!    zeps:      The inverse aspect ratio of the local flux surface
!    zrhostar:  The local value of rho-star, defined as the ion gyroradius
!                    rho_i divided by the minor radius, a, of the plasma
!    zrlt:      R/L_T, where R is the major radius of the plasma,
!                    and L_T is the local ion temperature scale length
!    zrltcrit:  The value of R/L_T associated with the local critical
!                    temperature threshold (for ITG transport); used
!                    only used with the OHEm model-- switch=1
!    zra:       The minor radius of the local flux surface, divided by
!                    the minor radius of the plasma (r/a)
!    zkappa:    Elongation of the local flux surface
!    zkapexp:   Power-law exponent on zkappa; this geometric factor
!                    multiplies the thermal diffusivities
!    zbound:    When multiplied by q*rho_s, the lower bound on the
!                    radial correlation length for the turbulence
!    zcoef1:    Coefficient for empirical hydrogen diffusivity
!    zcoef2:    Coefficient for empirical hydrogen diffusivity
!    zcoef3:    Coefficient for empirical impurity diffusivity
!    zcoef4:    Coefficient for empirical impurity diffusivity
!    switch:    Integer switch.  =0 selects the thresholdless OHE model
!                                =1 selects the thresholded OHEm model
!
! Outputs:
!    chii:      The ion thermal diffusivity, in chibohm's units
!    chie:      The electron thermal diffusivity, in chibohm's units
!    dh:        The hydrogenic ion particle diffusivity, in [m^2/sec]
!    dimp:      The impurity ion particle diffusivity, in [m^2/sec]
!    ierr:      If input data incorrect, output value of ierr is returned
!
!
! Declare variables
! NAMING CONVENTION: Dimensionless variables begin with a 'z'
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8                                                             &  
     &   chibohm,     zq,         zeps,       zrhostar,   zrlt           &  
     & , zrltcrit,    zra,        chii,       chie,       dh             &  
     & , dimp,        zkappa,     zkapexp,    zbound,     zfkappa        &  
     & , zalpha,      zrltb,      zchii,      zcoef1,     zcoef2         &  
     & , zcoef3,      zcoef4
 
      INTEGER switch, ierr
!
! check input for validity
!
      ierr = 0
      if ((zq .lt. 0.0_R8) .or. (zq .gt. 100.0_R8)) then
         ierr=1
         return
      elseif ((zeps .lt. 0.0_R8) .or. (zeps .gt. 1.0_R8)) then
         ierr=2
         return
      elseif ((zrhostar .lt. 0.0_R8) .or. (zrhostar .gt. 1.0_R8)) then
         ierr=3
         return
      elseif ((zbound .lt. 0.0_R8).or. (zbound .gt. 100.0_R8)) then
         ierr=4
         return
      elseif ((zra .le. 0.0_R8) .or. (zra .gt. 1.0_R8)) then
         ierr=5
         return
      elseif ((zkappa .lt. 0.0_R8) .or. (zkappa .gt. 10.0_R8)) then
         ierr=6
         return
      elseif (abs(zkapexp) .gt. 10.0_R8) then
         ierr=7
         return
      elseif ((zrltcrit .lt. 0.0_R8) .or. (zrltcrit .gt. 100.0_R8))      &  
     & then
         ierr=8
         return
      elseif (abs(zrlt) .gt. 100.0_R8) then
         ierr=9
         return
      endif
      zfkappa = zkappa**zkapexp
      zrlt = abs(zrlt)
! *
! * The Thresholdless Ottaviani-Horton-Erba ITG/TEM model
! *
      if ( switch .eq. 0) then
!
! Calculate dimensionless zchii, with lower bound on lambda_c
!
         zrltb = max( zbound, zrlt )
!
         zchii = zq * zq * (zrlt*zeps*zrhostar/zra) * zrltb**1.5_R8
! *
! * The Thresholded Ottaviani-Horton-Erba ITG/TEM model (OHEm)
! *
      else
!
! Calculate dimensionless zchii, with lower bound on lambda_c
!
         zrltb = max( zbound, zrlt )
!
         zalpha = max( 0.0_R8, zrltb-zrltcrit )
!
         zchii = zq * zq * (zrlt*zeps*zrhostar/zra) * zalpha
      end if
!
! Now determine the actual thermal and particle diffusivities.  The factor
! that depends on plasma elongation, zfkappa, was not included in the
! original OHE model.
!
         chii = (0.014_R8) * chibohm * zchii * zfkappa
         chie = chii * sqrt(zeps)
!
! The hydrogen and impurity diffusivities below are not included in the OHE
! model described in Ref. [1] but were included (with zcoef1=1.0, zcoef2=0.25
! and zcoef3=zcoef4=1.25) in obtaining the results in Ref. [2].
!
         dh   = zcoef1 + zcoef2 * zra * zra
         dimp = zcoef3 + zcoef4 * zra * zra
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
