      MODULE emitter
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
!     emitter.inc -------------------------------------------------------------
      INTEGER nr_source, nz_source
      REAL*8    mu_0, mu_width,                                          &  
     &        Z_bound_min, Z_bound_max, R_bound_min, R_bound_max,        &  
     &        Z_plasm_min, Z_plasm_max, R_plasm_min, R_plasm_max,        &  
     &        r_source(NRDIM), z_source(NZDIM),                          &  
     &        R_bound_min_sq, R_bound_max_sq
      REAL*8    PusherMajor, PusherMinor
      REAL*8    source_profile(NRDIM, NZDIM)

!cj   DATA Z_bound_min, Z_bound_max /                                    &  
!cj  &           -1.00_R8,        1.00_R8/
!cj   DATA R_bound_min, R_bound_max /                                    &  
!cj  &            1.10_R8,        2.00_R8/
!cj   DATA mu_0 / 1._R8/ mu_width / .2_R8/
!cj   DATA nr_source / NRDIM / nz_source / NZDIM /
!cj   DATA  PusherMajor, PusherMinor / 1.152_R8, 0.190_R8/



 
!     COMMON / emcom1 / nr_source, nz_source
!     COMMON / emcom2 / mu_0, mu_width,
!    ^        Z_bound_min, Z_bound_max, R_bound_min, R_bound_max,
!    ^        Z_plasm_min, Z_plasm_max, R_plasm_min, R_plasm_max,
!    ^     r_source, z_source, R_bound_min_sq,
!    ^     R_bound_max_sq, PusherMajor, PusherMinor
!     COMMON / emcom3 / source_profile
!     emitter.inc -------------------------------------------------------------
!
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE emitter
