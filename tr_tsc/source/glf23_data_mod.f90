!c@glf.m 11-Apr-01 J. Kinsey, General Atomics
! 03mar2007 pankin converted into F90 module
! 05-mar-01 changed 20 to nmode
! 23-aug-00 aligned common block, added xky_gf
! 14-june-00 added ngrow_k_gf
! 13-june-00 added ipert_gf
! 03-aug-99 added jeigen
!---------------------------------------------------------------------
    MODULE glf23_data_mod
!*FD This module contains data structure necessary for glf.f
        use spec_kind_mod

          IMPLICIT NONE
          INTEGER nmode
          PARAMETER (nmode=20)
          INTEGER ::iflagin_gf(30)=0, ngrow_k_gf(0:nmode)=0
          INTEGER :: nroot_gf=0
          INTEGER :: jeigen=0 &
                  ,lprint_gf=0,ikymax_gf=0,eigen_gf=0 &
                  ,i_err=0,first_order_gf=0,ipert_gf=0
!***NOTE: necessary for have quantities with no spaces and comma leading the
!line for f90doc.py to work.
!
!        DOUBLE PRECISION :: yparam_k_gf(nmode,nmode)=0. &
          REAL(r8) :: yparam_k_gf(nmode,nmode)=0. &
          ,gamma_k_gf(1:4,nmode)=0.,freq_k_gf(1:4,nmode)=0. &
          ,phi_norm_k_gf(1:4,nmode)=0. &
          ,xparam_gf(30)=0. &
          ,xkyf_k_gf(nmode)=0.,diff_k_gf(nmode)=0. &
          ,diff_im_k_gf(nmode)=0.,chii_k_gf(nmode)=0. &
          ,chie_k_gf(nmode)=0.,exch_k_gf(nmode)=0. &
          ,eta_par_k_gf(nmode)=0.,eta_per_k_gf(nmode)=0.,eta_phi_k_gf(nmode)=0. &
          ,chie_e_k_gf(nmode)=0., yparam_gf(nmode)=0. &
          ,gamma_gf(1:4)=0.,freq_gf(1:4)=0.,phi_norm_gf(1:4)=0.,xky_gf(1:4)=0. &
          ,xky0_gf=0.,rms_theta_gf=0.,rlti_gf=0. &
          ,rlte_gf=0.,rlne_gf=0.,rlni_gf=0.,rlnimp_gf=0.,dil_gf=0.,apwt_gf=0. &
          ,aiwt_gf=0.,taui_gf=0.,rmin_gf=0.,rmaj_gf=0.,q_gf=0.,xnu_gf=0.,betae_gf=1.e-6 &
          ,shat_gf=0.,alpha_gf=0.,elong_gf=0.,xwell_gf=0.,park_gf=0.,ghat_gf=0. &
          ,gchat_gf=0.,adamp_gf=0.,alpha_star_gf=0.,gamma_star_gf=0. &
          ,alpha_e_gf=0.,gamma_e_gf=0.,alpha_mode_gf=0.,gamma_mode_gf=0. &
          ,alpha_p_gf=0.,gamma_p_gf=0.,xkdamp_gf=0.,xkyf_gf=0. &
          ,diff_gf=0.,diff_im_gf=0.,chii_gf=0.,chie_gf=0.,exch_gf=0.,eta_par_gf=0. &
          ,eta_per_gf=0.,eta_phi_gf=0.,chie_e_gf=0. &
          ,cnorm_gf=0.,xkymin_gf=0.,xkymax_gf=0.,amassgas_gf=0.,amassimp_gf=0.,zimp_gf=0. 

!      double complex :: zevec_k_gf(nmode,12,12)=0. &
          COMPLEX(r8) :: zevec_k_gf(nmode,12,12)=0. &
          ,zomega_k_gf(nmode,12)=0.

        END MODULE glf23_data_mod
