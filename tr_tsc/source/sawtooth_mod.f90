      MODULE sawtooth_mod
      USE kind_spec_mod

      IMPLICIT NONE
      LOGICAL :: saw_on = .false., saw_tr = .false.
      CHARACTER*128 :: saw0_filename="saw0_state.cdf"
      CHARACTER*128 :: saw1_filename="saw1_state.cdf"
      REAL(kind=rspec) :: sawtooth_time
      REAL(kind=rspec) :: saw_t0=-1.0, saw_t1=-1.0

      CONTAINS

      SUBROUTINE tr_saw_on
      USE plasma_state_mod
      USE trtsc
      USE clinam, only : ineg

      IMPLICIT none
      INTEGER :: istat !, jsystem
      CHARACTER*128 :: extension

#ifdef PARTIAL_PLASMA_STATE
      integer :: cclist(ps_ccount)
      CHARACTER*128 :: saw0_partial_ps_filename="saw0_update_state.cdf"
      CHARACTER*128 :: saw1_partial_ps_filename="saw1_update_state.cdf"
#endif

      call wrgeqdsk
      call wrxpls
      extension="cp "//trim(suffix2)//"_ps.geq"//" "//"saw0_state.geq"
      istat=jsystem(trim(extension))
      extension="cp "//trim(suffix2)//"_ps.geq"//" "//"saw1_state.geq"
      istat=jsystem(trim(extension))
      call put_plasma_state_saw0
      if(ierr .ne. 0) then
      ineg=70
      return
      endif

      call ps_copy_plasma_state(saw0,saw1,ierr)
#ifdef PARTIAL_PLASMA_STATE
      write(*,*) &
      "JCdebug sawtooth_mod: tr_saw_on: partial ps store"

!     call PS_RESTORE_HASH_CODES(ierr, state=saw1)

      call ps_cclist_remove('*',cclist,ierr)
      call ps_cclist_add('PLASMA',cclist,ierr)
!     call ps_cclist_add('EQ',cclist,ierr)
!     call ps_cclist_add('NBI',cclist,ierr)
!     call ps_cclist_add('IC',cclist,ierr)
!     call ps_cclist_add('LH',cclist,ierr)
!     call ps_cclist_add('EC',cclist,ierr)
!     call ps_cclist_add('RAD',cclist,ierr)
!     call ps_cclist_add('GAS',cclist,ierr)

      !  cclist contains the PLASMA component only.
      call ps_update_report(ierr,state=saw1,cclist=cclist,              &
                        report_filename='myupdate.list')
!       call ckerr('ps_update_report (5)')

!     call ps_copy_plasma_state(saw1, aux, ierr, cclist=cclist)
!     call ckerr('ps_copy_plasma_state')

      call ps_write_update_file(trim(saw1_partial_ps_filename), ierr,   &
                                state=saw1, cclist=cclist)
!     call ps_store_plasma_state(ierr, filename=trim(saw1_filename), state=saw1)
!     call ckerr('ps_store_plasma_state') 
!#else
#endif
      call ps_store_plasma_state(ierr,trim(saw1_filename),state=saw1)

      sawtooth_time=saw1%t0
      extension="sawtooth_begin"
      open(51,file="tmpfile",status="unknown",iostat=istat)
      write(51,"(f12.6,2x,A,2x,A)") sawtooth_time, trim(saw1_filename), &
     &                              trim(extension)
      write(51,"(f12.6,2x,A)") ps%t1, trim(xplasma_filename)
      close(51)
      extension="mv tmpfile transp_step.dat"
      istat=jsystem(trim(extension))
      extension=""
      call flush(6)

      extension="saw1"
      call wait_transp(extension)
      saw_tr=.true.
      return

      END SUBROUTINE tr_saw_on

      SUBROUTINE tr_saw_off
      USE plasma_state_mod
      USE trtsc
      USE clinam, only : ineg
      IMPLICIT none
      INTEGER :: istat !, jsystem
      CHARACTER*128 :: extension

#ifdef PARTIAL_PLASMA_STATE
      integer :: cclist(ps_ccount)
      CHARACTER*128 :: saw0_partial_ps_filename="saw0_update_state.cdf"
      CHARACTER*128 :: saw1_partial_ps_filename="saw1_update_state.cdf"
#endif

      call wrgeqdsk
      call wrxpls
      extension="cp "//trim(suffix2)//"_ps.geq"//" "//"saw1_state.geq"
      istat=jsystem(trim(extension))
      call put_plasma_state_saw1
      if(ierr .ne. 0) then 
      ineg=70
      return  
      endif   

#ifdef PARTIAL_PLASMA_STATE
      write(*,*) &
      "JCdebug sawtooth_mod: tr_saw_off: partial ps store"

!     call PS_RESTORE_HASH_CODES(ierr, state=saw1)

      call ps_cclist_remove('*',cclist,ierr)
      call ps_cclist_add('PLASMA',cclist,ierr)
!     call ps_cclist_add('EQ',cclist,ierr)
!     call ps_cclist_add('NBI',cclist,ierr)
!     call ps_cclist_add('IC',cclist,ierr)
!     call ps_cclist_add('LH',cclist,ierr)
!     call ps_cclist_add('EC',cclist,ierr)
!     call ps_cclist_add('RAD',cclist,ierr)
!     call ps_cclist_add('GAS',cclist,ierr)

      !  cclist contains the PLASMA component only.
      call ps_update_report(ierr,state=saw1,cclist=cclist,              &
                        report_filename='myupdate.list')
!       call ckerr('ps_update_report (6)')

!     call ps_copy_plasma_state(saw1, aux, ierr, cclist=cclist)
!     call ckerr('ps_copy_plasma_state')

      call ps_write_update_file(trim(saw1_partial_ps_filename), ierr,   &
                                state=saw1, cclist=cclist)
!     call ps_store_plasma_state(ierr, filename=trim(saw1_filename), state=saw1)
!     call ckerr('ps_store_plasma_state') 
!#else
#endif
      call ps_store_plasma_state(ierr,trim(saw1_filename),state=saw1)
   
      return
      END SUBROUTINE tr_saw_off
    

      SUBROUTINE put_plasma_state_saw0

      USE plasma_state_mod
      USE tsc_ps_mod, ONLY : pa, eq, imp
      USE trtsc
      USE EZspline_obj 
      USE EZspline

      IMPLICIT NONE

      INTEGER :: istat
      INTEGER :: i, j, k, l, m, ii, nout=6
      LOGICAL :: ex, l_first=.true.
      CHARACTER*256 :: cmd, template, modfile, trfile
!     INTEGER :: jsystem
      INTEGER :: iz, izatom, ia

      INTEGER :: nrho_zc
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: rho_zc, tmp, tmp2
      REAL(KIND=rspec) :: yp1, ypn
      TYPE (EZspline1_r8) :: spln

#ifdef PARTIAL_PLASMA_STATE
      integer :: cclist(ps_ccount)
      CHARACTER*128 :: saw0_partial_ps_filename="saw0_update_state.cdf"
      CHARACTER*128 :: saw1_partial_ps_filename="saw1_update_state.cdf"
#endif

      INQUIRE(file=trim(saw0_filename),exist=ex)
      cmd="rm -f -r "//trim(saw0_filename)

      call ps_copy_plasma_state(ps,saw0,ierr)
      if(ierr .ne. 0) then
      write(6,*) "Error when copying ps to saw0"
      istat=jsystem("date > transp_kill.dat")
      return
      endif

      saw0%t0 = saw_t0
      saw0%t1 = saw_t1

      saw0%eqdsk_file = trim(saw0_filename)//".geq"
      CALL ps_update_equilibrium(ierr,state=saw0)
      if(ierr .ne. 0) then
      write(6,*) "plasma state eq update error"
      cmd="date > transp_kill.dat"
      istat = jsystem(trim(cmd))
      return
      endif

!     scalar
      saw0%vsur = pa%vsur
!     saw0%taup = pa%taup
!     saw0%edge_model = 1
!     saw0%imp_edge_model = 0

!     do i=1, saw0%nspec_gas
!     saw0%gasfl(i)=pa%gasfl_ion(i)
!     saw0%recyc(i)=pa%recyc_ion(i)
!     enddo
      do i=1, saw0%nspec_gas
      saw0%sc0(i)=pa%recyc_ion(i)
      saw0%sc0(i+ps%nspec_gas)=pa%gasfl_ion(i)
      enddo

      do i=1, pa%nbeam
      saw0%power_nbi(i)=pa%power_nbi(i)
      enddo
      do i=1, pa%nicrf
      saw0%power_ic(i)=pa%power_ic(i)
      enddo


!     zone centered arrays
      nrho_zc = saw0%nrho -1
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = (saw0%rho(1:nrho_zc) + saw0%rho(2:saw0%nrho))/2
      ALLOCATE(tmp2(pa%nrho))


      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%ne,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%ns(1:nrho_zc,0) = tmp(1:nrho_zc)


      CALL EZspline_setup(spln,pa%nhy,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      j=0
      if(eq%tsc_spec(1) .gt. 0) then
      j=j+1
      saw0%ns(1:nrho_zc,j) = pa%frac_h*tmp(1:nrho_zc)        !H
      endif
      if(eq%tsc_spec(2) .gt. 0) then
      j=j+1
      saw0%ns(1:nrho_zc,j) = pa%frac_d*tmp(1:nrho_zc)        !D
      endif
      if(eq%tsc_spec(3) .gt. 0) then
      j=j+1
      saw0%ns(1:nrho_zc,j) = pa%frac_t*tmp(1:nrho_zc)        !T
      endif

      if(eq%tsc_spec(4) .gt. 0) then
      j=j+1
      CALL EZspline_setup(spln,pa%nhe,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%ns(1:nrho_zc,j) = tmp(1:nrho_zc)                   !He
      endif

      do i=1,imp%nimpchg
      j=saw0%nspec_th-imp%nimpchg+i
      tmp2(:pa%nrho)=pa%nimp(:pa%nrho,i)
      CALL EZspline_setup(spln,tmp2,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%ns(1:nrho_zc,j) = tmp(1:nrho_zc) 
      enddo

      CALL EZspline_setup(spln,pa%te,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%ts(1:nrho_zc,0) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%ti,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      DO i = 1, saw0%nspec_th
         saw0%ts(1:nrho_zc,i) = tmp(1:nrho_zc)
      ENDDO

      CALL EZspline_setup(spln,pa%zeff,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
!     saw0%zeff_in(1:nrho_zc) = tmp(1:nrho_zc)
      saw0%zeff(1:nrho_zc) = tmp(1:nrho_zc)
!     need to separate out thermal vs fast ions for zeff
      saw0%zeff_th(1:nrho_zc) = tmp(1:nrho_zc)


      CALL EZspline_setup(spln,pa%vpars,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      DO i = 1, saw0%nspec_th
         saw0%v_pars(1:nrho_zc,i) = tmp(1:nrho_zc)
      ENDDO

!     CALL EZspline_setup(spln,pa%qimp,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw0%q_impurity(1:nrho_zc) = tmp(1:nrho_zc)*ps_xe

!     CALL EZspline_setup(spln,pa%aimp,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw0%m_impurity(1:nrho_zc) = tmp(1:nrho_zc)*ps_mp

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)

      DEALLOCATE(rho_zc,tmp)

!     interval boundary oriented data
      nrho_zc = saw0%nrho
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = saw0%rho(1:nrho_zc)

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

!     CALL EZspline_setup(spln,pa%adn,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw0%diffN_s(1:nrho_zc,1) = tmp(1:nrho_zc)

!     CALL EZspline_setup(spln,pa%adp,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw0%total_entropy(1:nrho_zc) = tmp(1:nrho_zc)      

!     CALL EZspline_setup(spln,pa%ade,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw0%e_entropy(1:nrho_zc) = tmp(1:nrho_zc)      

      CALL EZspline_setup(spln,pa%adi,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%iota(1:nrho_zc) = tmp(1:nrho_zc)      

!     CALL EZspline_setup(spln,pa%chie,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw0%chi_e(1:nrho_zc) = tmp(1:nrho_zc)      

!     CALL EZspline_setup(spln,pa%chii,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw0%chi_i(1:nrho_zc) = tmp(1:nrho_zc)

!     CALL EZspline_setup(spln,pa%chio,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw0%chi_omegat(1:nrho_zc) = tmp(1:nrho_zc)

!     CALL EZspline_setup(spln,pa%d_s,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     do i= 1, saw0%nspec_th
!     saw0%d_s(1:nrho_zc,i) = tmp(1:nrho_zc)
!     enddo

      CALL EZspline_setup(spln,pa%eta_p,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%eta_parallel(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%omegat,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%omegat(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%vloop,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%v_loop(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)

      DEALLOCATE(rho_zc,tmp)

!     radiated power, unit is watts/zone
      nrho_zc = saw0%nrho_rad 
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = saw0%rho_rad(1:nrho_zc)

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%prad,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%prad(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_br,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%prad_br(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_cy,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%prad_cy(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_li,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%prad_li(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
      DEALLOCATE(rho_zc,tmp)

      nrho_zc = saw0%nrho_eq -1
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = (saw0%rho_eq(1:nrho_zc) + saw0%rho_eq(2:saw0%nrho_eq))/2

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%pres,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%p_eq(1:nrho_zc) = tmp(1:nrho_zc)
 
      CALL EZspline_setup(spln,pa%qprof,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw0%q_eq(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
      DEALLOCATE(rho_zc,tmp,tmp2)
!
#ifdef PARTIAL_PLASMA_STATE
      write(*,*) &
      "JCdebug sawtooth_mod: put_plasma_state_saw0: partail ps store"

!     call PS_RESTORE_HASH_CODES(ierr, state=saw0)

      call ps_cclist_remove('*',cclist,ierr)
      call ps_cclist_add('PLASMA',cclist,ierr)
!     call ps_cclist_add('EQ',cclist,ierr)
!     call ps_cclist_add('NBI',cclist,ierr)
!     call ps_cclist_add('IC',cclist,ierr)
!     call ps_cclist_add('LH',cclist,ierr)
!     call ps_cclist_add('EC',cclist,ierr)
!     call ps_cclist_add('RAD',cclist,ierr)
!     call ps_cclist_add('GAS',cclist,ierr)

      !  cclist contains the PLASMA component only.
      call ps_update_report(ierr,state=saw0,cclist=cclist,              &
                        report_filename='myupdate.list')
!       call ckerr('ps_update_report (7)')

!     call ps_copy_plasma_state(saw0, aux, ierr, cclist=cclist)
!     call ckerr('ps_copy_plasma_state')

      call ps_write_update_file(trim(saw0_partial_ps_filename), ierr,   &
                                state=saw0, cclist=cclist)
!     call ps_store_plasma_state(ierr, filename=trim(saw0_filename), state=saw0)
!     call ckerr('ps_store_plasma_state') 
!#else
#endif
      call ps_store_plasma_state(ierr,trim(saw0_filename),state=saw0)
      if(ierr .ne. 0) then
      write(6,*) "Error when storing saw0"
      istat=jsystem("date > transp_kill.dat")
      endif

      return
      END SUBROUTINE put_plasma_state_saw0


      SUBROUTINE put_plasma_state_saw1

      USE plasma_state_mod
      USE tsc_ps_mod, ONLY : pa, eq, imp
      USE trtsc
      USE EZspline_obj 
      USE EZspline

      IMPLICIT NONE

      INTEGER :: istat
      INTEGER :: i, j, k, l, m, ii, nout=6
      LOGICAL :: ex, l_first=.true.
      CHARACTER*256 :: cmd, template, modfile, trfile
!     INTEGER :: jsystem
      INTEGER :: iz, izatom, ia

      INTEGER :: nrho_zc
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: rho_zc, tmp, tmp2
      REAL(KIND=rspec) :: yp1, ypn
      TYPE (EZspline1_r8) :: spln

#ifdef PARTIAL_PLASMA_STATE
      integer :: cclist(ps_ccount)
      CHARACTER*128 :: saw0_partial_ps_filename="saw0_update_state.cdf"
      CHARACTER*128 :: saw1_partial_ps_filename="saw1_update_state.cdf"
#endif

      INQUIRE(file=trim(saw1_filename),exist=ex)
      cmd="rm -f -r "//trim(saw1_filename)

      call ps_copy_plasma_state(ps,saw1,ierr)
      if(ierr .ne. 0) then
      write(6,*) "Error when copying ps to saw1"
      istat=jsystem("date > transp_kill.dat")
      return
      endif

      saw1%t0 = saw_t0
      saw1%t1 = saw_t1

      saw1%eqdsk_file = trim(saw1_filename)//".geq"
      CALL ps_update_equilibrium(ierr,state=saw1)
      if(ierr .ne. 0) then
      write(6,*) "plasma state eq update error"
      cmd="date > transp_kill.dat"
      istat = jsystem(trim(cmd))
      return
      endif

!     scalar
      saw1%vsur = pa%vsur
!     saw1%taup = pa%taup
!     saw1%edge_model = 1
!     saw1%imp_edge_model = 0

!     do i=1, saw1%nspec_gas
!     saw1%gasfl(i)=pa%gasfl_ion(i)
!     saw1%recyc(i)=pa%recyc_ion(i)
!     enddo
      do i=1, saw1%nspec_gas
      saw1%sc0(i)=pa%recyc_ion(i)
      saw1%sc0(i+ps%nspec_gas)=pa%gasfl_ion(i)
      enddo

      do i=1, pa%nbeam
      saw1%power_nbi(i)=pa%power_nbi(i)
      enddo
      do i=1, pa%nicrf
      saw1%power_ic(i)=pa%power_ic(i)
      enddo


!     zone centered arrays
      nrho_zc = saw1%nrho -1
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = (saw1%rho(1:nrho_zc) + saw1%rho(2:saw1%nrho))/2
      ALLOCATE(tmp2(pa%nrho))


      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%ne,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%ns(1:nrho_zc,0) = tmp(1:nrho_zc)


      CALL EZspline_setup(spln,pa%nhy,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      j=0
      if(eq%tsc_spec(1) .gt. 0) then
      j=j+1
      saw1%ns(1:nrho_zc,j) = pa%frac_h*tmp(1:nrho_zc)        !H
      endif
      if(eq%tsc_spec(2) .gt. 0) then
      j=j+1
      saw1%ns(1:nrho_zc,j) = pa%frac_d*tmp(1:nrho_zc)        !D
      endif
      if(eq%tsc_spec(3) .gt. 0) then
      j=j+1
      saw1%ns(1:nrho_zc,j) = pa%frac_t*tmp(1:nrho_zc)        !T
      endif

      if(eq%tsc_spec(4) .gt. 0) then
      j=j+1
      CALL EZspline_setup(spln,pa%nhe,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%ns(1:nrho_zc,j) = tmp(1:nrho_zc)                   !He
      endif

      do i=1,imp%nimpchg
      j=saw1%nspec_th-imp%nimpchg+i
      tmp2(:pa%nrho)=pa%nimp(:pa%nrho,i)
      CALL EZspline_setup(spln,tmp2,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%ns(1:nrho_zc,j) = tmp(1:nrho_zc) 
      enddo

      CALL EZspline_setup(spln,pa%te,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%ts(1:nrho_zc,0) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%ti,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      DO i = 1, saw1%nspec_th
         saw1%ts(1:nrho_zc,i) = tmp(1:nrho_zc)
      ENDDO

      CALL EZspline_setup(spln,pa%zeff,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
!     saw1%zeff_in(1:nrho_zc) = tmp(1:nrho_zc)
      saw1%zeff(1:nrho_zc) = tmp(1:nrho_zc)
      saw1%zeff_th(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%vpars,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      DO i = 1, saw1%nspec_th
         saw1%v_pars(1:nrho_zc,i) = tmp(1:nrho_zc)
      ENDDO

!     CALL EZspline_setup(spln,pa%qimp,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw1%q_impurity(1:nrho_zc) = tmp(1:nrho_zc)*ps_xe

!     CALL EZspline_setup(spln,pa%aimp,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw1%m_impurity(1:nrho_zc) = tmp(1:nrho_zc)*ps_mp

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)

      DEALLOCATE(rho_zc,tmp)

!     interval boundary oriented data
      nrho_zc = saw1%nrho
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = saw1%rho(1:nrho_zc)

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

!     CALL EZspline_setup(spln,pa%adn,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw1%diffN_s(1:nrho_zc,1) = tmp(1:nrho_zc)

!     CALL EZspline_setup(spln,pa%adp,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw1%total_entropy(1:nrho_zc) = tmp(1:nrho_zc)      

!     CALL EZspline_setup(spln,pa%ade,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw1%e_entropy(1:nrho_zc) = tmp(1:nrho_zc)      

      CALL EZspline_setup(spln,pa%adi,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%iota(1:nrho_zc) = tmp(1:nrho_zc)      

!     CALL EZspline_setup(spln,pa%chie,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw1%chi_e(1:nrho_zc) = tmp(1:nrho_zc)      

!     CALL EZspline_setup(spln,pa%chii,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw1%chi_i(1:nrho_zc) = tmp(1:nrho_zc)

!     CALL EZspline_setup(spln,pa%chio,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     saw1%chi_omegat(1:nrho_zc) = tmp(1:nrho_zc)

!     CALL EZspline_setup(spln,pa%d_s,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     do i= 1, saw1%nspec_th
!     saw1%d_s(1:nrho_zc,i) = tmp(1:nrho_zc)
!     enddo

      CALL EZspline_setup(spln,pa%eta_p,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%eta_parallel(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%omegat,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%omegat(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%vloop,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%v_loop(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)

      DEALLOCATE(rho_zc,tmp)

!     radiated power, unit is watts/zone
      nrho_zc = saw1%nrho_rad 
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = saw1%rho_rad(1:nrho_zc)

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%prad,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%prad(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_br,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%prad_br(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_cy,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%prad_cy(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_li,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%prad_li(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
      DEALLOCATE(rho_zc,tmp)

      nrho_zc = saw1%nrho_eq -1
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = (saw1%rho_eq(1:nrho_zc) + saw1%rho_eq(2:saw1%nrho_eq))/2

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%pres,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%p_eq(1:nrho_zc) = tmp(1:nrho_zc)
 
      CALL EZspline_setup(spln,pa%qprof,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      saw1%q_eq(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
      DEALLOCATE(rho_zc,tmp,tmp2)
!
#ifdef PARTIAL_PLASMA_STATE
      write(*,*) &
      "JCdebug sawtooth_mod: put_plasma_state_saw1: partail ps store"

!     call PS_RESTORE_HASH_CODES(ierr, state=saw1)

      call ps_cclist_remove('*',cclist,ierr)
      call ps_cclist_add('PLASMA',cclist,ierr)
!     call ps_cclist_add('EQ',cclist,ierr)
!     call ps_cclist_add('NBI',cclist,ierr)
!     call ps_cclist_add('IC',cclist,ierr)
!     call ps_cclist_add('LH',cclist,ierr)
!     call ps_cclist_add('EC',cclist,ierr)
!     call ps_cclist_add('RAD',cclist,ierr)
!     call ps_cclist_add('GAS',cclist,ierr)

      !  cclist contains the PLASMA component only.
      call ps_update_report(ierr,state=saw1,cclist=cclist,              &
                        report_filename='myupdate.list')
!       call ckerr('ps_update_report (8)')

!     call ps_copy_plasma_state(saw1, aux, ierr, cclist=cclist)
!     call ckerr('ps_copy_plasma_state')

      call ps_write_update_file(trim(saw1_partial_ps_filename), ierr,   &
                                state=saw1, cclist=cclist)
!     call ps_store_plasma_state(ierr, filename=trim(saw1_filename), state=saw1)
!     call ckerr('ps_store_plasma_state') 
!#else
#endif
      call ps_store_plasma_state(ierr,trim(saw1_filename),state=saw1)
      if(ierr .ne. 0) then
      write(6,*) "Error when storing saw1"
      istat=jsystem("date > transp_kill.dat")
      endif

      return
      END SUBROUTINE put_plasma_state_saw1

      END MODULE sawtooth_mod
