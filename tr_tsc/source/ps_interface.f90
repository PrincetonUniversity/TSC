      MODULE trtsc

      IMPLICIT NONE
      CHARACTER*256 :: xplasma_filename, xplasma_filesave 
      CHARACTER*256 :: xplasma_filelist(1)
      CHARACTER*256 :: suffix2
!cj   June 8, 2010 add transp_number, mover shot_number from char to integer
!cj   for the convenience of importing trxpl profiles.
      CHARACTER*256 :: transp_machine, transp_year, transp_number  !cj, shot_number
      INTEGER :: shot_number
      CHARACTER*20 :: swim_init =" "
      LOGICAL :: from_xplasma = .false.
      LOGICAL :: tr_init=.false.
      LOGICAL :: from_xplasma_p = .false. , read_xplasma_p=.true.
      LOGICAL :: from_xplasma_j = .false. , read_xplasma_j=.true.
      LOGICAL :: use_transp= .false., use_plasma_state=.false.
      LOGICAL :: use_swim= .false. 
      LOGICAL :: use_tsc=.true., use_tsc_beam=.true.
      LOGICAL :: use_tsc_fp=.true., use_tsc_ich=.true.
      LOGICAL :: use_tsc_lhh=.true., use_tsc_ech=.true.
      LOGICAL :: use_user_beam=.false., use_user_lhh=.false.
      LOGICAL :: use_user_radi=.false.,use_user_ech=.false.
      LOGICAL :: use_user_fw=.false., use_user_rot=.false.
      LOGICAL :: use_user_ne=.false.,use_user_chie=.false.
      LOGICAL :: use_user_pres=.false.
      LOGICAL :: wrtrst=.false.
      INTEGER :: iseq=0, file_seq=0, istep = 0, nstep_save = 5
      INTEGER :: ipos=4950, ierr
      REAL*8 :: trtimes(10)
      REAL*8 :: tinit, tfinal
      REAL*8 :: tcpu, tcpu_past=0.0, tcpu_now, tdiff
      CHARACTER*20, DIMENSION(10) :: trevent

      END MODULE trtsc


      MODULE tsc_ps_mod
      USE kind_spec_mod
      IMPLICIT NONE

      LOGICAL :: nb_present, rf_present, lh_present, ec_present
      LOGICAL, DIMENSION(10) :: imp_present

      TYPE :: eq_ps
      INTEGER :: nrho_eq, nchi_eq, nr, nz
      INTEGER :: nspec, nspec_tsc
      INTEGER, DIMENSION(8) :: tsc_spec
      REAL(KIND=rspec) :: t0, t1
      END TYPE eq_ps

      TYPE :: pa_ps
      INTEGER :: nrho, nbeam, nicrf, nlh, necrf
      REAL(KIND=rspec) :: frac_h, frac_d, frac_t
      REAL(KIND=rspec) :: vsur, taup, q95 
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: power_nbi, power_ic, power_lh, power_ec
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: rho
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: ne, nhy, nhe
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:,:) :: nimp
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: te, ti
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: pres, qprof
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: qimp, aimp
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: gasfl_ion,recyc_ion
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: adn,adp,ade,adi
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: chie,chii,chio,d_s
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: zeff, vpars
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: prad, prad_br
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: prad_cy, prad_li
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: omegat, eta_p 
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: vloop 
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: pohm, curbs, qie
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: curoh
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: pelh
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: curlh
      END TYPE pa_ps

      TYPE :: imp_ps
      INTEGER :: nimpchg
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: chg_imp,           &
     &                                               chg_imp_atom,      &
     &                                               mass_imp
      END TYPE imp_ps

      TYPE(eq_ps) :: eq
      TYPE(pa_ps) :: pa
      TYPE(imp_ps) :: imp
      END MODULE tsc_ps_mod



      MODULE plasma_state_in_mod

      CONTAINS

      subroutine check_auxiliary_heating(nb,ic,lh,ec)

      USE tsc_ps_mod
      USE clinam, ONLY : beamp, picrh, plhamp, pecrh

      LOGICAL :: nb, ic, lh, ec

      nb=.false.
      if(any(beamp > 0.0)) nb=.true.
      ic=.false.
      if(any(picrh > 0.0)) ic=.true.
      lh=.false.
      if(any(plhamp > 0.0)) lh=.true.
      ec=.false.
      if(any(pecrh > 0.0)) ec=.true.

      return
      end subroutine check_auxiliary_heating

      subroutine check_impurities(imp_status)

      USE tsc_ps_mod
      USE clinam

      LOGICAL, DIMENSION(*) :: imp_status
      integer :: j

      imp_status(1:pimp) = .false.

      do j = 1, pimp
!
!       the following line is no longer needed since acoef(j+853) 
!       was transfered to fraciv in inpt
!       if(acoef(j+853) .gt. 1.0d-10) imp_status(j)=.true.
        do i=1,ntpts
           if(fraciv(j,i).gt. 1.0d-10) imp_status(j)=.true.
        enddo
      enddo

      return
      end subroutine check_impurities

      subroutine get_tsc_data(a, nimp, nchrgmx, imp, zmp, nch)
      use clinam
      use radtab

      integer :: nimp, nchrgmx, imp, i
      real*8 :: zmp
      real*8 :: nch(*),a(*)

      do i=1, 5000
        a(i)=acoef(i)
      enddo
      do i=1, pimp
        nch(i)=nchrgsr(i)
      enddo
      nimp=pimp
      nchrgmx = pchrgmx
      imp = iimp
      zmp = zimp
      return
      end subroutine get_tsc_data

      subroutine machine_config

      USE clinam, ONLY : acoef
      USE trtsc, ONLY: transp_machine
      USE plasma_state_mod
      USE kind_spec_mod 
      USE clinam, ONLY : ambeam, ebeamkev, ineg

      IMPLICIT NONE

      INTEGER :: istat, ierr, i, j
!cj nov-29-2010 acoef(4994>0) extract nubeam power and voltage
      type (plasma_state) :: gog
      character*80 gog_filename 

!cj dec-07-2016 change machine.mdescr reading to be more general
!cj dec-07-2016 user must provide it in the format of
!cj dec-07-2016 XX.mdescr, where XX represents D3D, NSTX, ITER, or a machine name
      character*20 machine_mdescr

      
!     NBI variables
      INTEGER, DIMENSION(20) :: nbshapa
      write(6,*) "call to machine_config", ebeamkev, ambeam
      
      if(allocated(ps%nbi_src_name)) return
      write(6,*) "ps_interface: beam configuration is now being defined"
      call flush(6)
!
      ierr = 0
      call ps_alloc_plasma_state(ierr)
      if( ierr .ne. 0 ) then
        write(6,*) "ps_interface: Error allocating nbi arrays"
        ineg=61
        return
      endif
!
!cj dec-07-2016
!cj dec-07-2016      if (trim(transp_machine) .eq. "d3d" .or.                          &
!cj dec-07-2016     &    trim(transp_machine) .eq. "D3D" ) then
!cj dec-07-2016!
!cj dec-07-2016!...........new method to read in Neutral Beam parameters
!cj dec-07-2016!cj nov 2010 get from mdescr           ps%nbeam = 7
!cj dec-07-2016            call ps_mdescr_read("D3D.mdescr",ierr,state=ps)
!cj dec-07-2016      ps%lock_machine_descr = ps_unlocked
!cj dec-07-2016            write(6,*) "ps_interface: read D3D.mdescr", ps%nbeam
!cj dec-07-2016      endif
!cj dec-07-2016!
!cj dec-07-2016      if (trim(transp_machine) .eq. "nstx" .or.                         &
!cj dec-07-2016     &    trim(transp_machine) .eq. "NSTX" ) then
!cj dec-07-2016!
!cj dec-07-2016!...........new method to read in Neutral Beam parameters
!cj dec-07-2016!cj nov 2010 get from mdescr           ps%nbeam = 3
!cj dec-07-2016            call ps_mdescr_read("NSTX.mdescr",ierr,state=ps)
!cj dec-07-2016      ps%lock_machine_descr = ps_unlocked
!cj dec-07-2016            write(6,*) "ps_interface: read NSTX.mdescr", ps%nbeam
!cj dec-07-2016      endif
!cj dec-07-2016!
!cj dec-07-2016      if (trim(transp_machine) .eq. "iter" .or.                         &
!cj dec-07-2016     &    trim(transp_machine) .eq. "ITER" ) then
!cj dec-07-2016!
!cj dec-07-2016!...........new method to read in Neutral Beam parameters
!cj dec-07-2016!cj nov 2010 get from mdescr           ps%nbeam = 2
!cj dec-07-2016            call ps_mdescr_read("ITER_SJ.mdescr",ierr,state=ps)
!cj dec-07-2016      ps%lock_machine_descr = ps_unlocked
!cj dec-07-2016            write(6,*) "ps_interface: read ITER_SJ.mdescr", ps%nbeam
!cj dec-07-2016      endif
      machine_mdescr=trim(transp_machine)//".mdescr"
            call ps_mdescr_read(trim(machine_mdescr),ierr,state=ps)
      ps%lock_machine_descr = ps_unlocked
            write(6,*) "ps_interface: read ", trim(machine_mdescr), ps%nbeam
!
      if( ierr .ne. 0) then
         write(6,*) "ps_interface: Error loading NBI configuration"
         ineg=61
         return
      else
         write(6,*) "ps_interface: Successful loading NBI configuration"
         if( acoef(4994) .gt. 0 ) then
            !cj write wall_data for plasma state
            write(6,*) "ps_interface: write nubeam wall data into file"
            open(51,file="wall_data",status="unknown",iostat=istat)
               write(51,"('&wall_data')")
               write(51,"('nwall=',I3)") ps%num_rzlim
               write(51,"('rwall=')")
               write(51,7991) (ps%rlim(i),i=1,ps%num_rzlim)
               write(51,"('zwall=')")
               write(51,7991) (ps%zlim(i),i=1,ps%num_rzlim)
               write(51,"('/')")
            close(51)
         endif
      endif
 7991 format(1p3(2x,e15.7)) 


!     shot configuration fixed

      if( acoef(4994) .gt. 0 ) then
         gog_filename="xyz/gog_ps.cdf"
         call ps_get_plasma_state(ierr,trim(gog_filename),state=gog)
            if(ierr .ne. 0) then
               write(*,*) "trxpl: err ps_get_plasma_state gog"
               ineg=62
               return
            endif
         ps%nbion(1:ps%nbeam) = gog%nbion(1:ps%nbeam) 
         write(*,*) "trxpl: ok nbion gog"
      else
      do i=1, ps%nbeam
      if( int(ambeam+0.5) .eq. 1) then
      ps%nbion(i)="H"
      endif
      if( int(ambeam+0.5) .eq. 2) then
      ps%nbion(i)="D"
      endif
      if( int(ambeam+0.5) .eq. 3) then
      ps%nbion(i)="T"
      endif
      enddo
      endif


!     simulation_init fixed

      ps%nspec_beam = 1
      ierr = 0
      call ps_alloc_plasma_state(ierr)
      if( ierr .ne. 0 ) then
      write(6,*) "ps_interface: Error allocating nbi arrays"
      return
      endif

 
      do i = 1, ps%nspec_beam
        CALL ps_species_convert( 1,  1, 2,                              &
     &       ps%qatom_snbi(i), ps%q_snbi(i), ps%m_snbi(i), ierr)
      enddo
      CALL ps_label_species(ierr)

!     state variables fixed

      if( acoef(4994) .gt. 0 ) then
         ps%kvolt_nbi(1:ps%nbeam) = gog%kvolt_nbi(1:ps%nbeam) 
         write(*,*) "trxpl: ok kvolt_nbi gog"
      else
      do i = 1, ps%nbeam
      ps%kvolt_nbi(i) = ebeamkev
      enddo
      endif

!     done
      return
      end subroutine machine_config


      subroutine init_plasma_state
      USE plasma_state_mod
      USE tsc_ps_mod
      USE trtsc
      USE kind_spec_mod 
!     USE clinam, ONLY : acoef, iimp, zimp
!     USE radtab, ONLY : nchrgsr, pimp, pchrgmx
      USE clinam, ONLY : ambeam, ebeamkev, ineg, nbeami, q_snbi, nspec_beam, npsit, ppsi

      IMPLICIT NONE
      INTEGER :: istat 
      INTEGER :: i, j, k, l, m, ii, nout 
      LOGICAL :: ex, l_first=.true.
      CHARACTER*256 :: cmd, template, modfile, trfile
      INTEGER :: iz, izatom, ia
      INTEGER :: nspec, nchrgs
      INTEGER :: pimp, pchrgmx, iimp
      REAL*8 :: zimp 
      REAL*8 :: acoef(6000), nchrgsr(10)

      logical :: smooth = .true. , norm=.false.  
      INTEGER nrho,new_nrho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: new_rho, tmp_nbeami

#ifdef PARTIAL_PLASMA_STATE
      integer :: cclist(ps_ccount)
      CHARACTER*128 ::                                                  &
          xplasma_partial_ps_filename="ps_update_state_init.cdf"
#endif

      xplasma_filename = xplasma_filesave
      write(6,*) "ps%id_pbe", ps%id_pbe, ex, trim(xplasma_filename)
      INQUIRE(file=trim(xplasma_filename),exist=ex)
       if(ex) then
        CALL ps_get_plasma_state(ierr, filename=trim(xplasma_filename),state=ps)
        if(ierr .ne. 0) then
          write(6,*) "ps_interface: Plasma State initial read error"
          call flush(6)
          ineg=61
        else
          from_xplasma=.true.
          from_xplasma_p=.true.
          from_xplasma_j=.true.
          tr_init=.true.

          ! not using trxpl profile, get nbeami from nubeam
          if(allocated(ps%nbeami)) then 
             write(*,*) "ps_interface: read nbeami from nubeam ps"
             nrho=ps%nrho
             nspec_beam=ps%nspec_beam
             ALLOCATE(nbeami(ppsi,nspec_beam),q_snbi(nspec_beam)) 
             do i=1,nspec_beam
               q_snbi(i)=ps%q_snbi(i)/(1.6022e-19)
             enddo

             new_nrho = npsit
             ALLOCATE( new_rho(new_nrho), tmp_nbeami(new_nrho-1))
               if (istat .ne. 0) stop 'Allocation Error : new_rho' 
             DO i = 1, new_nrho
                !cj square root of toroidal flux, compatible with trxpl output
                new_rho(i) = sqrt( (real(i)-1.)/(real(new_nrho)-1.) )
             ENDDO
             new_rho(1) = 0.
             new_rho(new_nrho) = 1.

             DO i = 1,nspec_beam
                ii = ps%id_nbeami(i)
                call ps_rho_rezone(new_rho,ii,tmp_nbeami,ierr,             &
                        state=ps,zonesmoo=smooth, nonorm=norm)
                IF(ierr.ne.0) write(*,*) "nbeami ps_rho_rezone1 error"
                nbeami(2:new_nrho,i) = tmp_nbeami(1:new_nrho-1)
                write(*,*) "i, q_snbi(i)", i, q_snbi(i)
                write(*,*) "nbeami"
                write(*,1001) (nbeami(j,i),j=1,npsit)
       1001     format(1p10e12.4)
             ENDDO   ! i
             DEALLOCATE(new_rho, tmp_nbeami)
          endif !(allocated(ps%nbeami)) then 

        endif
       endif
      write(6,*) "ps%id_pbe", ps%id_pbe, ex, trim(xplasma_filename)

!#ifdef PARTIAL_PLASMA_STATE
!       call PS_SAVE_HASH_CODES(ierr, state=ps)
!       write(*,*) "JCdebug ps_interface: init_plasma_state partial ps"
!        call ps_read_update_file(trim(xplasma_filename),ierr,state=ps,  &
!     &                           update_complete=.true.) 
!#else
!#endif

      ps%runid = trim(suffix2)
      ps%tokamak_id = trim(transp_machine)
      ps%global_label="plasma state"//" "//trim(suffix2)
!cj   June 8, 2010 add transp_number, mover shot_number from char to integer
!cj   for the convenience of importing trxpl profiles.
!     write(51,"(a)") trim(shot_number)
!     rewind(51)
!     read(51,"(i10)") ps%shot_number
!     close(51,status="delete")
      ps%shot_number = shot_number
      ps%tinit = tinit
      ps%tfinal = tfinal
      ps%t0=tinit
      ps%t1=ps%t0

!     ps_debug = 2

      eq%nspec_tsc = 8
      eq%tsc_spec = 0

      call get_tsc_data(acoef,pimp,pchrgmx,iimp,zimp,nchrgsr)
      if(acoef(4955) .ge. 1.d-4) eq%tsc_spec(1)=1                       !H
      if(acoef(4955) .lt. 0.999999) then
      if(acoef(113) .gt. 0.0001) eq%tsc_spec(2) = 1                     !D
      if(acoef(113) .lt. 0.999999) eq%tsc_spec(3)=1                     !T
      endif
      if(eq%tsc_spec(3) .eq. 1) eq%tsc_spec(4)=1                        !He
      nspec = 0
      do i=1, eq%nspec_tsc
      nspec=nspec+eq%tsc_spec(i)
      enddo

      if(.not.allocated(imp%chg_imp))                                   &
     &                       allocate(imp%chg_imp(pimp*pchrgmx))
      if(.not.allocated(imp%chg_imp_atom))                              &
     &                       allocate(imp%chg_imp_atom(pimp*pchrgmx))
      if(.not.allocated(imp%mass_imp))                                  &
     &                       allocate(imp%mass_imp(pimp*pchrgmx))
      if (iimp .eq. 0) then
      imp%nimpchg = 1
      imp%chg_imp(1) = zimp
      imp%chg_imp_atom(1) = zimp
      imp%mass_imp(1) = 2.*zimp
      else

      call check_impurities(imp_present)
      i = 0
      do j=1, pimp
      nchrgs = nchrgsr(j)
      if(.not. imp_present(j)) cycle
      do k = 2, nchrgs
      i=i+1
      imp%chg_imp(i) = k-1
      imp%chg_imp_atom(i) = nchrgs-1
      imp%mass_imp(i) = 2.*(nchrgs-1)
      enddo
      enddo
      imp%nimpchg = i
      endif
      eq%nspec = nspec + imp%nimpchg

      call check_auxiliary_heating(nb_present,rf_present,lh_present,ec_present)

      if( (ps%nspec_th .eq. ps_uninit) .or. (.not. allocated(ps%q_s)))     &
     & then

      ps%nspec_th = eq%nspec

!       CALL ps_label_plasma_state(TRIM(suffix2), ierr)

        CALL ps_alloc_plasma_state(ierr)

        CALL ps_species_convert(-1, -1, 0,                              &
     &       ps%qatom_s(0), ps%q_s(0), ps%m_s(0), ierr)

        j=0
        if( eq%tsc_spec(1) .gt. 0) then
        j=j+1
        CALL ps_species_convert( 1,  1, 1,                              &
     &       ps%qatom_s(j), ps%q_s(j), ps%m_s(j), ierr)
        endif

        if( eq%tsc_spec(2) .gt. 0) then
        j=j+1
        CALL ps_species_convert( 1,  1, 2,                              &
     &       ps%qatom_s(j), ps%q_s(j), ps%m_s(j), ierr)
        endif

        if( eq%tsc_spec(3) .gt. 0) then
        j=j+1
        CALL ps_species_convert( 1,  1, 3,                              &
     &       ps%qatom_s(j), ps%q_s(j), ps%m_s(j), ierr)
        endif

        if( eq%tsc_spec(4) .gt. 0) then
        j=j+1
        CALL ps_species_convert( 2,  2, 4,                              &
     &       ps%qatom_s(j), ps%q_s(j), ps%m_s(j), ierr)

        endif


        if(imp%nimpchg .ne. 0) then
        j = 0
        do i = ps%nspec_th-imp%nimpchg+1, ps%nspec_th
        j = j + 1
        iz = imp%chg_imp(j) + 0.01
        izatom = imp%chg_imp_atom(j) +0.01
        ia = imp%mass_imp(j) + 0.01
        call ps_species_convert(izatom, iz, ia,                         &
     &                          ps%qatom_s(i),                          &
     &                          ps%q_s(i), ps%m_s(i), ierr )

        write(6,*) "impurity, iz, iza, ia ", iz, izatom, ia
        enddo
        endif
        write(6,*) "no. of thermal species : ", ps%nspec_th
        call flush(6)

        CALL ps_label_species(ierr)


!        ierr=0
!#ifdef PARTIAL_PLASMA_STATE
!        write(*,*) &
!        "JCdebug ps_interface: init_plasma_state: partial ps store I"
!
!       call PS_RESTORE_HASH_CODES(ierr, state=ps)
!
!        call ps_cclist_remove('*',cclist,ierr)
!        call ps_cclist_add('PLASMA',cclist,ierr)
!        call ps_cclist_add('EQ',cclist,ierr)
!        call ps_cclist_add('NBI',cclist,ierr)
!        call ps_cclist_add('IC',cclist,ierr)
!        call ps_cclist_add('LH',cclist,ierr)
!        call ps_cclist_add('EC',cclist,ierr)
!        call ps_cclist_add('RAD',cclist,ierr)
!        call ps_cclist_add('GAS',cclist,ierr)
!        call ps_cclist_add('LMHD',cclist,ierr)
!        call ps_cclist_add('RIPPLE',cclist,ierr)
!        call ps_cclist_add('ANOM',cclist,ierr)
!        call ps_cclist_add('FUS',cclist,ierr)
!        call ps_cclist_add('RUNAWAY',cclist,ierr)
!
!        !  cclist contains the PLASMA component only.  
!        call ps_update_report(ierr,state=aux,cclist=cclist,              &
!                        report_filename='myupdate_1.list')

!       call ps_copy_plasma_state(ps, aux, ierr, cclist=cclist) 
!       call ps_write_update_file(trim(xplasma_filename),ierr,state=ps, &
!        call ps_write_update_file(trim(xplasma_partial_ps_filename),    &
!                                  ierr,state=ps, cclist=cclist)
!       call ps_store_plasma_state(ierr, filename=xplasma_filename, state=ps)
!#else
!#endif
!        CALL ps_store_plasma_state(ierr,filename=xplasma_filename)
!        if(ierr.ne.0) then
!        write(6,*) "ps_interface: store plasma state error 1"
!        call flush(6)
!        return
!        endif

        if(use_transp) then
        template = "ptransp_template.dat"
        modfile  = "ptransp_modify.dat"
        trfile = trim(suffix2)//"TR.DAT"
        cmd= "namelist_adapt "//trim(template)//" "//                   &
     &        trim(modfile)//" "//trim(xplasma_filename)//" "           &
     &        //trim(trfile)
        ierr = jsystem(trim(cmd))
        if( ierr .ne.0) then
        write(6,*) "ps_interface: TRANSP namelist adapt error"
        call flush(6)
        cmd="date > transp_kill.dat"
        istat = jsystem(trim(cmd))
        return
        endif
        endif

!#ifdef PARTIAL_PLASMA_STATE
!       call PS_SAVE_HASH_CODES(ierr, state=ps)
!        write(*,*) &
!        "JCdebug ps_interface: init_plasma_state: partial ps read II"
!        call ps_read_update_file(trim(xplasma_partial_ps_filename),     &
!                                 ierr,state=ps, update_complete=.true.) 
!#else
!#endif
!        CALL ps_get_plasma_state(ierr, filename=xplasma_filename)
!        if(ierr .ne. 0) then
!        write(6,*) "ps_interface: PS file read error after adapt"
!        call flush(6)
!        return
!        endif

        if(.not.use_transp) then

        if( eq%tsc_spec(4) .gt. 0) then
        ps%nspec_fusion = 1
        call ps_alloc_plasma_state(ierr)
        CALL ps_species_convert(2, 2, 4,                                &
     &       ps%qatom_sfus(1), ps%q_sfus(1), ps%m_sfus(1), ierr)
        CALL ps_label_species(ierr)
        endif

!       call ps_merge_species_lists(ierr)
        ps%ngsc0 = 0
        CALL ps_neutral_species(ierr)
        endif

!
!     new neutral gas definition
      if( ps%ngsc0 .eq. 0) then
        if( ps%nspec_gas .gt. 0 ) then
        ps%ngsc0 = 2*ps%nspec_gas
        CALL ps_alloc_plasma_state(ierr)
        do j=1, ps%nspec_gas
        cmd = ps%sgas_name(j)
        i = len(trim(cmd))
        if( cmd(i-1:i) .eq. "_0" ) then
         ii = i-2
        else if (cmd(i:i) .eq. "0" ) then
         ii = i-1
        else
         ii = i
        endif
        ps%gs_name(j)=cmd(1:i)//"rcy"
        ps%gas_atom(j) = cmd(1:ii)
        ps%gs_name(j+ps%nspec_gas)=cmd(1:i)//"gf"
        ps%gas_atom(j+ps%nspec_gas)=cmd(1:ii)
        enddo

        call ps_gsc0_species_map(ierr)
        if(ierr .ne. 0) then
        write(6,*) "ps load neutral data, error in ps_gsc0_species_map"
        return
        endif
        endif ! cj oct 13 2010 moved from L551, if( ps%nspec_gas .gt. 0 ) then

      else if(ps%ngsc0 .ne. 2*ps%nspec_gas) then
        write(6,*) "nsgs0 .ne. 2%nspec_gas"
        return
      endif

!     initialize machine configuration if swim
      if( use_swim ) then
      call machine_config
!     temporary for average gas energy
      do i=1, ps%nspec_gas
      ps%e0_av(i) = 0.015_dp
      ps%e0_av(i+ps%nspec_gas)=0.005_dp
      enddo
      if (ps%dn0out .le. 0.0) ps%dn0out = 1.0e17
      endif

!
        ierr=0
#ifdef PARTIAL_PLASMA_STATE
!        call PS_RESTORE_HASH_CODES(ierr, state=ps)
!
        call ps_cclist_remove('*',cclist,ierr)
        call ps_cclist_add('*',cclist,ierr)
!       call ps_cclist_add('PLASMA',cclist,ierr)
!       call ps_cclist_add('EQ',cclist,ierr)
!       call ps_cclist_add('NBI',cclist,ierr)
!       call ps_cclist_add('IC',cclist,ierr)
!       call ps_cclist_add('LH',cclist,ierr)
!       call ps_cclist_add('EC',cclist,ierr)
!       call ps_cclist_add('RUNAWAY',cclist,ierr)
!       call ps_cclist_add('FUS',cclist,ierr)
!       call ps_cclist_add('RAD',cclist,ierr)
!       call ps_cclist_add('GAS',cclist,ierr)
!       call ps_cclist_add('LMHD',cclist,ierr)
!       call ps_cclist_add('RIPPLE',cclist,ierr)
!       call ps_cclist_add('ANOM',cclist,ierr)
 
        !  cclist contains the PLASMA component only.
        call ps_update_report(ierr,state=ps,cclist=cclist,             &
                        report_filename='ps_update_state_init.list')

!       call ps_copy_plasma_state(ps, aux, ierr, cclist=cclist) 
!       call ps_write_update_file(trim(xplasma_filename),ierr,state=ps, &
        call ps_write_update_file(trim(xplasma_partial_ps_filename),    &
                                  ierr, state=ps, cclist=cclist)
!       call ps_store_plasma_state(ierr, filename=xplasma_filename, state=ps)
!#else
#endif
        CALL ps_store_plasma_state(ierr,filename=xplasma_filename)
        if(ierr.ne.0) then
        write(6,*) "ps_interface: store plasma state error 2"
        call flush(6)
        return
        endif


      write(6,*) "initial plasma state created successfully"
      call flush(6)
      endif
!cj end of if( (ps%nspec_th .eq. ps_uninit) .or. (.not. allocated(ps%q_s)))

      ps%tfinal = tfinal

      RETURN

      end subroutine init_plasma_state
      END MODULE plasma_state_in_mod

      subroutine set_grid_size(nrho)
      use clinam, ONLY : acoef

      integer :: nrho

      nrho = int(acoef(4953))
      if(nrho .le. 51) nrho=51
      return
      end subroutine set_grid_size


      MODULE plasma_state_eq_mod
      IMPLICIT NONE

      CONTAINS

      subroutine put_plasma_state_eq
      USE plasma_state_mod
      USE tsc_ps_mod
      USE trtsc
      USE kind_spec_mod 
      IMPLICIT NONE


      INTEGER :: istat 
      INTEGER :: i, j, k, l, m, ii, nout, nrho
      LOGICAL :: ex, l_first=.true.
      CHARACTER*256 :: cmd, template, modfile, trfile
      INTEGER :: iz, izatom, ia

#ifdef PARTIAL_PLASMA_STATE
      integer :: cclist(ps_ccount)
      CHARACTER*128 ::                                                  &
          xplasma_partial_ps_filename="ps_update_state_eq.cdf"
#endif

      xplasma_filename = xplasma_filesave
      INQUIRE(file=trim(xplasma_filename),exist=ex)


      IF( .NOT. ex ) THEN
      write(6,*) "failed to find plasma state file"
      call flush(6)
      return
      ENDIF

#ifdef PARTIAL_PLASMA_STATE
      call PS_SAVE_HASH_CODES(ierr, state=ps)
      write(*,*) "JCdebug ps_interface: put_plasma_state_eq partial ps"
#endif

      IF( l_first ) then
        IF(.not. allocated(ps%rho)) then

        write(6,*) "allocating plasma and eq rho grids"
        call set_grid_size(nrho)
        write(6,*) "number of grids : ", nrho
!       ps%nrho = 51
!       ps%nrho_eq = 51
        ps%nrho = nrho
        ps%nrho_eq = nrho
        ps%nrho_rad = ps%nrho
        ps%nrho_lhrf = ps%nrho
        ps%nth_eq = 101
        ps%nR = 0; ps%nZ = 0
        ps%nrho_eq_geo = ps%nrho_eq

        CALL ps_alloc_plasma_state(ierr)

        DO i = 1, ps%nrho
           ps%rho(i) = (FLOAT(i)-1.0_dp)/(FLOAT(ps%nrho)-1.0_dp)
        ENDDO
        DO i = 1, ps%nrho_eq
           ps%rho_eq(i) = (FLOAT(i)-1.0_dp)/(FLOAT(ps%nrho_eq)-1.0_dp)
        ENDDO
        DO i = 1, ps%nrho_rad
           ps%rho_rad(i) = (FLOAT(i)-1.0_dp)/(FLOAT(ps%nrho_rad)-1.0_dp)
        ENDDO
        DO i = 1, ps%nrho_lhrf
           ps%rho_lhrf(i) = (FLOAT(i)-1.0_dp)/(FLOAT(ps%nrho_lhrf)-1.0_dp)
        ENDDO
        DO i = 1, ps%nth_eq
           ps%th_eq(i) = (FLOAT(i)-1.0_dp)/(FLOAT(ps%nth_eq)-1.0_dp)
        ENDDO
           ps%th_eq = 2.0_dp*acos(-1.0_dp)*ps%th_eq
        ENDIF
        ps%rho_eq_geo = ps%rho_eq

        if(.not. use_transp) then
        if(.not.allocated(ps%rho_gas)) then
        ps%nrho_gas = 51
        CALL ps_alloc_plasma_state(ierr)
        DO i = 1, ps%nrho_gas
           ps%rho_gas(i) = (FLOAT(i)-1.0_dp)/(FLOAT(ps%nrho_gas)-1.0_dp)
        ENDDO
        endif
        endif
!
!.....DEBUG
      write(6,*) "before check fast ion species"

!      check fast ion species and allocation one more time
!      if(nb_present) then
!      if(.not.allocated(ps%power_nbi)) then
!      write(6,*) "nbi used but nbeam undefined"
!      ps%nbeam = 2
!      write(6,*) "nbeam set to : ",ps%nbeam
!      call flush(6)
!      CALL ps_alloc_plasma_state(ierr)
!      ps%nbi_src_name(1)="NB1"
!      ps%nbi_src_name(2)="NB2"
!      ps%nbion(1)="D"
!      ps%nbion(2)="D"
!      endif
!      endif

       if(rf_present) then
       if(.not.allocated(ps%power_ic )) then
       write(6,*) "icrf used but nicrf_src undefined"
       ps%nicrf_src = 1
       write(6,*) "nicrf_src set to : ",ps%nicrf_src
       call flush(6)
       CALL ps_alloc_plasma_state(ierr)
       ps%icrf_src_name(1)="IC1"
       endif
       if(ps%nspec_rfmin .eq. ps_uninit .or.                                       &
     &    ps%nspec_rfmin .eq. 0) then
        write(6,*) "icrf minority species undefined"
        ps%nspec_rfmin = 1
        write(6,*) "number of minority species set to : ",ps%nspec_rfmin
        call flush(6)
       CALL ps_alloc_plasma_state(ierr)
        write(6,*) "minority species set to He-3"
       CALL ps_species_convert(2, 2, 3,                                           &
     &   ps%qatom_rfmin(1), ps%q_rfmin(1), ps%m_rfmin(1), ierr)
       CALL ps_label_species(ierr)
!      CALL ps_merge_species_lists(ierr)
        if(ierr .gt. 0) then
        write(6,*) "convert error for rfmin"
        call flush(6)
        ps%rfmin_name="He3"
        endif
        ps%fracmin = 0.001
       endif
       endif
!
!>>>debug
      write(6,*) "before check if LH is present=", lh_present

       if(lh_present) then
       if(.not.allocated(ps%power_lh )) then
       write(6,*) "put_plasma_state_eq ps%nlhrf_src=", ps%nlhrf_src
       ps%nlhrf_src = 1 
       write(6,*) "nlhrf_src set to : ",ps%nlhrf_src
       call flush(6)
       CALL ps_alloc_plasma_state(ierr)
       ps%lhrf_src_name(1)="LH1"
       endif
       endif

       if(ec_present) then
       if(.not.allocated(ps%power_ec )) then
       write(6,*) "ecrf used but necrf_src undefined"
       ps%necrf_src = 1
       write(6,*) "necrf_src set to : ",ps%necrf_src
       call flush(6)
       CALL ps_alloc_plasma_state(ierr)
       ps%ecrf_src_name(1)="EC"
       endif
       endif

!     final merge of species, if has not been done
      write(6,*) "before check ps_merge_species_lists"
        if(.not.use_transp) then
        call ps_merge_species_lists(ierr)
        endif
!
      if(.not. use_swim) then
       cmd="cp "//trim(xplasma_filename)//" "//trim(suffix2)//"_psp.cdf"
       istat = jsystem(trim(cmd))
      endif
 
      ELSE !cj the opposite of IF( l_first ) 
      write(6,*) "before call to ps_get_plasma_state"

!#ifdef PARTIAL_PLASMA_STATE
!!     call PS_SAVE_HASH_CODES(ierr, state=ps)
!!     write(*,*) &
!!     "JCdebug ps_interface: put_p_s_eq: partial ps read"
!        call ps_read_update_file(trim(xplasma_filename),ierr,state=ps,  &
!     &                           update_complete=.true.) 
!#else
      CALL ps_get_plasma_state(ierr, filename=trim(xplasma_filename),state=ps)
!#endif
      if(ierr .ne. 0) then
      write(6,*) "ps_interface: plasma_state read error"
      cmd="date > transp_kill.dat"
      istat = jsystem(trim(cmd))
      return
      endif


      if(.not. use_swim) then
      cmd=trim(suffix2)//"_psp.cdf"
      CALL ps_commit_plasma_state(ierr, trim(cmd))
      if(ierr .ne. 0) then
      write(6,*) "ps_interface: plasma state commit error"
      cmd="date > transp_kill.dat"
      istat = jsystem(trim(cmd))
      return
      endif
      endif

      ps%tfinal = tfinal
      ENDIF
!cj end of      IF( l_first )

      if(use_transp) then
      ps%t0 = eq%t0
      ps%t1 = eq%t1
      endif


      write(6,*) "before check ps_update_equilibrium"
      ps%eqdsk_file = trim(suffix2)//"_ps.geq"
      CALL ps_update_equilibrium(ierr)
      if(ierr .ne. 0) then
      write(6,*) "ps_interface: plasma state eq update error"
      cmd="date > transp_kill.dat"
      istat = jsystem(trim(cmd))
      stop
      endif

      write(6,*) "before check ps_mhdeq_derive"
      call ps_mhdeq_derive('Everything',ierr)
#ifdef PARTIAL_PLASMA_STATE
      call PS_RESTORE_HASH_CODES(ierr, state=ps)

      call ps_cclist_remove('*',cclist,ierr)
      call ps_cclist_add('*',cclist,ierr)
!     call ps_cclist_add('PLASMA',cclist,ierr)
!     call ps_cclist_add('EQ',cclist,ierr)
!     call ps_cclist_add('NBI',cclist,ierr)
!     call ps_cclist_add('IC',cclist,ierr)
!     call ps_cclist_add('LH',cclist,ierr)
!     call ps_cclist_add('EC',cclist,ierr)
!     call ps_cclist_add('RUNAWAY',cclist,ierr)
!     call ps_cclist_add('FUS',cclist,ierr)
!     call ps_cclist_add('RAD',cclist,ierr)
!     call ps_cclist_add('GAS',cclist,ierr)
!     call ps_cclist_add('LMHD',cclist,ierr)
!     call ps_cclist_add('RIPPLE',cclist,ierr)
!     call ps_cclist_add('ANOM',cclist,ierr)

      !  cclist contains the PLASMA component only.
      call ps_update_report(ierr,state=ps,cclist=cclist,               &
                        report_filename='ps_update_state_eq.list')

!     call ps_copy_plasma_state(ps, aux, ierr, cclist=cclist) 
!     call ps_write_update_file(trim(xplasma_filename), ierr, state=ps, &
      call ps_write_update_file(trim(xplasma_partial_ps_filename),      &
                                ierr, state=ps, cclist=cclist)
!     call ps_store_plasma_state(ierr, filename=xplasma_filename, state=ps)
!#else
#endif
      CALL ps_store_plasma_state(ierr,xplasma_filename)
      l_first = .false.
      if(ierr .ne. 0) then
      write(6,*) "ps_interface: plasma state store error"
      cmd="date > transp_kill.dat"
      istat = jsystem(trim(cmd))
      return
      endif
 
      end subroutine put_plasma_state_eq

      END MODULE plasma_state_eq_mod


      MODULE plasma_state_pa_mod
      IMPLICIT NONE


      CONTAINS

      subroutine put_plasma_state_pa
      USE clinam, ONLY : acoef, ineg, npsit
      USE plasma_state_mod
      USE tsc_ps_mod , ONLY :  pa, imp, eq
      USE EZspline_obj
      USE EZspline
      USE trtsc
      USE SPDMOD, ONLY : csuma, ncurrentg

      IMPLICIT NONE

      INTEGER :: istat !,jsystem
      INTEGER :: i, j, k, l, m, nout=6, nrho_zc
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: rho_zc, tmp, tmp2,tmp_ni
      REAL(KIND=rspec) :: yp1, ypn, curtot, powtot
!     INTEGER :: jsystem
      CHARACTER*256 :: cmd
      TYPE (EZspline1_r8) :: spln

#ifdef PARTIAL_PLASMA_STATE
      integer :: cclist(ps_ccount)
      CHARACTER*128 ::                                                  &
          xplasma_partial_ps_filename="ps_update_state_pa.cdf"
#endif

!cj nov-29-2010 acoef(4994>0) extract nubeam power and voltage
      type (plasma_state) :: gog
      character*80 gog_filename

      REAL*8, ALLOCATABLE, DIMENSION(:) :: x, xh
      logical :: norm=.true.

!#ifdef PARTIAL_PLASMA_STATE
!       call PS_SAVE_HASH_CODES(ierr, state=ps)
!       write(*,*) "JCdebug ps_interface: put_plasma_state_pa partial ps"
!#endif

!     scalar
      if(use_transp) then
      ps%t0 = eq%t0
      ps%t1 = eq%t1
      endif

      if(use_swim) then
      if(trim(swim_init) .eq. "init") then
      ps%t0 = eq%t0
      ps%t1 = ps%t0
      endif
!     ps%t0 = eq%t0
!     ps%t1 = eq%t1
      endif

      ps%vsur = pa%vsur
!     ps%q95  = ps%q95
!     ps%taup = pa%taup
!     ps%edge_model = 1
!     ps%imp_edge_model = 0

!     do i=1, ps%nspec_th
!     do i=1, ps%nspec_gas
!     ps%gasfl(i)=pa%gasfl_ion(i)
!     ps%recyc(i)=pa%recyc_ion(i)
!     enddo
      do i=1, ps%nspec_gas
      ps%sc0(i)=pa%recyc_ion(i)
      ps%sc0(i+ps%nspec_gas)=pa%gasfl_ion(i)
      enddo

      if(allocated(ps%power_nbi)) then
         if( acoef(4994) .gt. 0 ) then
            gog_filename="xyz/gog_ps.cdf"
            call ps_get_plasma_state(ierr,trim(gog_filename),state=gog)
               if(ierr .ne. 0) then
                  write(*,*) "trxpl: err ps_get_plasma_state gog"
                  ineg=62
                  return
               endif
            ps%power_nbi(1:ps%nbeam) = gog%power_nbi(1:ps%nbeam) 

            ! read ne and zeff from plasma state trxpl gog
            ! nrho_zc_gog = gog%nrho -1
            ! pa%ne(1:nrho_zc_gog) = gog%ns(1:nrho_zc_gog,0)
            ! pa%zeff(1:nrho_zc_gog) = gog%ne(1:nrho_zc_gog)

         else
      do i=1, pa%nbeam
      ps%power_nbi(i)=pa%power_nbi(i)
      write(*,1810) i,pa%nbeam,ps%power_nbi(i)
 1810 format("i,pa%nbeam,ps%power_nbi(i)",2i5,1pe12.4)
      enddo
         endif
      endif
      if (allocated(ps%power_ic)) then
      do i=1, pa%nicrf
      ps%power_ic(i)=pa%power_ic(i)
      enddo
      endif
      if (allocated(ps%power_lh)) then
      do i=1, pa%nlh
      ps%power_lh(i)=pa%power_lh(i)
      enddo
      endif
      if (allocated(ps%power_ec)) then
      do i=1, pa%necrf
      ps%power_ec(i)=pa%power_ec(i)
      enddo
      endif


!     zone centered arrays
      nrho_zc = ps%nrho -1
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = (ps%rho(1:nrho_zc) + ps%rho(2:ps%nrho))/2
      ALLOCATE(tmp2(pa%nrho))
      ALLOCATE(tmp_ni(pa%nrho))

      tmp_ni = 0.D0

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

         if( acoef(4994) .le. 0 ) then
      CALL EZspline_setup(spln,pa%ne,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
         else
            write(*,*) "zeff read from trxpl", gog%nrho
            tmp(1:nrho_zc) = gog%ns(1:nrho_zc,0)
         endif
      ps%ns(1:nrho_zc,0) = tmp(1:nrho_zc)


      CALL EZspline_setup(spln,pa%nhy,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      j=0
      if(eq%tsc_spec(1) .gt. 0) then
        j=j+1
        ps%ns(1:nrho_zc,j) = pa%frac_h*tmp(1:nrho_zc)        !H
        tmp_ni = tmp_ni+pa%frac_h*tmp(1:nrho_zc)
      endif
      
      if(eq%tsc_spec(2) .gt. 0) then
        j=j+1
        ps%ns(1:nrho_zc,j) = pa%frac_d*tmp(1:nrho_zc)        !D
        tmp_ni = tmp_ni+pa%frac_d*tmp(1:nrho_zc)
      endif
      if(eq%tsc_spec(3) .gt. 0) then
        j=j+1
        ps%ns(1:nrho_zc,j) = pa%frac_t*tmp(1:nrho_zc)        !T
        tmp_ni = tmp_ni+pa%frac_t*tmp(1:nrho_zc)
      endif

      if(eq%tsc_spec(4) .gt. 0) then
        j=j+1
        CALL EZspline_setup(spln,pa%nhe,ierr)
        CALL EZspline_error(ierr)
        CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
        CALL EZspline_error(ierr)
        ps%ns(1:nrho_zc,j) = tmp(1:nrho_zc)                   !He
        tmp_ni = tmp_ni+tmp(1:nrho_zc)
      endif

      do i=1,imp%nimpchg
        j=ps%nspec_th-imp%nimpchg+i
        tmp2(:pa%nrho)=pa%nimp(:pa%nrho,i)
        CALL EZspline_setup(spln,tmp2,ierr)
        CALL EZspline_error(ierr)
        CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
        CALL EZspline_error(ierr)
        ps%ns(1:nrho_zc,j) = tmp(1:nrho_zc) 
        tmp_ni = tmp_ni+tmp(1:nrho_zc)
      enddo
      ps%ni(1:nrho_zc) = tmp_ni 

      CALL EZspline_setup(spln,pa%te,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%ts(1:nrho_zc,0) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%ti,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      DO i = 1, ps%nspec_th
         ps%ts(1:nrho_zc,i) = tmp(1:nrho_zc)
      ENDDO
      ps%ti(1:nrho_zc) = ps%ts(1:nrho_zc,1)

         if( acoef(4994) .le. 0 ) then
      CALL EZspline_setup(spln,pa%zeff,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
         else
            write(*,*) "zeff read from trxpl", gog%nrho
            tmp(1:nrho_zc) = gog%zeff(1:nrho_zc)
         endif
      ps%zeff(1:nrho_zc) = tmp(1:nrho_zc)
!
!     need to separate out thermal vs fast ions for zeff
      ps%zeff_th(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%vpars,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      DO i = 1, ps%nspec_th
         ps%v_pars(1:nrho_zc,i) = tmp(1:nrho_zc)
      ENDDO

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)

      DEALLOCATE(rho_zc,tmp)

!     surface centered data
      nrho_zc = ps%nrho
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = ps%rho(1:nrho_zc)

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%adi,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%iota(1:nrho_zc) = tmp(1:nrho_zc)      

!     CALL EZspline_setup(spln,pa%chie,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     ps%chi_e(1:nrho_zc) = tmp(1:nrho_zc)      

!     CALL EZspline_setup(spln,pa%chii,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     ps%chi_i(1:nrho_zc) = tmp(1:nrho_zc)

!     CALL EZspline_setup(spln,pa%chio,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     ps%chi_omegat(1:nrho_zc) = tmp(1:nrho_zc)

!     CALL EZspline_setup(spln,pa%d_s,ierr)
!     CALL EZspline_error(ierr)
!     CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
!     CALL EZspline_error(ierr)
!     do i= 1, ps%nspec_th
!     ps%d_s(1:nrho_zc,i) = tmp(1:nrho_zc)
!     enddo

      CALL EZspline_setup(spln,pa%eta_p,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%eta_parallel(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%omegat,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%omegat(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_setup(spln,pa%vloop,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%v_loop(1:nrho_zc) = tmp(1:nrho_zc)


      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)

      DEALLOCATE(rho_zc,tmp)

!     radiated power, unit is watts/zone
      nrho_zc = ps%nrho_rad 
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = ps%rho_rad(1:nrho_zc)*ps%rho_rad(1:nrho_zc)

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)*pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%prad,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%prad(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_br,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%prad_br(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_cy,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%prad_cy(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%prad_li,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%prad_li(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
      DEALLOCATE(rho_zc,tmp)

!cj aug-15_2011 pilh==0, pelh (W, electron heating by LH), curlh (A, LH current drive)
#ifdef EZSPL
      nrho_zc = ps%nrho_lhrf-1
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = (ps%rho_lhrf(1:nrho_zc)+ps%rho_lhrf(2:nrho_zc+1))/2.

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

        if(from_xplasma_p .and. .not. use_tsc_lhh) then
           write(6,*) "from_xplasma_p pass pelh", from_xplasma_p, use_tsc_lhh
           write(*,*) "from_xplasma_p pass pelh", from_xplasma_p, use_tsc_lhh
        else
           write(6,*) "writing LH to plasma state", from_xplasma_p, use_tsc_lhh
      CALL EZspline_setup(spln,pa%pelh,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%pelh(1:nrho_zc) = tmp(1:nrho_zc)
        endif

        if( (from_xplasma_j .and. .not. use_tsc_lhh) ) then
           write(6,*) "from_xplasma_j pass curlh", from_xplasma_j, use_tsc_lhh
           write(*,*) "from_xplasma_j pass curlh", from_xplasma_j, use_tsc_lhh
        else
           write(6,*) "writing LH to plasma state", from_xplasma_p, use_tsc_lhh
      CALL EZspline_setup(spln,pa%curlh,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%curlh(1:nrho_zc) = tmp(1:nrho_zc)
        endif

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
      DEALLOCATE(rho_zc,tmp)
#else
!cj mar14_2012 default interpolation
      if( (from_xplasma_j .and. use_tsc_lhh) ) then
      nrho_zc = ps%nrho_lhrf-1
      if(allocated(x)) deallocate(x)
      if(allocated(xh)) deallocate(xh)
      if(allocated(tmp)) deallocate(tmp)
      allocate(x(npsit),stat=ierr)
      allocate(xh(nrho_zc),stat=ierr)
      allocate(tmp(npsit),stat=ierr)
      do i = 1, npsit
         x(i) = float(i-1)/float(npsit-1)
         x(i) = sqrt(x(i))
      enddo
      x(1) = 0.D0
      x(npsit) = 1.D0
      tmp(1:npsit)=pa%pelh(1:npsit)
      call ps_user_rezone1(x,ps%rho_lhrf,tmp, ps%pelh, ierr,   &
                                        nonorm=norm,zonediff=xh)
      tmp(1:npsit)=pa%curlh(1:npsit)
      call ps_user_rezone1(x,ps%rho_lhrf,tmp, ps%curlh, ierr,   &
                           curdens=norm,nonorm=norm,zonediff=xh)
      deallocate(x)
      deallocate(xh)
      deallocate(tmp)
      endif
#endif
!cj mar14_2012 debug
      powtot=0.
      curtot=0.
      do i=1,nrho_zc
         powtot=powtot+ps%pelh(i)
         curtot=curtot+ps%curlh(i)
         !write(*,*) "mypscurlh=", i, ps%pelh(i), ps%curlh(i)
      enddo
      write(*,*) "mypscurlh powtot curtot=", powtot, curtot

!     more zone integrated data
      nrho_zc = ps%nrho
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = ps%rho(1:nrho_zc)*ps%rho(1:nrho_zc)

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)*pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%pohm,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%pohme(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%qie,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%qie(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

!     current also in zone data
      CALL EZspline_setup(spln,pa%curbs,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%curr_bootstrap(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_setup(spln,pa%curoh,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
       ps%curr_ohmic(1:nrho_zc-1) = tmp(2:nrho_zc)-tmp(1:nrho_zc-1)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
      DEALLOCATE(rho_zc,tmp)

!     more step profile data
      nrho_zc = ps%nrho_eq -1
      ALLOCATE(rho_zc(nrho_zc),tmp(nrho_zc))
      rho_zc = (ps%rho_eq(1:nrho_zc) + ps%rho_eq(2:ps%nrho_eq))/2

      CALL EZlinear_init(spln,pa%nrho,ierr)
      CALL EZspline_error(ierr)
      spln%x1(1:pa%nrho) = pa%rho(1:pa%nrho)

      CALL EZspline_setup(spln,pa%pres,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%p_eq(1:nrho_zc) = tmp(1:nrho_zc)
 
      CALL EZspline_setup(spln,pa%qprof,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_interp(spln,nrho_zc,rho_zc,tmp,ierr)
      CALL EZspline_error(ierr)
      ps%q_eq(1:nrho_zc) = tmp(1:nrho_zc)

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
      DEALLOCATE(rho_zc,tmp,tmp2)

!cj jan 12, 2010 coil current
      ps%ncircuits = ncurrentg
      CALL ps_alloc_plasma_state(ierr)
      ps%coil_apt(1:ncurrentg) = csuma(1:ncurrentg)
      write(6,*) "ps_interface: coil current in pa"
!
!     put profiles into plasma state
#ifdef PARTIAL_PLASMA_STATE
!       call PS_RESTORE_HASH_CODES(ierr, state=ps)

      call ps_cclist_remove('*',cclist,ierr) 
      call ps_cclist_add('*',cclist,ierr)
!     call ps_cclist_add('PLASMA',cclist,ierr)
!     call ps_cclist_add('EQ',cclist,ierr)
!     call ps_cclist_add('NBI',cclist,ierr)
!     call ps_cclist_add('IC',cclist,ierr)
!     call ps_cclist_add('LH',cclist,ierr)
!     call ps_cclist_add('EC',cclist,ierr)
!     call ps_cclist_add('RUNAWAY',cclist,ierr)
!     call ps_cclist_add('FUS',cclist,ierr)
!     call ps_cclist_add('RAD',cclist,ierr)
!     call ps_cclist_add('GAS',cclist,ierr)
!     call ps_cclist_add('LMHD',cclist,ierr)
!     call ps_cclist_add('RIPPLE',cclist,ierr)
!     call ps_cclist_add('ANOM',cclist,ierr)

      !  cclist contains the PLASMA component only.
      call ps_update_report(ierr,state=ps,cclist=cclist,               &
                        report_filename='ps_update_state_pa.list')

!     call ps_copy_plasma_state(ps, aux, ierr, cclist=cclist) 
!     call ps_write_update_file(trim(xplasma_filename), ierr, state=ps, &
      call ps_write_update_file(trim(xplasma_partial_ps_filename),      &
                                ierr, state=ps, cclist=cclist)
!     call ps_store_plasma_state(ierr, filename=xplasma_filename, state=ps)
!#else
#endif
!>>>>debug
      CALL ps_store_plasma_state(ierr,xplasma_filename,state=ps)
      if(ierr .ne. 0) then
      write(6,*) "ps_interface: plasma state store error in pa"
      cmd="date > transp_kill.dat"
      ierr = jsystem(trim(cmd))
      stop
      endif
      

      end subroutine put_plasma_state_pa
      END MODULE plasma_state_pa_mod

