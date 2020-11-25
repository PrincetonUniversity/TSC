       MODULE TRTSCP
       IMPLICIT NONE

       REAL*8 :: tr_ts1, tr_ts2
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_pfi, tr_pfe
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_pbi, tr_pbe
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_pilh, tr_pelh
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_piicrf, tr_peicrf
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_piecrf, tr_peecrf
       
       CONTAINS


       subroutine get_dataP_from_xplasma(filenamx,nxin,vary,read_error)
       USE trtsc
       USE plasma_state_mod
       USE sawtooth_mod
       USE clinam, only : kcycle, ineg, npsit

!      use clinam, only : vary
!      use scr3
!      use saprop
!      use wallcl



       implicit none
       INTEGER, intent (in) :: nxin
       REAL*8, intent (in) :: vary(*)
       REAL*8, ALLOCATABLE, DIMENSION(:) :: x, xh 
       REAL*8, ALLOCATABLE, DIMENSION(:) :: pint, pdens
       REAL*8, ALLOCATABLE, DIMENSION(:) :: tmp, tmp2
       CHARACTER*32 :: list_name
       CHARACTER*32, ALLOCATABLE, DIMENSION(:) :: enames
       CHARACTER*80, ALLOCATABLE, DIMENSION(:) :: elabels
       REAL*8, ALLOCATABLE, DIMENSION(:) :: rvals
       INTEGER, ALLOCATABLE, DIMENSION(:) :: itype, ivals
       CHARACTER*(*) :: filenamx

       INTEGER :: ix, iwant, id, list_id, inum, istat
       INTEGER :: it, ii, len_ename
       integer :: icount = 0, read_error
       logical :: smooth = .true. , norm=.true.
      real*8 :: sum1, sum2



       call ps_get_plasma_state(ierr,filenamx)
       if(ierr.ne.0) then
         write(6,*) " in get_dataP_from_xplasma : failed to read ",  trim(filenamx)
         from_xplasma_p=.false.
         from_xplasma=.false.
         ineg=61
         return
       endif

       if(allocated(tr_pfe)) deallocate(tr_pfe)
       if(allocated(tr_pfi)) deallocate(tr_pfi)
       if(allocated(tr_pbe)) deallocate(tr_pbe)
       if(allocated(tr_pbi)) deallocate(tr_pbi)
       if(allocated(tr_pelh)) deallocate(tr_pelh)
       if(allocated(tr_pilh)) deallocate(tr_pilh)
       if(allocated(tr_peicrf)) deallocate(tr_peicrf)
       if(allocated(tr_piicrf)) deallocate(tr_piicrf)
       if(allocated(tr_peecrf)) deallocate(tr_peecrf)
       if(allocated(tr_piecrf)) deallocate(tr_piecrf)
       allocate(tr_pfe(nxin,2),stat=ierr)
       allocate(tr_pfi(nxin,2),stat=ierr)
       allocate(tr_pbe(nxin,2),stat=ierr)
       allocate(tr_pbi(nxin,2),stat=ierr)
       allocate(tr_pelh(nxin,2),stat=ierr)
       allocate(tr_pilh(nxin,2),stat=ierr)
       allocate(tr_peicrf(nxin,2),stat=ierr)
       allocate(tr_piicrf(nxin,2),stat=ierr)
       allocate(tr_peecrf(nxin,2),stat=ierr)
       allocate(tr_piecrf(nxin,2),stat=ierr)
       tr_pfe = 0.0d0
       tr_pfi = 0.0d0
       tr_pbi = 0.0d0
       tr_pbe = 0.0d0
       tr_pelh = 0.0d0
       tr_pilh = 0.0d0
       tr_peicrf = 0.0d0
       tr_piicrf = 0.0d0
       tr_peecrf = 0.0d0
       tr_piecrf = 0.0d0
!
!      define time interval
!      if(.not.allocated(tmp)) allocate(tmp(2),stat=istat)
!      call eq_glistnum('TR_TIMESTEP',id)
!      if(id .ne. 0) then
!      call eq_list_r8vals(id,2,tmp,ierr)
!      tr_ts1 = tmp(1)
!      tr_ts2 = tmp(2)
!      else
!      stop "error in getting times"
!      endif

       if(allocated(x)) deallocate(x)
       if(allocated(xh)) deallocate(xh)
       if(allocated(pint)) deallocate(pint)
       if(allocated(pdens)) deallocate(pdens)
       allocate(x(nxin),stat=ierr)
       allocate(xh(nxin-1),stat=ierr)
       allocate(pint(nxin-1),stat=ierr)
       allocate(pdens(nxin),stat=ierr)

       do ix = 1, nxin
         x(ix) = float(ix-1)/float(nxin-1)
         x(ix) = sqrt(x(ix))
       enddo
       x(1) = 0.D0
       x(nxin) = 1.D0
!
!      xh array is zone volume
       xh(1:nxin-1)=vary(2:nxin)-vary(1:nxin-1)

       iwant=0
       it = 1
!
!      use tsc differential volume and zone total heating
!      to calculate heating density for TSC so as to preserve the
!      total heating calculated by TRANSP
!
!---fusion energy (electron power)
       if(.not. use_tsc_fp) then
         id=ps%id_pfuse
         if(saw_tr) id=saw1%id_pfuse
         call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
         if(ierr .ne. 0) then
           write(6,*) "---> pfuse read error " 
           pint = 0.0
!          from_xplasma_p = .false.
           use_tsc_fp=.true.
           read_error = read_error + 1
         endif
         pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
         pdens(1) = pdens(2)
         do ix = 1, nxin
           tr_pfe(ix,it)=pdens(ix)
         enddo
       endif
!
!---neutral beams (electron power)
       if(.not. use_tsc_beam) then
         id=ps%id_pbe
         if(saw_tr) id=saw1%id_pbe
         call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
         if(ierr .ne. 0) then
!          from_xplasma_p = .false.
!          use_tsc_beam=.true.
!          if(kcycle.gt.0) ineg=61
           read_error = read_error + 1
           write(6,*) "---> pbe read error ", ps%id_pbe, kcycle, npsit, ierr, read_error
           pint = 0.0
!>>>>>debug
         else
           write(6,*) "---> pbe read -ok- ",id,kcycle,npsit
         endif
         pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
         pdens(1) = pdens(2)
         do ix = 1, nxin
           tr_pbe(ix,it)=pdens(ix)
         enddo
!
!.....debug printout
!       write(6,4001) (pdens(ix),ix=1,nxin)
!4001 format(1p10e12.4)

!
!.....debug:  calculate and print total beam power to electrons
         sum1=0.0
         do ix=1, ps%nrho_nbi-1
           sum1=sum1+ps%pbe(ix)
         enddo
         sum2=0.0
         do ix=2, nxin
           sum2=sum2+pdens(ix)*(vary(ix)-vary(ix-1))
         enddo
         write(6,*) "beam power to electrons: ", sum1, sum2
       endif

!
!....ion cyclotron (ICRF) power to electrons
       if(.not. use_tsc_ich) then
         if(allocated(ps%id_picrf_totals)) then
           id=ps%id_picrf_totals(0)
           if(saw_tr) id=saw1%id_picrf_totals(0)
           call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
           if(ierr .ne. 0) then
             write(6,*) "---> picrf_e read error "
             pint = 0.0
!            from_xplasma_p = .false.
             use_tsc_ich=.true.
             read_error = read_error + 1
           endif
           pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
         else
           ierr = 1
           write(6,*) "---> picrf_e read error "
           pint = 0.0
!          from_xplasma_p = .false.
           use_tsc_ich=.true.
           read_error = read_error + 1
         endif

!
         id=ps%id_pmine
         if(saw_tr) id=saw1%id_pmine
         call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
         if(ierr .ne. 0) then
           pint = 0.0
           write(6,*) "---> pmine read error "
!          from_xplasma_p = .false.
           use_tsc_ich=.true.
           read_error = read_error + 1
         endif
         pdens(2:nxin)=pdens(2:nxin)+pint(1:nxin-1)/xh(1:nxin-1)
         pdens(1) = pdens(2)
         do ix = 1, nxin
           tr_peicrf(ix,it)=pdens(ix)
         enddo
         sum1=0.0
         do ix=1, ps%nrho_icrf-1
           sum1=sum1+ps%picrf_totals(ix,0)+ps%pmine(ix)
         enddo
         sum2=0.0
         do ix=2, nxin
           sum2=sum2+pdens(ix)*(vary(ix)-vary(ix-1))
         enddo
         write(6,*) "icrf power to electrons: ", sum1, sum2
       endif
!
!...lower hybrid
       if(.not. use_tsc_lhh) then
          if(allocated(ps%pelh)) then
            id=ps%id_pelh
            if(saw_tr) id=saw1%id_pelh
            call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
            if(ierr .ne. 0) then
              write(6,*) "---> pelh read error "
              pint = 0.0
!             from_xplasma_p = .false.
              use_tsc_lhh=.true.
              read_error = read_error + 1
            endif
            pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
            pdens(1) = pdens(2)
            do ix = 1, nxin
              tr_pelh(ix,it)=pdens(ix)
            enddo
          else
            write(6,*) "---> pelh not exits in ps use_tsc_lhh=", use_tsc_lhh
          endif
!
         sum1=0.0
         do ix=1, ps%nrho_lhrf-1
           sum1=sum1+ps%pelh(ix)
         enddo
         sum2=0.0
         do ix=2, nxin
           sum2=sum2+pdens(ix)*(vary(ix)-vary(ix-1))
         enddo
         write(6,*) "LH power to electrons: ",sum1,sum2
       endif
!
!..electron cyclotron heating
       if(.not. use_tsc_ech) then
         id=ps%id_peech
         if(saw_tr) id=saw1%id_peech
         call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
         if(ierr .ne. 0) then
           write(6,*) "---> peech read error "
           pint = 0.0
!          from_xplasma_p = .false.
           use_tsc_ech=.true.
           read_error = read_error + 1
         endif
         pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
         pdens(1) = pdens(2)
         do ix = 1, nxin
           tr_peecrf(ix,it)=pdens(ix)
         enddo 
         sum1=0.0
         do ix=1, ps%nrho_ecrf-1
           sum1=sum1+ps%peech(ix)
         enddo
         sum2=0.0
         do ix=2, nxin
           sum2=sum2+pdens(ix)*(vary(ix)-vary(ix-1))
         enddo
         write(6,*) "ecrf power to electrons: ", sum1,sum2
       endif

!
!...fusion power (ion power)
       if(.not. use_tsc_fp) then
         id=ps%id_pfusi
         if(saw_tr) id=saw1%id_pfusi
         call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
         if(ierr .ne. 0) then
           write(6,*) "---> pfusi read error "
           pint = 0.0
!          from_xplasma_p = .false.
           use_tsc_fp=.true.
           read_error = read_error + 1
         endif
         pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
         pdens(1) = pdens(2)
         do ix = 1, nxin
           tr_pfi(ix,it)=pdens(ix)
         enddo
       endif
!
!...neutral beam power to ions
       if(.not. use_tsc_beam) then
         id=ps%id_pbi
         if(saw_tr) id=saw1%id_pbi
         call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
         if(ierr .ne. 0) then
!          from_xplasma_p = .false.
           use_tsc_beam=.true.   !fmp. uncommented for consistency with other sources
!          if(kcycle.gt.0) ineg=61
           read_error = read_error + 1
           write(6,*) "---> pbi read error ", id,kcycle,npsit,ierr, read_error
!          pint = 0.0
!>>>>>debug
         else
           write(6,*) "---> pbi read -ok- ",id,kcycle,npsit
         endif
         pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
         pdens(1) = pdens(2)
         do ix = 1, nxin
           tr_pbi(ix,it)=pdens(ix)
         enddo
!
!.....debug printout
!         write(6,4001) (pdens(ix),ix=1,nxin)

         sum1=0.0
         do ix=1, ps%nrho_nbi-1
           sum1=sum1+ps%pbi(ix)
         enddo
         sum2=0.0
         do ix=2, nxin
           sum2=sum2+pdens(ix)*(vary(ix)-vary(ix-1))
         enddo
         write(6,*) "beam power to ion: ", sum1,sum2
       endif
!
!....ich power to ions
       if(.not. use_tsc_ich) then
         id=ps%id_picth
         if(saw_tr) id=saw1%id_picth
         call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
         if(ierr .ne. 0) then
           write(6,*) "---> picth read error "
           pint = 0.0
!          from_xplasma_p = .false.
           use_tsc_ich=.true.
           read_error = read_error + 1
         endif
         pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
         id=ps%id_pmini
         if(saw_tr) id=saw1%id_pmini
         call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
         if(ierr .ne. 0) then
           write(6,*) "---> pmini read error "
           pint = 0.0
!          from_xplasma_p = .false.
           use_tsc_ich=.true.
           read_error = read_error + 1
         endif
         pdens(2:nxin)=pdens(2:nxin)+pint(1:nxin-1)/xh(1:nxin-1)
         pdens(1) = pdens(2)
         do ix = 1, nxin
           tr_piicrf(ix,it)=pdens(ix)
         enddo

         sum1=0.0
         do ix=1, ps%nrho_icrf-1
           sum1=sum1+ps%picth(ix)+ps%pmini(ix)
         enddo
         sum2=0.0
         do ix=2, nxin
           sum2=sum2+pdens(ix)*(vary(ix)-vary(ix-1))
          enddo
         write(6,*) "icrf power to ion: ", sum1,sum2
       endif
!
!....lower hybrid power to ions
       if(.not. use_tsc_lhh) then
         if(allocated(ps%pilh)) then
           id=ps%id_pilh
           if(saw_tr) id=saw1%id_pilh
           call ps_rho_rezone1(x,id,pint,ierr,zonesmoo=smooth,nonorm=norm)
           if(ierr .ne. 0) then
             write(6,*) "---> pilh read error "
             pint = 0.0
!            from_xplasma_p = .false.
             use_tsc_lhh=.true.
             read_error = read_error + 1
           endif
           pdens(2:nxin)=pint(1:nxin-1)/xh(1:nxin-1)
           pdens(1) = pdens(2)
           do ix = 1, nxin
             tr_pilh(ix,it)=pdens(ix)
           enddo 
         else
           write(6,*) "---> pilh not exits in ps use_tsc_lhh=", use_tsc_lhh
         endif
       endif

       tr_pfe(:,2) = tr_pfe(:,1)
       tr_pfi(:,2) = tr_pfi(:,1)
       tr_pbi(:,2) = tr_pbi(:,1)
       tr_pbe(:,2) = tr_pbe(:,1)
       tr_pelh(:,2) = tr_pelh(:,1)
       tr_pilh(:,2) = tr_pilh(:,1)
       tr_peicrf(:,2) = tr_peicrf(:,1)
       tr_piicrf(:,2) = tr_piicrf(:,1)
       tr_peecrf(:,2) = tr_peecrf(:,1)
       tr_piecrf(:,2) = tr_piecrf(:,1)

       if(read_error .gt. 0 .or. (.not. from_xplasma_p) ) then
!      tr_pfe = 0.0d0
!      tr_pfi = 0.0d0
!      tr_pbi = 0.0d0
!      tr_pbe = 0.0d0
!      tr_pelh = 0.0d0
!      tr_pilh = 0.0d0
!      tr_peicrf = 0.0d0
!      tr_piicrf = 0.0d0
!      tr_peecrf = 0.0d0
!      tr_piecrf = 0.0d0
          print *, " >>> Error detected when reading data from PS <<<"
          if(use_tsc_fp)                                                     &
            print *, " >>> Fusion Power, if any, is now from TSC <<<"
          if(use_tsc_beam)                                                     &
            print *, " >>> Beam Power, if any, is  now from TSC <<<"
          if(use_tsc_ich)                                                     &
            print *, " >>> ICRF Power, if any, is  now from TSC <<<"
          if(use_tsc_lhh)                                                     &
            print *, " >>> LHH Power, if any, is  now from TSC <<<"
          if(use_tsc_ech)                                                     &
            print *, " >>> ECH Power, if any, is  now from TSC <<<"
!         from_xplasma_p=.false.
!         from_xplasma=.false.
       endif

       return

       end subroutine get_dataP_from_xplasma


       end MODULE TRTSCP

!============================================================
      subroutine auxheat
!
!.....defines source terms for auxiliary heating and alpha heating
!
      USE trtsc
      USE trtscp
      USE CLINAM
      USE SAPROP
      USE SCR3
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ifrstah,j,l,lsav,ns,jj,jstop
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sumhe,tion,exparg,svdt,zlnzti,rdtcgs
      REAL*8 ddens,tdens,savtotj,fac,pval,ppval,sumne,rval,rpval
      REAL*8 sum7,palphat,auxmax,sum1,uintt,term,termd,bpfeed,zeffx
      REAL*8 facden,ar,f1,expr,f2,ecrit,ebeam,facibm,denom,ffact,f3
      REAL*8 afw,dfw,a1fw,a2fw,bedge,bmag0,akpar,omega,rho,blocal
      REAL*8 andu,alphe,cnorm,psiint,fnorm,befo,term1,term2,term3
      REAL*8 omegb,testop,veth,warg
      REAL*8 sum, aint
!============
!     dimension form(ppsi)
!     dimension seun(ppsi),siun(ppsi),alphaea(ppsi)
      data ifrstah/ 1 /
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: form
      REAL*8, ALLOCATABLE, DIMENSION(:) :: seun
      REAL*8, ALLOCATABLE, DIMENSION(:) :: siun
      REAL*8, ALLOCATABLE, DIMENSION(:) :: alphaea
!============      
      character*256 :: filenamx
      
      integer :: ifirst_write=1
      integer :: npsit_old, read_error_pass=0, read_error=0
      real*8 :: dum, bpowers_i, bpowers_e
      save :: npsit_old, read_error_pass, read_error
      logical :: first_read_beam=.true., first_read_fw=.true.
      logical :: first_read_lhh=.true., first_read_ec=.true.
      logical :: ex_beam=.false., ex_lhh=.false., ex_fw=.false., ex_ec=.false.
      logical :: ex_radi=.false.,ex_totradi=.false.
      real*8 :: tint, rint, pbe1, pbe2, pbi1, pbi2, ple1, ple2
      real*8 :: pli1, pli2, pfe1, pfe2, pfi1, pfi2
      real*8 :: prad1,prad2,psum
      real*8, dimension(:),allocatable :: time_beam, rad_beam, tot_beam
      real*8, dimension(:),allocatable :: time_lhh, rad_lhh,tot_lhh
      real*8, dimension(:),allocatable :: time_fw, rad_fw, tot_pfw
      real*8, dimension(:),allocatable :: time_ec, rad_ec, tot_pec
      real*8, dimension(:),allocatable :: time_radi, rad_radi,tot_radi,facpow
      real*8, dimension(:,:),allocatable :: pbe, pbi, ple, pli
      real*8, dimension(:,:),allocatable :: pfe, pfi
      real*8, dimension(:,:),allocatable :: pee, pei
      real*8, dimension(:,:),allocatable :: pow_radi
      integer :: j1, j2, kk, ii
      integer :: ntime_beam, nrad_beam
      integer :: ntime_lhh, nrad_lhh
      integer :: ntime_fw, nrad_fw
      integer :: ntime_ec, nrad_ec
      integer :: ntime_radi, nrad_radi
      real*8 :: tfluxd
      real*8 :: tfluxn, rhox
!============      
      IF(.not.ALLOCATED(form)) ALLOCATE( form(ppsi), STAT=istat)
      IF(.not.ALLOCATED(seun)) ALLOCATE( seun(ppsi), STAT=istat)
      IF(.not.ALLOCATED(siun)) ALLOCATE( siun(ppsi), STAT=istat)
      IF(.not.ALLOCATED(alphaea)) ALLOCATE( alphaea(ppsi), STAT=istat)
!============      
      IF(ALLOCATED(pow_radi)) DEALLOCATE(pow_radi)
      IF(ALLOCATED(time_radi)) DEALLOCATE(time_radi)
      IF(ALLOCATED(rad_radi)) DEALLOCATE(rad_radi)
      IF(ALLOCATED(tot_radi)) DEALLOCATE(tot_radi)
      IF(ALLOCATED(facpow)) DEALLOCATE(facpow)
!============      
      if (istat .ne. 0) stop 'Allocation Error : auxheat  ' 
!============      
      ifirst_write=0
      if(read_error_pass .eq. 0) read_error=0
!     allowed first pass to go through for initialization
      if(read_error .gt. 1 ) then
      if(read_error_pass .eq. 1 .and. int(acoef(4980)+0.1) .eq. 1) then
        print *,                                                              &
       &  " >>> Error reading PS detected, TSC is asked to stop <<<"
        ineg=61
      endif
      endif
      read_error_pass = read_error_pass + 1

      if(from_xplasma) then
        if( .not. use_tsc_fp .or. .not. use_tsc_beam .or.                     &
            .not. use_tsc_ich .or. .not. use_tsc_lhh .or.                     &
            .not. use_tsc_ech ) then
          filenamx=trim(xplasma_filename)
          if(read_xplasma_p .or. (npsit .ne. npsit_old)) then
            call get_dataP_from_xplasma(filenamx, npsit, vary, read_error)
!
!.........set flags so get_dataP_from_xplasma not called again unless npsit changes
            npsit_old = npsit
            read_xplasma_p = .false.
          endif
        endif
        ifirst_write=1
      endif

      if(ifrstah .eq. 1) then
        ifrstah = 0
        if(irst1.eq.0) then
          toldaux = times
        endif
      endif
!
      do 10 j=1,npsit
        savee(j) = 0._R8
        savei(j) = 0._R8
        savia(j) = 0._R8
        savea(j) = 0._R8
        sraveb(j) = 0._R8
        savelh(j) = 0._R8
        savefw(j) = 0._R8
        savifw(j) = 0._R8
        savebm(j) = 0._R8
        savibm(j) = 0._R8
   10 continue
!
!.(1)...............................................................
!...local alpha heating
!...expressions from   hively,  nuc. fus. vol 17 p873  1977
!
!.....calculate alpha power for all values of ialpha, but include in
!     rest of calculation only for ialpha=1
!
      sumhe=0.0_R8
      do 50 j=1,npsit
      tion=ti(j)*.001_R8
      if(tion.lt.0.5_R8) tion = 0.5_R8
      if(tion.gt.240._R8) tion = 240._R8
      if(tion.gt.80._R8) goto 60
      exparg= -23.836_R8- 22.712_R8*tion**(-0.275_R8) - 9.393E-2_R8*     &  
     & tion                                                              &  
     &     +7.994E-4_R8*tion**2 - 3.144E-6_R8*tion**3
      svdt=exp(exparg)
      goto 61
60    zlnzti=log(tion)
      exparg=-.187_R8*zlnzti**2+1.447_R8*zlnzti-37.43_R8
      svdt=1.027_R8*exp(exparg)
61    continue
!
!.....note:  assumes 50%/50% DT mix, includes dilution due to alpha-ash
      rdtcgs=0.5_R8*(anhy(j))*1.E-6_R8
!
!.....energy due to alpha-particles assuming all heat is
!     deposited at same surface
      ddens = acoef(113)*rdtcgs/.5_R8
      tdens = (1._R8-acoef(113))*rdtcgs/.5_R8
      savtotj =5.6E-7_R8*ddens*tdens*svdt*usdp/usdt
!.....70% to electrons and 30% to ions
      savea(j) = 0.70_R8*savtotj
      savia(j) = 0.30_R8*savtotj
!
!.....number of alpha-particles born on this surface per sec
      alphanu(j) = rdtcgs**2*svdt*vp(j)*dpsi*1.E6_R8
!
!.....add to density of helium ash
      if(acoef(93).le.0._R8) acoef(93)=1.E6_R8
      if(ialpha.ge.1) then
        anheb(j) = anheb(j) + (times-toldaux)*rdtcgs**2*svdt*1.E6_R8       &
                 - (times-toldaux)*anheb(j)/acoef(93)
      endif
      if(anheb(j).le.0._R8) anheb(j) = 0._R8
      sumhe = sumhe + anheb(j)*vp(j)*dpsi
!
!.....define alpha density and alpha pressure here
      alphade(j) = 0._R8
      alphapr(j) = 0._R8
      if(tion.lt.3.7_R8) go to 50
      fac = .029_R8*((anhy(j))/ane(j))**2*(tion-3.7_R8)
      call peval(xsv(j),2,pval,ppval,imag,jmag)
      alphapr(j) = fac*pval
      alphade(j) = anheb(j)
   50 continue

      if(from_xplasma_p .and. .not. use_tsc_fp) then
      fac = usdp/usdt
!     dum = (times-tr_ts1)/(tr_ts2-tr_ts1)
      dum = 0.0
      do j=1,npsit
      savia(j)=tr_pfi(j,1)+dum*(tr_pfi(j,2)-tr_pfi(j,1))
      savea(j)=tr_pfe(j,1)+dum*(tr_pfe(j,2)-tr_pfe(j,1))
      savia(j)=fac*savia(j)
      savea(j)=fac*savea(j)
      enddo
      endif

!
!
!.....thermal helium density redefined to have electron profile shape
!.....NOTE: anhe is thermal helium density, anheb is helium birth
!.....      density, the integrals of these are always assumed equal
      sumne=0.0_R8
      do 52 j=1,npsit
      call reval(xsv(j),idens,isurf,rval,rpval,imag,jmag)
      sumne = sumne + rval*vp(j)*dpsi*udsd
   52 continue
      do 53 j=1,npsit
      call reval(xsv(j),idens,isurf,rval,rpval,imag,jmag)
      anhe(j) = (sumhe/sumne)*rval*udsd
      anhe(j) = 0.25_R8*anheb(j) + 0.75_R8*anhe(j)
  53  continue
      toldaux = times
51    continue
      if(irippl .gt. 0 .and. iskipsf.le.1) call ripple
!
!.....calculate total alpha power
      sum7 = 0._R8
      do 71 j=2,npsit
   71 sum7 = sum7 + (vary(j)-vary(j-1))*(savea(j) + savia(j))
      palphat = sum7*ialpha*udsp/udst
!
!.....limit max auxialliary heating for burn control
      auxmax = acoef(47)*1.E6_R8- (palphat+pohmic)
      if(auxmax .le.0._R8 ) auxmax = 0._R8
      if(acoef(47) .eq. 0) auxmax = 1.E12_R8
!
!
!.....new feedback on DT species ratio added 7/28/94
      if(acoef(114) .le. 0 .or. kcycle.le.0) go to 55
      sum1 = 0._R8
      do 6609 l=2,npsit
      call peval(xsv(l),2,pval,ppval,imag,jmag)
 6609 sum1 = sum1 + pval*vp(l)*dpsi
      uintt = 1.5_R8*sum1*udsi
      term = 1.0_R8*(acoef(114)*1.E6_R8- uintt)
      termd = -acoef(115)*(uintt-uintold)/(dts)
      bpfeed = bpowers+term+termd
      if(bpfeed .lt. 0._R8) bpfeed=0._R8
!     if(acoef(113) .lt. 0.01 .or. acoef(113) .gt. 0.49) go to 54
!     fnt = uintt / (acoef(113)*(1.-acoef(113)))
!     denom = (1. - 2.*acoef(113))*fnt
!     if(denom .le. 0 .or. dts .le. 0) go to 54
!     term = 0.1*( acoef(114)*1.e6
!    1                           - (uintt)) / denom
!     termd = -acoef(115)*(uintt - uintold)/(denom*dts)
!c    tmax = 0.10*acoef(113)
!c    if(abs(term) .gt. tmax) term = sign(tmax,term)
!     acoef(113) = acoef(113) + term + termd
!  54 continue
!     if(acoef(113) .gt. 0.49) acoef(113) = 0.49
!     if(acoef(113) .lt. 0.01) acoef(113) = 0.01
      uintold=uintt
   55 continue
!
!.....end of alpha heating coding
!..................................................................
!##############################################################################
!#                                                                            #
!#               * * *   NOTE   * * *                                         #
!#   SPECIAL FORMS FOR THE SOURCE FUNCTIONS FOR ACOEF(296)=5.                 #
!#                                                                            #
!############################################################################## #
      if(acoef(296) .eq. 5._R8) go to 6600
!
!.(2).bremstrahlung and cyclotron radiation
!
!         formula from NRL Plasma Formulary
!
!                                    (assumes Gaunt factor of 1.0)
      do 81 j=1,npsit
!
      zeffx = zeff
!	bottom line - commented out on Jan, 19th, 2011 (F.M. Poli)
!	The radiated power was underestimated
!      if(iimp.gt.0 .and. zeffb.gt.0._R8) zeffx = zeffb
      savebre(j) = (1.6E-38_R8)*(ane(j))**2*zeffx*sqrt(te(j))*usdp/usdt  &  
     &    * acoef(877)
!
!....temporary fix for hydrogen radiative recombination
!   (suggested by Paul Barks and Brad Merrill  9/01/99)
      if(te(j).le.100._R8) then
      savebre(j) = savebre(j)*(1._R8+33._R8/te(j))
                        endif
!
      if(rmajor.le.0._R8 .or. rminor.le.0._R8) go to 81
!
!....Reference:   Trubnikov, Reviews of Plasma Physics, Vol 7 (1979)
      savecyc(j) = (1.3E-18_R8)*sqrt(ane(j))*(te(j))**2.5_R8             &  
     &             *(gzero/rmajor)**2.5_R8*(1._R8-acoef(55))**0.5_R8/    &  
     & rminor**0.5_R8                                                    &  
     &             *(1._R8+570._R8*rminor/rmajor/te(j)**0.5_R8)**0.5_R8*  &  
     & usdp/usdt
   81 continue
!
!.......................................................................
! read line radiation from external file
!
      if(use_user_radi) then
        inquire(file="user_radi_heating", exist=ex_radi)
        if (ex_radi) then
          open(55,file="user_radi_heating",form="formatted",status="old")
          read(55,*) ntime_radi, nrad_radi
          allocate(time_radi(ntime_radi), rad_radi(nrad_radi))
          allocate(pow_radi(nrad_radi,ntime_radi))
          allocate(tot_radi(ntime_radi))
          read(55,*) time_radi(1:ntime_radi)
          read(55,*) rad_radi(1:nrad_radi)
          read(55,*) tot_radi(1:ntime_radi)    
          read(55,*) pow_radi(1:nrad_radi,1:ntime_radi)
          close (55)
        endif
      endif
      if(use_user_radi .and. iimp.ne.0.0) then
        ineg=21
        write(nout,1011)
 1011 format(" * * * error * * *..acoef(4978).gt.0.0 and iimp.ne.0.0=")
        return
      endif
 
      if(use_user_radi .and. ex_radi) then
        do  kk = 1,ntime_radi-1
          if(times .gt. time_radi(kk) .and. times .le. time_radi(kk+1)) then
             j1 = kk
             j2 = kk + 1
             tint = (times-time_radi(j1))/(time_radi(j2)-time_radi(j1))
             exit
          endif
        enddo
        if(times .le. time_radi(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif
        if(times .gt. time_radi(ntime_radi)) then
           j1 = ntime_radi-1
           j2 = ntime_radi
           tint = 1.
        endif   
      
        sum = 0._R8
        fac = 0._R8
        do j=2,npsit
           tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
           do ii=1,nrad_radi-1
              if(tfluxd .le. rad_radi(1)) then
                 prad1 = pow_radi(1,j1)
                 prad2 = pow_radi(1,j2)
              endif
              if(tfluxd .gt. rad_radi(ii) .and. tfluxd .le. rad_radi(ii+1)) then
                 rint = (tfluxd-rad_radi(ii))/(rad_radi(ii+1)-rad_radi(ii))
                 prad1 = pow_radi(ii,j1) + (pow_radi(ii+1,j1)-pow_radi(ii,j1))*rint
                 prad2 = pow_radi(ii,j2) + (pow_radi(ii+1,j2)-pow_radi(ii,j2))*rint
                 exit
              endif
           enddo             
           sradion(j) = prad1 + (prad2-prad1)*tint
           sum = sum + sradion(j)*(vary(j)-vary(j-1))
        enddo
        sradion(1)=sradion(2)
        fac = tot_radi(j1)+(tot_radi(j2)-tot_radi(j1))*tint
        if (sum .gt. 0._R8) fac = fac/sum
        do j=1,npsit
           sradion(j) = sradion(j)*fac
        end do
        go to 8131 
      endif      
!	end of radiation data
!----------------------------------------------------------------------

 8131 continue


!
!.......................................................................
!.(3).auxiliary heating from beamp array
      do 20 l=1,ntpts
!cj   write(*,*) "--- auxheat L694 beamp(l) =", l, beamp(l)
        if(beamp(l) .gt. 0._R8) go to 30
   20 continue
      go to 512
   30 continue

      if(from_xplasma_p .and. .not. use_tsc_beam) then
!     dum = (times-tr_ts1)/(tr_ts2-tr_ts1)
        dum=0.0
        do j=1,npsit
          savei(j)=tr_pbi(j,1)+dum*(tr_pbi(j,2)-tr_pbi(j,1))
          savee(j)=tr_pbe(j,1)+dum*(tr_pbe(j,2)-tr_pbe(j,1))
        enddo
        bpowers_i = 0.0
        bpowers_e = 0.0
        do j=1,npsit-1
          bpowers_i=bpowers_i+savei(j+1)*(vary(j+1)-vary(j))
          bpowers_e=bpowers_e+savee(j+1)*(vary(j+1)-vary(j))
        enddo
        fac = usdp/usdt
        do j=1,npsit
          savei(j)=fac*savei(j)
          savee(j)=fac*savee(j)
        enddo
        bpowers = bpowers_i+bpowers_e
        dum = bpowers
        if(acoef(4973) .gt. 0.0) then
          if(acoef(114) .ne. 0.0 .and. bpfeed .lt. bpowers)             &
     &              bpowers=bpfeed
          if(bpowers .gt. auxmax) bpowers = auxmax
        endif
        facden = 0.0
        if(ebeamkev.ne.0)                                               &
     &      facden = (bpowers/(ebeamkev*1.6E-16))*usdd/usdt
        if(dum.gt.0.0d0) then
          do j=2, npsit
            sraveb(j) = facden*(savei(j)+savee(j))/fac/(bpowers_i+bpowers_e)
          enddo
          sraveb(1) = sraveb(2)
          fac = bpowers/dum
          do j=2, npsit
            savei(j)=fac*savei(j)
            savee(j)=fac*savee(j)
          enddo
          savei(1)=savei(2)
          savee(1)=savee(2)
        endif

!     if(from_xplasma) then
!     if(ifirst_write .eq. 1 ) then
!        write(nterm,*) "bpowers, savei, savee, tr_pbi, tr_pbe start", bpowers
!!       dvol=0.0
!        do j=1,npsit
!!          if(j.gt.1) dvol=vary(j)-vary(j-1)
!           write(nterm,"(i5,1p6e14.5)") j, savei(j), savee(j), &
!                 tr_pbi(j,1), tr_pbi(j,2), tr_pbe(j,1), tr_pbe(j,2)
!        enddo
!        write(nterm,*) "bpowers, savei, savee, tr_pbi, tr_pbe end", bpowers
!     endif
!     endif
        go to 512
      endif

      if(use_user_beam .and. first_read_beam) then
        inquire(file="user_beam_heating", exist=ex_beam)
        if (ex_beam) then
          open(55,file="user_beam_heating",form="formatted",status="old")
          read(55,*) ntime_beam, nrad_beam
          allocate(time_beam(ntime_beam), rad_beam(nrad_beam))
          allocate(tot_beam(ntime_beam))
          allocate(pbe(nrad_beam,ntime_beam),pbi(nrad_beam,ntime_beam))
          read(55,*) time_beam(1:ntime_beam)
          read(55,*) rad_beam(1:nrad_beam)
          read(55,*) tot_beam(1:ntime_beam)
          read(55,*) pbe(1:nrad_beam,1:ntime_beam)
          read(55,*) pbi(1:nrad_beam,1:ntime_beam)
          close (55)
        endif
        first_read_beam=.false.
      endif

      if(use_user_beam .and. ex_beam) then
        do  kk = 1,ntime_beam-1
          if(times .gt. time_beam(kk) .and. times .le. time_beam(kk+1)) then
             j1 = kk
             j2 = kk + 1
             tint = (times-time_beam(j1))/(time_beam(j2)-time_beam(j1))
             exit
          endif
        enddo
        if(times .le. time_beam(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif
        if(times .gt. time_beam(ntime_beam)) then
           j1 = ntime_beam-1
           j2 = ntime_beam
           tint = 1.
        endif
       
        do l=1,ntpts
           if(tpro(l).le.time .and. tpro(l+1).gt.time) lsav=l 
        end do 
        fac =beamp(lsav)*usdp/usdt
        if (acoef(4990) .GT. 0) fac = tot_beam(j1)*usdp/usdt
        bpowers_e = 0._R8
        bpowers_i = 0._R8
        do j = 2, npsit
           tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
           do  ii=1,nrad_beam-1
              if(tfluxd .gt. rad_beam(ii) .and.                         &
     &                   tfluxd .le. rad_beam(ii+1)) then
                 rint=(tfluxd-rad_beam(ii))/(rad_beam(ii+1)-rad_beam(ii))
                 pbe1 = pbe(ii,j1) + (pbe(ii+1,j1)-pbe(ii,j1))*rint
                 pbe2 = pbe(ii,j2) + (pbe(ii+1,j2)-pbe(ii,j2))*rint
                 pbi1 = pbi(ii,j1) + (pbi(ii+1,j1)-pbi(ii,j1))*rint
                 pbi2 = pbi(ii,j2) + (pbi(ii+1,j2)-pbi(ii,j2))*rint
                 savee(j) = (pbe1 + (pbe2-pbe1)*tint)
                 savei(j) = (pbi1 + (pbi2-pbi1)*tint)
                 bpowers_i=bpowers_i+savei(j)*(vary(j)-vary(j-1))
                 bpowers_e=bpowers_e+savee(j)*(vary(j)-vary(j-1))
                 exit
              endif
           enddo
           if(tfluxd .le. rad_beam(1)) then
              savee(j)=(pbe(1,j1)+(pbe(1,j2)-pbe(1,j1))*tint)
              savei(j)=(pbi(1,j1)+(pbi(1,j2)-pbi(1,j1))*tint)
           endif
        enddo
        savee(1)=savee(2)
        savei(1)=savei(2)
        bpowers = bpowers_e+bpowers_i
        if (bpowers .gt. 0._R8) fac = fac/bpowers
        do j = 1,npsit
           savee(j) = savee(j)*fac
           savei(j) = savei(j)*fac
        end do
        go to 512
      endif

      
      do 400 l=1,ntpts-1
        lsav = l
        if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 410
  400 continue
  410 continue
      bpowers = beamp(lsav)
      if(acoef(4973) .gt. 0.0) then
        if(acoef(114) .ne. 0._R8.and. bpfeed .lt. beamp(lsav))          &
     &             bpowers=bpfeed
        if(bpowers .gt. auxmax) bpowers = auxmax
      endif
      fac = bpowers*usdp/usdt
      facden = 0._R8
      if(ebeamkev.ne.0)                                                 &
     &       facden = (bpowers/(ebeamkev*1.6E-16_R8))*usdd/usdt
      sum = 0._R8
      do 510 j=2,npsit
        ar = (xsv(j)-psimin)/(psilim-psimin)
        if(abs(ar).ge.1.0_R8) go to 508
        f1 = dbeam**2/((ar-abeam)**2 + dbeam**2)
!     ITER off-axis heating
        if( int(acoef(4975)+0.1) .eq. 1 ) then
          f1 = (ar**1.5)*dbeam**2/((ar-abeam)**2 + dbeam**2)
        endif
        expr = 1._R8
        if(nebeam.ne.0) expr=nebeam
        f2 = (1._R8-ar**2)**expr
        form(j) = f1*f2
!     NSTX and KSTAR
!     f1 = dbeam*(1.-ar)**abeam
!     f2 = (1.-dbeam)*(1.-ar**nebeam)
!     form(j) = f1+f2
          go to 509
  508     form(j) = 0._R8
  509     continue
          sum = sum + form(j)*(vary(j)-vary(j-1))
  510   continue
        do 511 j=2,npsit
          ecrit = 14.8_R8*ambeam*te(j)*zeffa2(j)
          ebeam = ebeamkev*1000._R8
          facibm = (0.33_R8*ecrit/ebeam)*(log((ecrit-sqrt(ecrit*ebeam)       &  
     &           + ebeam)/(ecrit+2._R8*sqrt(ecrit*ebeam)+ebeam))+2._R8*sqrt(3._R8)*  &  
     &             atan((2._R8*sqrt(ebeam)-sqrt(ecrit))/sqrt(3._R8*ecrit))+           &  
     &             (sqrt(3._R8)*pi/3._R8))
          savei(j) = facibm*fac*form(j)/sum
          savee(j) = (1._R8-facibm)*fac*form(j)/sum
          sraveb(j) = facden*form(j)/sum
  511   continue
!
!......................................................................
!.(4).power from lower hybrid
!
  512 continue
      IF ( ilhcd .LT. 1) GOTO 612
      if(ifk.gt.0) go to 700
!
!.....use profiles from input file
      do 598 l=1,ntpts
         if(plhamp(l).gt.0) go to 599
  598 continue
      go to 612
  599 continue

      if(use_user_lhh .and. first_read_lhh) then
        inquire(file="user_lhh_heating", exist=ex_lhh)
        if (ex_lhh) then
          open(55,file="user_lhh_heating",form="formatted",status="old")
          read(55,*) ntime_lhh, nrad_lhh
          allocate(time_lhh(ntime_lhh), rad_lhh(nrad_lhh))
          allocate(tot_lhh(ntime_lhh))
          allocate(ple(nrad_lhh,ntime_lhh),pli(nrad_lhh,ntime_lhh))
          read(55,*) time_lhh(1:ntime_lhh)
          read(55,*) rad_lhh(1:nrad_lhh)
          read(55,*) tot_lhh(1:ntime_lhh)
          read(55,*) ple(1:nrad_lhh,1:ntime_lhh)
          read(55,*) pli(1:nrad_lhh,1:ntime_lhh)
          close (55)
        endif
        first_read_lhh=.false.
      endif

      if(use_user_lhh .and. ex_lhh) then
        do  kk = 1,ntime_lhh-1
          if(times .gt. time_lhh(kk) .and. times .le. time_lhh(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_lhh(j1))/(time_lhh(j2)-time_lhh(j1))
            exit
          endif
        enddo
        if(times .le. time_lhh(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif
        if(times .gt. time_lhh(ntime_beam)) then
          j1 = ntime_lhh-1
          j2 = ntime_lhh
          tint = 1.
        endif

        do l=1,ntpts
           if(tpro(l).le.time .and. tpro(l+1).gt.time) lsav=l 
        end do 
        fac = plhamp(lsav)*usdp/usdt
        if (acoef(4990) .GT. 0) fac = tot_lhh(j1)*usdp/usdt
        bpowers_e = 0._R8
        bpowers_i = 0._R8 
        do j = 2, npsit
           tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
           do  ii=1,nrad_lhh-1
              if(tfluxd .gt. rad_lhh(ii) .and.                          &
     &                  tfluxd .le. rad_lhh(ii+1)) then
                 rint = (tfluxd-rad_lhh(ii))/(rad_lhh(ii+1)-rad_lhh(ii))
                 ple1 = ple(ii,j1) + (ple(ii+1,j1)-ple(ii,j1))*rint
                 ple2 = ple(ii,j2) + (ple(ii+1,j2)-ple(ii,j2))*rint
                 pli1 = pli(ii,j1) + (pli(ii+1,j1)-pli(ii,j1))*rint
                 pli2 = pli(ii,j2) + (pli(ii+1,j2)-pli(ii,j2))*rint
!fmp             savelh(j) = (pbe1 + (pbe2-pbe1)*tint)*usdp/usdt
!fmp             savilh(j) = (pbi1 + (pbi2-pbi1)*tint)*usdp/usdt
                 savelh(j) = (ple1 + (ple2-ple1)*tint)
                 savilh(j) = (pli1 + (pli2-pli1)*tint)
                 bpowers_i = bpowers_i+savilh(j)*(vary(j)-vary(j-1))
                 bpowers_e = bpowers_e+savelh(j)*(vary(j)-vary(j-1))
                 exit
               endif
            enddo
            if(tfluxd .le. rad_lhh(1)) then
!fmp           savelh(j) = (ple(1,j1)+(ple(1,j2)-ple(1,j1))*tint)*usdp/usdt
!fmp           savilh(j) = (pli(1,j1) + (pli(1,j2)-pli(1,j1))*tint)*usdp/usdt
               savelh(j) = (ple(1,j1)+(ple(1,j2)-ple(1,j1))*tint)
               savilh(j) = (pli(1,j1) + (pli(1,j2)-pli(1,j1))*tint)
            endif
         enddo
         savelh(1)=savelh(2)
         savilh(1)=savilh(2)
         bpowers = bpowers_e+bpowers_i
         if (bpowers .gt. 0._R8) fac = fac/bpowers
         do j = 1,npsit
            savelh(j) = savelh(j)*fac
            savilh(j) = savilh(j)*fac
         end do
         go to 612
         
      endif

      
      if(from_xplasma_p .and. .not. use_tsc_lhh) then
        fac = usdp/usdt
!     dum=(times-tr_ts1)/(tr_ts2-tr_ts1)
        dum=0.0
        do j=1,npsit
          savelh(j)=tr_pelh(j,1)+dum*(tr_pelh(j,2)-tr_pelh(j,1))
          savilh(j)=tr_pilh(j,1)+dum*(tr_pilh(j,2)-tr_pilh(j,1))
        enddo
        bpowers_i = 0.0
        bpowers_e = 0.0
        do j=1,npsit-1
          bpowers_i=bpowers_i+savelh(j+1)*(vary(j+1)-vary(j))
          bpowers_e=bpowers_e+savilh(j+1)*(vary(j+1)-vary(j))
        enddo
        bpowers=bpowers_i+bpowers_e
        do j=2, npsit
          savelh(j)=fac*savelh(j)
          savilh(j)=fac*savilh(j)
        enddo
        savelh(1)=savelh(2)
        savilh(1)=savilh(2)
        go to 612
      endif

!

      do 600 l=1,ntpts-1
        lsav = l
        if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 607
  600 continue
  607 continue
!.....linear time point interpolation
      denom = tpro(lsav+1)-tpro(lsav)
      if(denom.eq.0) denom = 1._R8
      ffact = (time-tpro(lsav))/denom
      bpowers = plhamp(lsav)
      alh = alhd(lsav) + ffact * (alhd(lsav+1) - alhd(lsav))
      dlh = dlhd(lsav) + ffact * (dlhd(lsav+1) - dlhd(lsav))
      a1lh = a1lhd(lsav) + ffact * (a1lhd(lsav+1) - a1lhd(lsav))
      a2lh = a2lhd(lsav) + ffact * (a2lhd(lsav+1) - a2lhd(lsav))
      fac = bpowers*usdp/usdt
      sum = 0._R8
      do 610 j=2,npsit
        ar = (xsv(j)-psimin)/(psilim-psimin)
        if(ar.le.0) go to 608
        if(abs(ar).ge.1.0_R8) go to 608
        f1 = dlh**2/((ar-alh)**2 + dlh**2)
        f2 = 1._R8
        if(a1lh .gt. 0._R8) f2 = ar**a1lh
        f3 = 1._R8
        if(a2lh .gt. 0._R8) f3 = (1._R8-ar)**a2lh
        form(j) = f1*f2*f3
        go to 609
  608   form(j) = 0._R8
  609   continue
        sum = sum + form(j)*(vary(j)-vary(j-1))
  610 continue
      do 611 j=2,npsit
        savelh(j) = fac*form(j)/sum
  611 continue

!.....use values read from disk file tscouta
!                          powtsc is in watts/cc
  700 continue
      do 711 j=2,npsit
      savelh(j) = powtsc(j)*1.E6_R8*usdp/usdt
      if (nspc.lt.1) go to 711
      savilh(j) = 0._R8
      do 710 ns=1,nspc
      savilh(j) = savilh(j) + powtsci(ns,j)*1.E6_R8*usdp/usdt
  710 continue
  711 continue
  612 continue
  
!
!.(5).power from fast wave icrh
!
      IF ( iicrh .LT. 1) GOTO 812
!.....use profiles from input file
       do 898 l=1,ntpts
          if(picrh(l).gt.0) go to 899
  898 continue
       go to 812
  899 continue
      
      if(use_user_fw .and. first_read_fw) then
         inquire(file="user_fw_heating", exist=ex_fw)
         if (ex_fw) then
           open(55,file="user_fw_heating",form="formatted",status="old")
           read(55,*) ntime_fw, nrad_fw
           allocate(time_fw(ntime_fw), rad_fw(nrad_fw))
           allocate(pfe(nrad_fw,ntime_fw),pfi(nrad_fw,ntime_fw))
           allocate(tot_pfw(ntime_fw))
           read(55,*) time_fw(1:ntime_fw)
           read(55,*) rad_fw(1:nrad_fw)
           read(55,*) tot_pfw(1:ntime_fw)    
           read(55,*) pfe(1:nrad_fw,1:ntime_fw)
           read(55,*) pfi(1:nrad_fw,1:ntime_fw)
           close (55)
         endif
         first_read_fw=.false.
      endif

      if(use_user_fw .and. ex_fw) then
        do  kk = 1,ntime_fw-1
           if(times .gt. time_fw(kk) .and. times .le. time_fw(kk+1)) then
              j1 = kk
              j2 = kk + 1
              tint = (times-time_fw(j1))/(time_fw(j2)-time_fw(j1))
              exit
           endif
        enddo
        if(times .le. time_fw(1)) then
           j1 = 1
           j2 = 2
           tint = 0.
        endif
        if(times .gt. time_fw(ntime_fw)) then
           j1 = ntime_fw-1
           j2 = ntime_fw
           tint = 1.
        endif
 
        do l=1,ntpts
          if(tpro(l).le.time .and. tpro(l+1).gt.time) lsav=l
        end do
        fac = picrh(lsav)*usdp/usdt
        if (acoef(4990) .GT. 0) fac = tot_pfw(j1)*usdp/usdt
        bpowers_e = 0._R8
        bpowers_i = 0._R8
        do j = 2, npsit
          tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
          do  ii=1,nrad_fw-1
            if(tfluxd .gt. rad_fw(ii) .and. tfluxd .le. rad_fw(ii+1)) then
               rint = (tfluxd-rad_fw(ii))/(rad_fw(ii+1)-rad_fw(ii))
               pfe1 = pfe(ii,j1) + (pfe(ii+1,j1)-pfe(ii,j1))*rint
               pfe2 = pfe(ii,j2) + (pfe(ii+1,j2)-pfe(ii,j2))*rint
               pfi1 = pfi(ii,j1) + (pfi(ii+1,j1)-pfi(ii,j1))*rint
               pfi2 = pfi(ii,j2) + (pfi(ii+1,j2)-pfi(ii,j2))*rint
               savefw(j) = (pfe1 + (pfe2-pfe1)*tint)
               savifw(j) = (pfi1 + (pfi2-pfi1)*tint)
               bpowers_i=bpowers_i+savifw(j)*(vary(j)-vary(j-1))
               bpowers_e=bpowers_e+savefw(j)*(vary(j)-vary(j-1))
               exit
            endif
          enddo
          if(tfluxd .le. rad_fw(1)) then
             savefw(j) = (pfe(1,j1) + (pfe(1,j2)-pfe(1,j1))*tint)
             savifw(j) = (pfi(1,j1) + (pfi(1,j2)-pfi(1,j1))*tint)
          endif
        enddo
        savefw(1)=savefw(2)
        savifw(1)=savifw(2)
        bpowers = bpowers_e+bpowers_i
        if (bpowers .gt. 0._R8) fac = fac/bpowers
        do j = 1,npsit
           savefw(j) = savefw(j)*fac
           savifw(j) = savifw(j)*fac
        end do
        go to 812
      endif

!     if(.not. from_xplasma_p .or. use_tsc_ich) then
      if((.not. use_user_fw .and. use_tsc_ich) .or.                     &
     &      (use_user_fw .and. .not. ex_fw)) then
        do 800 l=1,ntpts-1
          lsav = l
          if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 807
  800   continue
  807   continue
!.....linear time point interpolation
        denom = tpro(lsav+1)-tpro(lsav)
        if(denom.eq.0) denom = 1._R8
        ffact = (time-tpro(lsav))/denom
        bpowers = picrh(lsav)
        if( acoef(4974) .gt. 0.0 ) then
          if(acoef(114) .ne. 0. .and. bpfeed .lt. picrh(lsav))              &
     &        bpowers=bpfeed
          if(bpowers .gt. auxmax) bpowers = auxmax
        endif
        afw = afwd(lsav) + ffact * (afwd(lsav+1) - afwd(lsav))
        dfw = dfwd(lsav) + ffact * (dfwd(lsav+1) - dfwd(lsav))
        a1fw = a1fwd(lsav) + ffact * (a1fwd(lsav+1) - a1fwd(lsav))
        a2fw = a2fwd(lsav) + ffact * (a2fwd(lsav+1) - a2fwd(lsav))
        fac = bpowers*usdp/usdt
        sum = 0._R8
        do 810 j=2,npsit
          if(int(acoef(4976)+0.1) .eq. 0) then
            ar = (xsv(j)-psimin)/(psilim-psimin)
          else
            tfluxd = float(j-1)*dpsi
            tfluxn = tfluxd/(float(npsit-1)*dpsi)
            rhox = sqrt(tfluxn)
            ar = rhox 
          endif
          if(ar.le.0) go to 808
          if(abs(ar).ge.1.0_R8) go to 808
          f1 = 1.0
          if(int(acoef(4976)+0.1) .eq. 0) then
            f1 = dfw**2/((ar-afw)**2 + dfw**2)
          endif
          f2 = 1._R8
          if(int(acoef(4976)+0.1) .eq. 0) then
            if(a1fw .gt. 0) f2 = ar**a1fw
          endif
          f3 = 1._R8
          if(int(acoef(4976)+0.1) .eq. 0) then
            if(a2fw .gt. 0) f3 = (1._R8-ar)**a2fw
          else
            if(a2fw .gt. 0 .and. a1fw .gt. 0) f3 = (1.-ar**a1fw)**a2fw
            if(a2fw .gt. 0 .and. a1fw .le. 0) f3 = (1.-ar)**a2fw
          endif
          form(j) = f1*f2*f3
          go to 809
  808     form(j) = 0._R8
  809     continue
          sum = sum + form(j)*(vary(j)-vary(j-1))
  810   continue
        do 811 j=2,npsit
!     savefw(j) = fac*form(j)/sum
          savefw(j) = acoef(4971)*fac*form(j)/sum
          savifw(j) = (1.0-acoef(4971))*fac*form(j)/sum
!     CMOD
!     savefw(j) = 0.55*fac*form(j)/sum
!     savifw(j) = 0.45*fac*form(j)/sum
  811   continue
      endif
      if(from_xplasma_p .and. .not. use_tsc_ich) then
!     dum=(times-tr_ts1)/(tr_ts2-tr_ts1)
        dum=0.0
        do j=2,npsit
          savefw(j)=tr_peicrf(j,1)+dum*(tr_peicrf(j,2)-tr_peicrf(j,1))
!     use savifw for storing icrf direct heating to ions
          savifw(j)=tr_piicrf(j,1)+dum*(tr_piicrf(j,2)-tr_piicrf(j,1))
        enddo
        bpowers_i = 0.0
        bpowers_e = 0.0
        do j=1,npsit-1
          bpowers_i=bpowers_i+savifw(j+1)*(vary(j+1)-vary(j))
          bpowers_e=bpowers_e+savefw(j+1)*(vary(j+1)-vary(j))
        enddo
        bpowers=bpowers_i+bpowers_e
        dum = bpowers
        if(acoef(4974) .gt. 0.0) then
          if(acoef(114) .ne. 0. .and. bpfeed .lt. bpowers)              &
     &             bpowers=bpfeed
          if(bpowers .gt. auxmax) bpowers = auxmax
        endif

        fac = usdp/usdt
        if(dum .gt. 0.0) fac = fac*bpowers/dum
        do j=2, npsit
          savefw(j)=fac*savefw(j)
          savifw(j)=fac*savifw(j)
        enddo
        savefw(1)=savefw(2)
        savifw(1)=savifw(2)
      endif

!cj...      if(from_xplasma_p .and. use_tsc_ech) then
!cj...      fac = usdp/usdt
!cj...!     dum = (times-tr_ts1)/(tr_ts2-tr_ts1)
!cj...      dum = 0.0
!cj...!     use savebm for storing ech heating to electrons
!cj...      do j=2,npsit
!cj...      savebm(j)=tr_peecrf(j,1)+dum*(tr_peecrf(j,2)-tr_peecrf(j,1))
!cj...      savebm(j)=fac*savebm(j)
!cj...      enddo
!cj...      savebm(1) = savebm(2)
!cj...      endif
  812 continue
!
!.(6).power from ecrh
!
      IF ( iecrh .LT. 1) GOTO 912
!
!.....use profiles from input file
      do 998 l=1,ntpts
         if(pecrh(l).gt.0) go to 999
  998 continue
      go to 912
  999 continue

      if( use_user_ech .and. first_read_ec) then
         inquire(file="user_ec_heating", exist=ex_ec)
         if (ex_ec) then
           open(66,file="user_ec_heating",form="formatted",status="old")
           read(66,*) ntime_ec, nrad_ec
!cj      allocate(time_ec(ntime_beam), rad_ec(nrad_beam))
           allocate(time_ec(ntime_ec), rad_ec(nrad_ec), tot_pec(ntime_ec))
           allocate(pee(nrad_ec,ntime_ec),pei(nrad_ec,ntime_ec))
           read(66,*) time_ec(1:ntime_ec)
           read(66,*) rad_ec(1:nrad_ec)
           read(66,*) tot_pec(1:ntime_ec)
           read(66,*) pee(1:nrad_ec,1:ntime_ec)
           read(66,*) pei(1:nrad_ec,1:ntime_ec)
           close (66)
         endif
         first_read_ec=.false.
      endif

      if( use_user_ech .and. ex_ec) then
        do  kk = 1,ntime_ec-1
           if(times .gt. time_ec(kk) .and. times .le. time_ec(kk+1)) then
             j1 = kk
             j2 = kk + 1
             tint = (times-time_ec(j1))/(time_ec(j2)-time_ec(j1))
             exit
           endif
        enddo
        if(times .le. time_ec(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif
        if(times .gt. time_ec(ntime_ec)) then
           j1 = ntime_ec-1
           j2 = ntime_ec
           tint = 1.
        endif
        
        do l=1,ntpts
           if(tpro(l).le.time .and. tpro(l+1).gt.time) lsav=l   
        end do
        fac = pecrh(lsav)*usdp/usdt
        if (acoef(4990) .GT. 0) fac = tot_pec(j1)*usdp/usdt
        bpowers_e = 0._R8
        bpowers_i = 0._R8
        do j = 2, npsit
           tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
           do  ii=1,nrad_ec-1
              if(tfluxd .gt. rad_ec(ii) .and. tfluxd .le. rad_ec(ii+1)) then
                 rint = (tfluxd-rad_ec(ii))/(rad_ec(ii+1)-rad_ec(ii))
                 pfe1 = pee(ii,j1) + (pee(ii+1,j1)-pee(ii,j1))*rint
                 pfe2 = pee(ii,j2) + (pee(ii+1,j2)-pee(ii,j2))*rint
!cj      pfi1 = pfi(ii,j1) + (pfi(ii+1,j1)-pfi(ii,j1))*rint
!cj      pfi2 = pfi(ii,j2) + (pfi(ii+1,j2)-pfi(ii,j2))*rint
                 pfi1 = pei(ii,j1) + (pei(ii+1,j1)-pei(ii,j1))*rint
                 pfi2 = pei(ii,j2) + (pei(ii+1,j2)-pei(ii,j2))*rint
!fmp             savebm(j) = (pfe1 + (pfe2-pfe1)*tint)*usdp/usdt
!fmp             savibm(j) = (pfi1 + (pfi2-pfi1)*tint)*usdp/usdt
                 savebm(j) = (pfe1 + (pfe2-pfe1)*tint)
                 savibm(j) = (pfi1 + (pfi2-pfi1)*tint)
                 bpowers_e = bpowers_e + savebm(j)*(vary(j)-vary(j-1))
                 bpowers_i = bpowers_i + savibm(j)*(vary(j)-vary(j-1))
                 exit
              endif
           enddo
           if(tfluxd .le. rad_ec(1)) then
!fmp          savebm(j) = (pee(1,j1) + (pee(1,j2)-pee(1,j1))*tint)*usdp/usdt
!fmp          savibm(j) = (pei(1,j1) + (pei(1,j2)-pei(1,j1))*tint)*usdp/usdt
              savebm(j) = (pee(1,j1) + (pee(1,j2)-pee(1,j1))*tint)
              savibm(j) = (pei(1,j1) + (pei(1,j2)-pei(1,j1))*tint)
           endif
        enddo
        savebm(1)=savebm(2)
        savibm(1)=savibm(2)
        bpowers = bpowers_e+bpowers_i
        if (bpowers .gt. 0._R8) fac = fac/bpowers
        do j = 1,npsit
           savebm(j) = savebm(j)*fac
           savibm(j) = savibm(j)*fac
        end do
        go to 912
      endif

!     if(.not. from_xplasma_p .or. use_tsc_ech) then
      if(.not. use_user_ech .and. use_tsc_ech) then
        do 900 l=1,ntpts-1
          lsav = l
          if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 907
  900   continue
  907   continue
!.....linear time point interpolation
        denom = tpro(lsav+1)-tpro(lsav)
        if(denom.eq.0) denom = 1._R8
        ffact = (time-tpro(lsav))/denom
        bpowers = pecrh(lsav)
        if( acoef(4974) .gt. 0.0 ) then
          if(acoef(114) .ne. 0. .and. bpfeed .lt. pecrh(lsav))              &
     &         bpowers=bpfeed
          if(bpowers .gt. auxmax) bpowers = auxmax
        endif 
        afw = aecd(lsav) + ffact * (aecd(lsav+1) - aecd(lsav))
        dfw = decd(lsav) + ffact * (decd(lsav+1) - decd(lsav))
        a1fw = a1ecd(lsav) + ffact * (a1ecd(lsav+1) - a1ecd(lsav))
        a2fw = a2ecd(lsav) + ffact * (a2ecd(lsav+1) - a2ecd(lsav)) 
        fac = bpowers*usdp/usdt
        sum = 0._R8
        do 910 j=2,npsit
          if(int(acoef(4976)+0.1) .eq. 0) then
             ar = (xsv(j)-psimin)/(psilim-psimin)
          else
            tfluxd = float(j-1)*dpsi
            tfluxn = tfluxd/(float(npsit-1)*dpsi)
            rhox = sqrt(tfluxn)
            ar = rhox 
          endif
          if(ar.le.0) go to 908
          if(abs(ar).ge.1.0_R8) go to 908
          f1 = 1.0
          if(int(acoef(4976)+0.1) .eq. 0) then
            f1 = dfw**2/((ar-afw)**2 + dfw**2)
          endif
          f2 = 1._R8
          if(int(acoef(4976)+0.1) .eq. 0) then
            if(a1fw .gt. 0) f2 = ar**a1fw
          endif
          f3 = 1._R8
          if(int(acoef(4976)+0.1) .eq. 0) then
            if(a2fw .gt. 0) f3 = (1._R8-ar)**a2fw
          else
            if(a2fw .gt. 0 .and. a1fw .gt. 0) f3 = (1.-ar**a1fw)**a2fw
            if(a2fw .gt. 0 .and. a1fw .le. 0) f3 = (1.-ar)**a2fw
          endif
          form(j) = f1*f2*f3
          go to 909
  908     form(j) = 0._R8
  909     continue
          sum = sum + form(j)*(vary(j)-vary(j-1))
  910   continue
        do 911 j=2,npsit
!     savefw(j) = fac*form(j)/sum
          savebm(j) = acoef(4971)*fac*form(j)/sum
          savibm(j) = (1.0-acoef(4971))*fac*form(j)/sum
!     CMOD
!     savefw(j) = 0.55*fac*form(j)/sum
!     savifw(j) = 0.45*fac*form(j)/sum
  911   continue
      endif
      if(from_xplasma_p .and. .not. use_tsc_ech) then
!     dum=(times-tr_ts1)/(tr_ts2-tr_ts1)
        dum=0.0
        do j=2,npsit
          savebm(j)=tr_peecrf(j,1)+dum*(tr_peecrf(j,2)-tr_peecrf(j,1))
!         use savibm for storing ecrf direct heating to ions
          savibm(j)=tr_piecrf(j,1)+dum*(tr_piecrf(j,2)-tr_piecrf(j,1))
        enddo
        bpowers_i = 0.0
        bpowers_e = 0.0
        do j=1,npsit-1
          bpowers_i=bpowers_i+savibm(j+1)*(vary(j+1)-vary(j))
          bpowers_e=bpowers_e+savebm(j+1)*(vary(j+1)-vary(j))
        enddo
        bpowers=bpowers_i+bpowers_e
        dum = bpowers
        if(acoef(4974) .gt. 0.0) then
          if(acoef(114) .ne. 0. .and. bpfeed .lt. bpowers)              &
     &             bpowers=bpfeed
          if(bpowers .gt. auxmax) bpowers = auxmax
        endif

        fac = usdp/usdt    
        if(dum .gt. 0.0) fac = fac*bpowers/dum
        do j=2, npsit
          savebm(j)=fac*savebm(j)
          savibm(j)=fac*savibm(j)
        enddo
        savebm(1)=savebm(2)
        savibm(1)=savibm(2)
      endif 
  912 continue
      ifirst_write=0
      return
 6600 continue
!
! --->  SPECIAL CODING FOR ACOEF(296)=5 TO GO WITH subroutine MISSIONC
!       DEFINES SOURCE TERMS NEEDED FOR OPTIMAL CONTROL
!
!
!.....initialize arrays to zero
      do 6599 j=1,npsi
      siun(j) = 0._R8
      seun(j) = 0._R8
!.....ICRH heating arrays
      savebm(j) = 0._R8
      savibm(j) = 0._R8
!.....Fast wave heating and current drive arrays
      savee(j) = 0._R8
      ajavcd(j) = 0._R8
!.....Lower-Hybrid heating and current drive arrays
      savelh(j) = 0._R8
      ajavlh(j) = 0._R8
      ajavlh2(j) = 0._R8
 6599 continue
      if(kcycle.le.0) go to 6699
!.............................................................................
!     Section 2.0 * * * Fast Wave Sources * * *
!.............................................................................
      bedge = gzero/(rmajor+rminor)
      bmag0 = gzero/xmag
!
!.....k-parallel
      akpar = 1.62E2_R8*bedge/sqrt(te(2))
      omega = 9.58E7_R8*bedge
      do 6610 j=2,npsit
!
!.....normalized square root of toroidal flux
      rho = sqrt((j-1.5_R8)/(npsit-1+fraclst))
!
!.....local magnetic field strength
      blocal = gzero / (rmajor + rho*rminor)
!
!.....deuterium density
      andu = 0.5_R8*anhy(j)
!
!.....absorption coefficient
      alphe = 4.44E-33_R8*bedge/blocal**3 * sqrt(te(2)*te(j))            &  
     &       * sqrt(andu+1.5_R8*andu)*ane(j)                             &  
     &       *exp(-(te(2)/te(j)))
      alphaea(j) = alphe
 6610 continue
      cnorm = 0._R8
      do 6620 j=2,npsit
      aint = 0._R8
      do 6630 l=j,npsit
      psiint = (l-1)*dpsi
      aint = aint + dpsi*(alphaea(l)/sqrt(psiint))
 6630 continue
      seun(j) = alphaea(j)*exp(-rminor/(2._R8*sqrt((npsit-1)*dpsi))*     &  
     & aint)
      cnorm = cnorm + seun(j)*vp(j)*dpsi
 6620 continue
      fnorm = pfwmc*1.E6_R8*usdp/usdt/cnorm
!
!.....source terms
      do 6611 j=2,npsit
      savee(j) = fnorm*seun(j)
      befo = 1.92E15_R8*rmajor*te(j)/ane(j)
      term1 = 8.42_R8*sqrt(te(j)/te(1))/(0.678_R8+zeffa(j))
      term2 = 4.13_R8/zeffa(j)**(0.707_R8)
      term3 = 8._R8/(5._R8+zeffa(j))*(te(1)/te(j))
      ajavfw(j)=befo*(term1+term2+term3)*(1._R8-0.5_R8*(ftrap(j)+        &  
     & ftrap(j-1)))                                                      &  
     &     *savee(j)*udsp/udst*usdi
 6611 continue
!
!.............................................................................
!     Section 3.0 * * * Lower Hybrid Sources * * *
!.............................................................................
!
!.....central deuterium density
      andu = 0.5_R8*anhy(2)
!
!.....broadcast frequency
      omegb = 2.79_R8*sqrt(0.5_R8*anhy(2))/ (1._R8+(1.03E-19_R8)*ane(2)/  &  
     & bmag0**2)
!.....k-parallel
      akpar = 9.31E-9_R8*sqrt(andu)*(1.0_R8+ 1._R8/sqrt(9.71E18_R8*      &  
     & bmag0**2/ane(2)                                                   &  
     &                                             + 1.0_R8))
!
      testop = 2.3E-5_R8*omegb/akpar
      if(testop .le. te(npsit)) go to 6699
!
      do 6612 j=2,npsit-1
      jj = 1+npsit-j
      if(te(jj).gt.testop) go to 6613
 6612 continue
      jj = 1
 6613 jstop = jj+1
      cnorm = 0._R8
      do 6614 j=jstop,npsit
      veth = 5.92E5_R8*sqrt(te(j))
      warg = omegb/(akpar*veth)
      seun(j) = ane(j)**2*sqrt(te(j))*exp(-warg**2)
      cnorm = cnorm + seun(j)*vp(j)*dpsi
 6614 continue
      fnorm = plhmc*1.E6_R8*usdp/usdt/cnorm
      do 6615 j=jstop,npsit
      savelh(j) = fnorm*seun(j)
      befo = 1.92E15_R8*rmajor*te(j)/ane(j)
      term1 = (1.26E6_R8/zeffa(j))*(akpar/omegb)*sqrt(te(j))
      term2 = 3.83_R8/zeffa(j)
      term3 = 2.28E-11_R8*(omegb/akpar)**2/(te(j)*(5._R8+zeffa(j)))
      ajavlh(j) = befo*(term1+term2+term3)*(1._R8-ftrap(j))*             &  
     &          savelh(j)*udsp/udst*usdi*xplas/gzero
 6615 continue
!
!.....define surface centered arrays
      do 720 j=1,npsit
      ajavlh2(j) = .5_R8*(ajavlh(j)+ajavlh(j+1))
  720 continue
 7630 format(" j       savelh       ajavlh")
 7631 format(i3,1p3e13.4)
 6699 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
