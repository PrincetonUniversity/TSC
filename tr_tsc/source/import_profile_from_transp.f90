      MODULE Trxpl_Ps 
      USE plasma_state_mod
      SAVE gog
      type (plasma_state) :: gog 

!...  passed from tsc
!     INTEGER npsit
!     real*8 times
!     character*80 transp_number,transp_machine
!     real*8 psimin, psilim
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: xsv2, xsv

!...  local to trxpl_ps
      real*8 fetched_t0, fetched_t1
      integer :: num_timeslice
      REAL*8, ALLOCATABLE, DIMENSION(:) :: timeslice
      character*40, ALLOCATABLE, DIMENSION(:) :: psfilename

      INTEGER nrho,new_nrho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rho, ns, ts
      REAL*8, ALLOCATABLE, DIMENSION(:) :: new_rho, new_ns, new_ts
      REAL*8, ALLOCATABLE, DIMENSION(:) :: new_nsi, new_tsi
      REAL*8, ALLOCATABLE, DIMENSION(:) ::          new_ns_t1, new_ts_t1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: new_nsi_t1, new_tsi_t1
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: rho_zc, new_rho_zc
      logical :: smooth = .true. , norm=.false.  
      logical :: first_fetch= .true.  

      END MODULE Trxpl_Ps 


      subroutine import_profile_from_transp

      USE CLINAM, ONLY: times,                                          &
     &                  npsit, psilim, psimin, xsv2,                    &
     &                  adp ,ade ,adn, vp, vpg,                         &
     &                  usdd, usdh,                                     &
     &                  acoef, ineg, kcycle
      USE trtsc, ONLY: transp_machine, transp_year, transp_number
      USE Trxpl_Ps

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100) 

      integer i, ierr, istat
!cj      real*8 tmp_adp, tmp_ade, tmp_adn


      if( (times .le. timeslice(1)            )   .or.                  &
          (times .ge. timeslice(num_timeslice)) ) then
#ifdef CJDebug
         write(*,*)                                                     &
         'trxpl: err tsc times too small',                              &
         '       no transp profiles exist for this time',               &
         '       check cards 14 and 18',                                &
         '       your tsc times=',times,                                &
         '       smallest tansp times=',timeslice(1)
#endif
!        call trxpl_kill
!        stop 'trxpl Error: tsc times too small, check cards 14 and 18' 
!        return
            first_fetch=.true.
      else
#ifdef CJDebug
         write( *,'(8x, a40, 1pe13.6)' )                                &
         "trxpl: extract profile at tsc times=", times
#endif

!.... start with an impossible times
!.... if plasma state has already been fetched in "times"
!.... we skip fetch and jump to rezoning
         if(first_fetch) then 
            fetched_t0=-10._R8
            fetched_t1=-10._R8

!.... if first_fetch set true  here, trxpl fetch ps every "times"
!.... if first_fetch set false here, trxpl fetch only transp timeslice. Then
!....                                we interpolate profiles at "times" below
            first_fetch=.false.
!2010july28 first_fetch=.true.
         endif 
      endif 

!     times=.25_R8
!     npsit=80 
!     psimin= -.19457_R8
!     psilim= -.15104_R8
!     write(*,*) "import_profile_from_transp: "                         &
!               ,"transp_machine = ", trim(transp_machine)              &
!               ,"transp_number = ", trim(transp_number)                &
!               ,"timeslice(1) = ", timeslice(1)
!               ,"npsit = ", npsit                                      &
!               ,"psimin = ", psimin                                    &
!               ,"psilim = ", psilim

!     call trxpl_init
      call trxpl_fetch
!     call trxpl_kill

!     IF(.not.ALLOCATED(xsv2)) ALLOCATE( xsv2(npsit), STAT=istat)
!        if (istat .ne. 0) stop 'Allocation Error : xsv2' 
!     IF(.not.ALLOCATED(xsv)) ALLOCATE( xsv(npsit), STAT=istat)
!        if (istat .ne. 0) stop 'Allocation Error : xsv' 
!     do i = 1, npsit
!     xsv2(i) = psimin + sqrt(float(i-1)/float(npsit-1)) * (psilim - psimin)
!     enddo
!     xsv(2:npsit) = ( xsv2(1:npsit-1) + xsv2(2:npsit) )/2._R8

      new_nrho = npsit
      ALLOCATE( new_rho(new_nrho), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_rho' 
      ALLOCATE( new_ns(new_nrho-1), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_ns' 
      ALLOCATE( new_ts(new_nrho-1), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_ts' 
      ALLOCATE( new_nsi(new_nrho-1), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_nsi' 
      ALLOCATE( new_tsi(new_nrho-1), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_tsi' 
      ALLOCATE( new_ns_t1(new_nrho-1), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_ns_t1' 
      ALLOCATE( new_ts_t1(new_nrho-1), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_ts_t1' 
      ALLOCATE( new_nsi_t1(new_nrho-1), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_nsi_t1' 
      ALLOCATE( new_tsi_t1(new_nrho-1), STAT=istat)
         if (istat .ne. 0) stop 'Allocation Error : new_tsi_t1' 

!.... rezone profile from transp grid to tsc grid
      call trxpl_rezone

!.... linear interpolation for profiles at times <- [t0,t1]
!.... xx = xx_t0 + (t - t0)/(t1 - t0) * ( xx_1 - xx_0)
!cj      tmp_adp=0.
!cj      tmp_ade=0.
!cj      tmp_adn=0.
      if(first_fetch) then 
       do i=2,npsit

         if( acoef(4993) .gt. .0 )                                      &
         adp(i) = ( new_ns(i-1) * usdd *                                &
                    new_ts(i-1) * 1000._R8 * usdh/usdd +                &
                    new_nsi(i-1)* usdd *                                &
                    new_tsi(i-1)* 1000._R8 * usdh/usdd    ) * vpg(i)
!cj         tmp_adp=tmp_adp + adp(i)

         if( acoef(4992) .gt. .0 )                                      &
         ade(i) = ( new_ns(i-1) * usdd *                                &
                    new_ts(i-1) * 1000._R8 * usdh/usdd    ) * vpg(i)
!cj         tmp_ade=tmp_ade + ade(i)

         if( acoef(4991) .gt. .0 )                                      &
         adn(i) = ( new_ns(i-1) *usdd                                   &
!cj                +new_nsi(i-1)*usdd                                   &
                  ) * vp(i)
!cj         tmp_adn=tmp_adn + adn(i)
!     write(93,*) "    i  ne          te          ne*usdd     te*usdh"
!     write(93,9661) i, new_ns(i-1),      new_ts(i-1)                   &
!                     , new_ns(i-1)*usdd, new_ts(i-1)*1000._R8*usdh/usdd
!9661 format(i5, 1p4e12.4) 
       enddo
!cj      write(*,9661) kcycle, tmp_adp, tmp_ade, tmp_adn
!cj9661 format("cjDebug kcycle=", i8, 1p3e12.4) 

      else
       do i=2,npsit
         if( acoef(4993) .gt. .0 )                                      &
         adp(i) = ( new_ns(i-1) * usdd *                                &
                    new_ts(i-1) * 1000._R8 * usdh/usdd +                &
                    new_nsi(i-1)* usdd *                                &
                    new_tsi(i-1)* 1000._R8 * usdh/usdd                  &
                  ) * vpg(i)                                            &
                 +                                                      &
                    (times - fetched_t0) / (fetched_t1 - fetched_t0) *  &
                  ( (new_ns_t1(i-1) - new_ns(i-1))*usdd *               &
                    (new_ts_t1(i-1) - new_ts(i-1))*1000._R8*usdh/usdd + &
                    (new_nsi_t1(i-1) - new_nsi(i-1))* usdd *            &
                    (new_tsi_t1(i-1) - new_tsi(i-1))*1000._R8*usdh/usdd &
                  ) * vpg(i)

         if( acoef(4992) .gt. .0 )                                      &
         ade(i) = ( new_ns(i-1) * usdd *                                &
                    new_ts(i-1) * 1000._R8 * usdh/usdd                  &
                  ) * vpg(i)                                            &
                 +                                                      &
                    (times - fetched_t0) / (fetched_t1 - fetched_t0) *  &
                  ( (new_ns_t1(i-1) - new_ns(i-1))*usdd *               &
                    (new_ts_t1(i-1) - new_ts(i-1))*1000._R8*usdh/usdd   &
                  ) * vpg(i)

         if( acoef(4991) .gt. .0 )                                      &
         adn(i) = ( new_ns(i-1) *usdd                                   &
!cj                +new_nsi(i-1)*usdd                                   &
                  ) * vp(i)                                             &
                 +                                                      &
                    (times - fetched_t0) / (fetched_t1 - fetched_t0) *  &
                  ( (new_ns_t1(i-1) - new_ns(i-1)) *usdd                &
!cj                +(new_nsi_t1(i-1) - new_nsi(i-1))*usdd               &
                  ) * vp(i) 
       enddo
      endif

!       call trxpl_kill
!     stop
!.... needed or allocate/deallocate everytime import_profile_from_transp
!.... is called since npsit changes at every kcycle
      DEALLOCATE (new_rho, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_rho'
      DEALLOCATE (new_ns, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_ns'
      DEALLOCATE (new_ts, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_ts'
      DEALLOCATE (new_nsi, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_nsi'
      DEALLOCATE (new_tsi, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_tsi'
      DEALLOCATE (new_ns_t1, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_ns_t1'
      DEALLOCATE (new_ts_t1, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_ts_t1'
      DEALLOCATE (new_nsi_t1, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_nsi_t1'
      DEALLOCATE (new_tsi_t1, STAT=istat) 
         if (istat .ne. 0) stop 'Deallocation Error : new_tsi_t1'

      return
      end


      subroutine trxpl_init

      USE CLINAM, ONLY: ineg
      USE trtsc, ONLY: transp_machine, transp_year, transp_number
      USE Trxpl_Ps

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

      character*4 :: mode_str
      character*250 dpath, opath
      character*200 cmd, mdescr
      integer i, ierr, istat !, jsystem
      integer :: iprotimes, ios
      character*120 :: protimes
      character*8 :: oprfix

      logical :: mhdeq_flag   ! -mhdeq indicator

      integer :: iBccw        ! B field direction option
      integer :: iJccw        ! current direction option


!...  init
!cj   transp_machine="NSTX"
!cj   transp_number="124379A10"
!...  jan 04, 2011
!     trxpl_server crashes after a certain amount of time fetch
!     Now we get all the data in plasma state before hand
!     and just do interpolation now 
#ifdef NO_TRXPL
      write(*,*) "no trxpl: read data from /p/tsc/transp_data/ start"
#else
#ifdef OLD
!... trxpl_server
      write(dpath,*)      "transpgrid.pppl.gov:"                        &
                          ,"transp_",trim(transp_machine)               &
                          ,"\(",trim(transp_machine)                    &
                          ,",",trim(transp_number),"\)"
      cmd = "trxpl_server"//                                            &
            " "//                                                       &
            "INIT"//                                                    &
            " "//                                                       &
            "MDS+"//                                                    &
            " "//                                                       &
            trim(dpath)//                                               &
            " "//                                                       &
            "xyz transp.mdescr >& log.trxpl_init"
      istat = jsystem(trim(cmd))
            if( istat.ne.0) then
              write(*,*) "trxpl_init: err trxpl_server init"
              ineg=62
              return
            else
              call mystop('trxpl_init: ok trxpl_server init')
            endif
#else
!... trxpl2ps
      write(dpath,*)      "transpgrid.pppl.gov:"                        &
                          ,"transp_",trim(transp_machine)               &
                          ,"(",trim(transp_machine)                     &
                          ,",",trim(transp_number),")"
      opath = "xyz"
      mode_str = "MDS+"
      mhdeq_flag = .FALSE.
      iBccw = 0
      iJccw = 0
      oprfix = "trps_"

      ! select a reference machine description 
      mdescr = trim(transp_machine)//                                   &
               "_05.mdescr"
      call trxplib_clear_mdescr
      call trxplib_add_mdescr(trim(mdescr),ierr)

      ! initialize Plasma State object (gog) and module 
      call ps_init_user_state(gog, "trxpl2ps", ierr)
      if(ierr.ne.0) then 
         write(*,*) "trxpl2ps_init: err ps_init_user_state"
         ineg=62
         return
         ! call exit(10)
      endif
  
      ! initialize trxplib library options 
      call trxplib_ps_reset
      call trxplib_mhdeq_opt(mhdeq_flag)
      call trxplib_ccw_options(iBccw,iJccw)
  
      ! select TRANSP run; error checks 
      call trxplib_ps_open(mode_str, dpath, opath, trim(oprfix), ierr)
      if(ierr.ne.0) then 
         write(*,*) "trxpl2ps_init: err trxplib_ps_open"
         ineg=62
         return
         ! call exit(11)
      else
              call mystop('trxpl_init: ok trxpl2ps init')
      endif 
#endif
#endif

!.... read in the transp timeslice 
      iprotimes=11
#ifdef NO_TRXPL
      protimes= "/p/tsc/transp_data/"//                                 &
                trim(transp_machine)//                                  &
                "_"//                                                   &
                trim(transp_number)//                                   &
                "/abc/protimes.dat"
#else
      protimes= "xyz/protimes.dat"
#endif
      write(*,*) " --- protimes =", protimes
      open(iprotimes,file=trim(protimes),status='old',iostat=ios)
      if(ios .ne. 0) then
         write(*,*) "trxpl_init: protimes.dat open err"
         ineg=62
         return
      else
         call mystop('trxpl_init: Extract TRANP timeslices ...')

         rewind (unit = iprotimes) 
         read(iprotimes,'(1x,i5)') num_timeslice
         write( *,'(8x, a20, i5, 2x, a15)' )                            &
         "protimes.dat has", num_timeslice, "time slices"
         IF(.not.ALLOCATED(timeslice)) ALLOCATE( timeslice(num_timeslice), STAT=istat)
            if (istat .ne. 0) stop 'Allocation Error : timeslice' 
#ifdef NO_TRXPL
         IF(.not.ALLOCATED(psfilename)) ALLOCATE( psfilename(num_timeslice), STAT=istat)
            if (istat .ne. 0) stop 'Allocation Error : psfilename' 
#endif
         do i = 1, num_timeslice
#ifdef NO_TRXPL
            read(iprotimes,'(1x,1pe13.6, A40)') timeslice(i), psfilename(i)
!           write(*,*) "timeslice(",i,")=", timeslice(i), trim(psfilename(i))
#else
            read(iprotimes,'(1x,1pe13.6)') timeslice(i)
!           write(*,*) "timeslice(",i,")=", timeslice(i)
#endif
         enddo
      endif 
      close(iprotimes)

      return
      end


      subroutine trxpl_kill

      USE CLINAM, ONLY: ineg
      USE Trxpl_Ps

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

      character*100 tmpchar 
      character*200 cmd
      integer ierr, istat !, jsystem


      DEALLOCATE (timeslice, STAT=istat) 
         if (istat .ne. 0) stop 'trxpl_kill: Deallocate Error timeslice'

#ifdef NO_TRXPL
!...  jan 04, 2011
!     trxpl_server crashes after a certain amount of time fetch
!     Now we get all the data in plasma state before hand
!     and just do interpolation now

      DEALLOCATE (psfilename, STAT=istat) 
         if (istat .ne. 0) stop 'trxpl_kill:Deallocate Error psfilename'
      write(*,*) "no trxpl: read data from /p/tsc/transp_data/ end"
      return
#endif

!...  kill
#ifdef OLD
!... trxpl_server
      cmd = "trxpl_server"//                                            &
            " "//                                                       &
            "KILL"//                                                    &
            " "//                                                       &
            "xyz >& log.trxpl_kill"
      istat = jsystem(trim(cmd))
            if( istat.ne.0) then
              write(6,*) "trxpl_kill: err trxpl_server kill"
              ineg=62
              return
            else
              call mystop('trxpl_kill: ok trxpl_server kill')
            endif
#else
!... trxpl2ps
      call ps_free_user_state(gog, ierr)
      if(ierr.ne.0) then
         write(6,*) "trxpl2ps_kill: err trxpl_server kill"
         ineg=62
         return
      else
          call mystop('trxplps_kill: ok trxpl_server kill')
      endif

      call trxplib_ps_close 
#endif

      return
      end


      subroutine trxpl_fetch

      USE CLINAM, ONLY: times, ineg
      USE trtsc, ONLY: transp_machine, transp_year, transp_number
      USE Trxpl_Ps

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

      character*100 tmpchar 
      character*200 cmd
      character*80 ps_filename
      integer i, ierr, istat, i0, i1 !, jsystem

      real*8 t0, t1, tmpreal


!.... and compare it to times, fetched_t0, fetched_t1
!.... then decide whether "fetch" a new plasma state 
!.... or use the one just fetched at last kcycle
      t0 = timeslice(1)
      t1 = timeslice(num_timeslice)
      do i = 1, num_timeslice
      if(timeslice(i) .le. times) then
         t0 = timeslice(i)
#ifdef NO_TRXPL
         i0=i  ! the first ps for interpolation
#endif
      endif
      if(timeslice(i) .ge. times) then 
         t1 = timeslice(i)
#ifdef NO_TRXPL
         i1=i  ! the second ps for interpolation
#endif
         exit
      endif
      enddo

      if(t0.gt.fetched_t0 .or. t1.gt.fetched_t1) then
         fetched_t0=t0
         fetched_t1=t1
         if(first_fetch) then 
            write( *,'(8x, a45, 1e13.6)' )                              &
            "trxpl_fetch: extract ps at times =", times
         else
            write( *,'(8x, a42, 2e13.6)' )                              &
            "trxpl_fetch: extract ps between=", t0, t1 
         endif
      else 
!        write( *,'(8x, a42, 2e13.6)' )                                 &
!        "trxpl_fetch: ps already extracted between=", t0, t1
         return
      endif

!...  fetch t0
      if(first_fetch) then 
      ! extract interpolated data
         if(times .le. timeslice(1)) then
         ! if tsc time is too small, use the first transp data slice
            write(tmpchar,*) timeslice(1)
            tmpreal = timeslice(1)   ! used in trxpl2ps
#ifdef NO_TRXPL
         i0=1  ! the first ps
#endif
         else if(times .ge. timeslice(num_timeslice)) then
         ! if tsc time is too big, use the last transp data slice
            write(tmpchar,*) timeslice(num_timeslice)
            tmpreal = timeslice(num_timeslice)   ! used in trxpl2ps
#ifdef NO_TRXPL
         i1=num_timeslice  ! the last ps
#endif
         else
            write(tmpchar,*) times
            tmpreal = times   ! used in trxpl2ps
         endif
      else
      ! extract first data to interpolat
         write(tmpchar,*) t0
         tmpreal = t0   ! used in trxpl2ps
#ifdef NO_TRXPL
      ! already get i0 and i1 above
#endif
      endif

#ifdef NO_TRXPL
      cmd = "cp"//                                                      &
            " "//                                                       &
            "/p/tsc/transp_data/"//                                     &
            trim(transp_machine)//                                      &
            "_"//                                                       &
            trim(transp_number)//                                       &
            "/abc/"//                                                   &
            trim(ADJUSTL(psfilename(i0)))//                             &
            " "//                                                       &
            "gog_ps.cdf"
      istat = jsystem(trim(cmd))
            if( istat.ne.0) then
              write(6,*) "cp: err cp ps0 gog_ps.cdf"
              ineg=62
              return
            endif
#else
#ifdef OLD
!... trxpl_server
      cmd = "trxpl_server"//                                            &
            " "//                                                       &
            "FETCH"//                                                   &
            " "//                                                       &
            trim(tmpchar)//                                             &
            " "//                                                       &
            "xyz"//                                                     &
            " "//                                                       &
            "gog_ps >& log.trxpl_fetch"
      istat = jsystem(trim(cmd))
            if( istat.ne.0) then
              write(6,*) "trxpl_fetch: err trxpl_server fetch t0"
              ineg=62
              return
            endif
#else
!... trxpl2ps
     call trxplib_ps_get(tmpreal, gog, ierr)
     if(ierr.ne.0) then
         write(6,*) "trxpl2ps_fetch: err trxplib_ps_get t0"
         ineg=62
         return
     endif

     call ps_store_plasma_state(ierr,state=gog,filename='gog_ps.cdf')
     if(ierr.ne.0) then
        write(6,*) "trxpl2ps_fetch: err ps_store_plasma_state t0"
        ineg=62
        return
        ! call exit(12)
     endif 
#endif
#endif

      if(first_fetch) return

!...  fetch t1
      ! extract second data to interpolat
      write(tmpchar,*) t1 
      tmpreal = t0   ! used in trxpl2ps
#ifdef NO_TRXPL
      cmd = "cp"//                                                      &
            " "//                                                       &
            "/p/tsc/transp_data/"//                                     &
            trim(transp_machine)//                                      &
            "_"//                                                       &
            trim(transp_number)//                                       &
            "/abc/"//                                                   &
            trim(ADJUSTL(psfilename(i1)))//                             &
            " "//                                                       &
            "gog_ps_t1.cdf"
      istat = jsystem(trim(cmd))
            if( istat.ne.0) then
              write(6,*) "cp: err cp ps0 gog_ps_t1.cdf"
              ineg=62
              return
            endif
#else
#ifdef OLD
!... trxpl_server
      cmd = "trxpl_server"//                                            &
            " "//                                                       &
            "FETCH"//                                                   &
            " "//                                                       &
            trim(tmpchar)//                                             &
            " "//                                                       &
            "xyz"//                                                      &
            " "//                                                       &
            "gog_ps_t1 >& log.trxpl_fetch_t1"
      istat = jsystem(trim(cmd))
            if( istat.ne.0) then
              write(6,*) "trxpl_fetch: err trxpl_server fetch t1"
              ineg=62
              return
            endif
#else
!... trxpl2ps
     call trxplib_ps_get(tmpreal, gog, ierr)
     if(ierr.ne.0) then
        write(6,*) "trxpl2ps_fetch: err trxplib_ps_get t1"
        ineg=62
        return
     endif

     call ps_store_plasma_state(ierr,state=gog,filename='gog_ps_t1.cdf')
     if(ierr.ne.0) then
        write(6,*) "trxpl2ps_fetch: err ps_store_plasma_state t1"
        ineg=62
        return
        ! call exit(12)
     endif 
#endif 
#endif 

      return
      end


      subroutine trxpl_rezone

      USE CLINAM, ONLY: times,                                          &
     &                  npsit, psilim, psimin, xsv2, ineg
      USE Trxpl_Ps

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

      character*100 tmpchar 
      character*80 ps_filename
      integer i, ierr, istat


#ifdef OLD
!... trxpl_server
      ps_filename="xyz/gog_ps.cdf"
#else
!... trxpl2ps
      ps_filename="gog_ps.cdf"
#endif
      call ps_get_plasma_state(ierr,trim(ps_filename),state=gog)
            if(ierr .ne. 0) then 
               write(*,*) "trxpl_rezone: err ps_get_plasma_state gog"
               ineg=62
               return
#ifdef CJDebug
            else
               write(*,*) "trxpl_rezone: ok ps_get_plasma_state gog"
#endif
            endif

!     write(*,*) "gog%t0=", gog%t0, "gog%t1=", gog%t1 
!     write(*,*) "gog%Ts(nrho=",gog%nrho-1,"nspec_th=",gog%nspec_th,")" 
!     do i=0, gog%nspec_th
!        write(*,*) "gog%id_ns(",i, ")=", gog%id_ns(i),                   &
!                   "gog%id_Ts(",i, ")=", gog%id_Ts(i)
!     enddo

      nrho=gog%nrho
      IF(.not.ALLOCATED(rho)) ALLOCATE( rho(nrho), STAT=istat) 
        if (istat .ne. 0) stop 'trxpl_rezone: Allocation err rho' 
      IF(.not.ALLOCATED(ns)) ALLOCATE( ns(nrho-1), STAT=istat) 
         if (istat .ne. 0) stop 'trxpl_rezone: Allocation err ns' 
      IF(.not.ALLOCATED(ts)) ALLOCATE( ts(nrho-1), STAT=istat) 
         if (istat .ne. 0) stop 'trxpl_rezone: Allocation err ts'
!     IF(.not.ALLOCATED(rho_zc)) ALLOCATE( rho_zc(nrho-1), STAT=istat)
!        if (istat .ne. 0) stop 'Allocation Error : rho_zc' 

      rho(:)=gog%rho(:)
      ns(:)=gog%ns(:,0)
      ts(:)=gog%Ts(:,0)
!     rho_zc = ( gog%rho(1:gog%nrho-1) + gog%rho(2:ps%nrho) )/2 
!cj   do i=1, nrho-1
!cj   write(*,*) "orig ns(",i,")=",rho_zc(i),ns(i),ts(i)
!cj   enddo
!cj   write(*,*)

!cj  SUBROUTINE ps_rho_rezone1(rho, id, ansv, ierr, &
!cj       state, curdens, nonorm, zonesmoo, zonediff) 
!     IF(.not.ALLOCATED(new_rho_zc)) ALLOCATE( new_rho_zc(new_nrho-1), STAT=istat)
!        if (istat .ne. 0) stop 'Allocation Error : new_rho_zc' 

      do i = 1, new_nrho
         !cj poloidal flux
!cj      new_rho(i) = (xsv2(i) - psimin) / (psilim - psimin)
!cj      write(*,*) "new_rho", i, new_rho(i), xsv2(i)

         !cj square root of toroidal flux, compatible with trxpl output
         new_rho(i) = sqrt( (real(i)-1.)/(real(new_nrho)-1.) ) 
      enddo
      new_rho(1) = 0._R8
      new_rho(new_nrho) = 1._R8
!     new_rho_zc = ( new_rho(1:new_nrho-1) + new_rho(2:new_nrho) )/2._R8

      call ps_rho_rezone1(new_rho,gog%id_ns(0),new_ns,ierr,             &
                          state=gog, zonesmoo=smooth, nonorm=norm)
         if(ierr .ne. 0) then 
           write(*,*) "trxpl_rezone: err ps_rho_rezone1 new_ns"
           ineg=62
           return
#ifdef CJDebug
         else
           write(*,*) "trxpl_rezone: ok ps_rho_rezone1 new_ns"
#endif
         endif

      call ps_rho_rezone1(new_rho,gog%id_Ts(0),new_ts,ierr,             &
                          state=gog, zonesmoo=smooth, nonorm=norm)
         if(ierr .ne. 0) then 
           write(*,*) "trxpl_rezone: err ps_rho_rezone1 new_ts"
           ineg=62
           return
#ifdef CJDebug
         else
           write(*,*) "trxpl_rezone: ok ps_rho_rezone1 new_ts"
#endif
         endif 

      call ps_ti_rezone(new_rho,ierr,                                   &
              state=gog,zonesmoo=smooth,ni=new_nsi,ti=new_tsi);
         if(ierr .ne. 0) then 
           write(*,*) "trxpl_rezone: err ps_rho_rezone1 new_nsi new_tsi"
           ineg=62
           return
#ifdef CJDebug
         else
           write(*,*) "trxpl_rezone: ok ps_rho_rezone1 new_nsi new_tsi"
#endif
         endif 

!cj   do i=1, new_nrho-1
!cj      write(*,*) "rezoned ns(",i,")=", new_rho_zc(i),new_ns(i),new_ts(i), new_nsi(i),new_tsi(i)
!cj   enddo
!cj   write(*,*)

      if(first_fetch) return

             !!!!!!!! .......... t1 .......... !!!!!!!!
#ifdef OLD
!... trxpl_server
      ps_filename="xyz/gog_ps_t1.cdf"
#else
!... trxpl2ps
      ps_filename="gog_ps_t1.cdf"
#endif
      call ps_get_plasma_state(ierr,trim(ps_filename),state=gog)
            if(ierr .ne. 0) then 
               write(*,*) "trxpl_rezone: err ps_get_plasma_state gog t1"
               ineg=62
               return
            endif

      if (nrho .ne. gog%nrho) then
         write(*,*) "trxpl_rezone: err ps t1 dimension nrho wrong"
         ineg=62
         return
      else
         rho(:)=gog%rho(:)
         ns(:)=gog%ns(:,0)
         ts(:)=gog%Ts(:,0)
      endif

      call ps_rho_rezone1(new_rho,gog%id_ns(0),new_ns_t1,ierr,          &
              state=gog, zonesmoo=smooth, nonorm=norm)
         if(ierr .ne. 0) then 
           write(*,*) "trxpl_rezone: err ps_rho_rezone1 new_ns"
           ineg=62
           return
         endif

      call ps_rho_rezone1(new_rho,gog%id_Ts(0),new_ts_t1,ierr,          &
              state=gog, zonesmoo=smooth, nonorm=norm)
         if(ierr .ne. 0) then 
           write(*,*) "trxpl_rezone: err ps_rho_rezone1 new_ts"
           ineg=62
           return
         endif 

      call ps_ti_rezone(new_rho,ierr,                                   &
              state=gog,zonesmoo=smooth,ni=new_nsi_t1,ti=new_tsi_t1);
         if(ierr .ne. 0) then 
           write(*,*) "trxpl_rezone: err ps_rho_rezone1 new_nsi new_tsi"
           ineg=62
           return
         endif 

      return
      end


      SUBROUTINE mystop ( ErrMsg )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*(*) ErrMsg
      write( *,'(8x, a40)' ) ErrMsg
!cj   write(*,'('''',a40 )') ErrMsg
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
