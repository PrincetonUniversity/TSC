      subroutine trtsci
      use trtsc
      use plasma_state_in_mod
      use clinam

      IMPLICIT NONE
      logical :: first_call= .true.
!     integer :: numargs
      logical :: ex
      character*30 :: temp_name
      integer :: get_shot_number
      INTEGER :: jsystem
      CHARACTER*256 :: cmd


      if(first_call) then
      call getcarg(1, suffix2, numargs)
      if(numargs .lt. 1) suffix2="I12345"
      call getcarg(2, temp_name, numargs)
      if(numargs .lt. 2) then
      transp_machine ="ITER"
      else
      transp_machine = trim(temp_name)
      endif
      call getcarg(3, temp_name, numargs)
      if(numargs .lt. 3) then
      transp_year ="09"
      else
      transp_year = trim(temp_name)
      endif
      call getcarg(4, temp_name, numargs)
!cj   June 8, 2010 add transp_number, mover shot_number from char to integer
!cj   for the convenience of importing trxpl profiles.
      if(numargs .lt. 4) then
!cj   shot_number ="10001"
      shot_number =10001
      else
      transp_number = trim(temp_name)
!cj   shot_number = trim(temp_name)
      shot_number = get_shot_number(transp_number)
      endif
      xplasma_filename = trim(suffix2)//"_ps.cdf"
      xplasma_filesave = trim(xplasma_filename)
      first_call = .false.
      tinit = times
      tfinal = acoef(29)
!
      if(acoef(4942) .gt. 0.0) use_user_pres=.true.
      if(acoef(4948) .gt. 0.0) use_user_chie=.true.
!
      if(acoef(4959) .gt. 0.0) use_swim=.true.
      if(acoef(4960) .gt. 0.0) use_transp=.true.
      if(acoef(4961) .gt. 0.0) use_plasma_state=.true.
      if(acoef(4962) .gt. 0.0) use_tsc_beam=.false.
      if(acoef(4963) .gt. 0.0) use_tsc_fp=.false.
      if(acoef(4964) .gt. 0.0) use_tsc_ich=.false.
      if(acoef(4965) .gt. 0.0) use_tsc_lhh=.false.
      if(acoef(4966) .gt. 0.0) use_tsc_ech=.false.
      if(.not.use_tsc_beam .and. .not.use_tsc_fp .and. .not.use_tsc_ich &
     &                   .and. .not.use_tsc_lhh .and. .not.use_tsc_ech) &
     &  use_tsc=.false.
      if(acoef(4967) .gt. 0.0) use_user_beam=.true.
      if(acoef(4968) .gt. 0.0) use_user_lhh=.true.
      if(acoef(4969) .gt. 0.0) use_user_fw=.true.
      if(acoef(4970) .gt. 0.0) use_user_ne=.true.
      if(acoef(4972) .gt. 0.0) use_user_rot=.true.
      if(acoef(4978) .gt. 0.0) use_user_radi=.true.
      if(acoef(4979) .gt. 0.0) use_user_ech=.true.
      if(use_transp) use_plasma_state=.true.
      if(use_swim) use_plasma_state=.true.

      if( acoef(4991) .gt. .0 .or.                                      &
          acoef(4992) .gt. .0 .or.                                      &
          acoef(4993) .gt. .0 .or.                                      &
          acoef(4994) .gt. .0      ) then
          call trxpl_init 

         if( acoef(4994) .gt. 0 ) then
            cmd= "cp"//                                                 &
                 " "//                                                  &
                 "xyz/init_test.mdescr"//                               &
                 " "//                                                  &
                 trim(transp_machine)//                                 &
                 ".mdescr"
            !cj ierr = jsystem(trim(cmd))
            call execute_command_line(trim(cmd),exitstat=ierr)
            if( ierr .ne.0) then
               write(6,*) "trtsci: copy nb machine mdescr err"
               ineg=70
            else
               write(6,*) "trtsci: copy nb machine mdescr ok"
            endif
         endif

         !cj nov 29 2010 added for acoef(4994) nubeam parameters
         call import_profile_from_transp 
      endif !cj nov-09-2010 moved from tsc, added nb machine mdescr

      if(use_plasma_state) then
      call init_plasma_state
      call getcarg(5, swim_init, numargs)
      if(numargs .ge. 5 .and. trim(swim_init) .eq. "init") then
      write(6,*) "init detected, TSC stops"
      call flush(6)
      acoef(29)=-1.0
      ncycle = kcycle
!     stop
      endif                                                              !init
      endif                                                              !ps
      endif                                                              !first call

      return
      end

      subroutine trtscm
      use trtsc
      use plasma_state_eq_mod
      use plasma_state_pa_mod
      use sawtooth_mod 
      use clinam
      IMPLICIT NONE
      logical :: tr_flag=.false.
      integer :: istat
      real*8 :: t
      logical :: ex
      integer :: i, jsystem
      character*1 :: nrec(10)
      character*256 :: extension
      CHARACTER*256 :: cmd, template, mod_file, trfile
      logical :: first_call = .true.

!

!     kill if there is an error in TSC
      if(ineg .ne. 0 .and. use_transp) then
        open(51,file="transp_save.dat",status="unknown",iostat=istat)
        write(51,*) "save transp for restart"
        close(51)

        open(51,file="transp_kill.dat",status="unknown",iostat=istat)
        write(51,*) "tsc error, job killed"
        close(51)
        return
      endif
      
      if(tr_flag) then
!
        inquire(file=trim(xplasma_filename),exist=ex)
        if(.not.ex) then
          extension=""
          xplasma_filename=trim(xplasma_filesave)//trim(extension) 
          inquire(file=trim(xplasma_filename),exist=ex)
!..THIS CODING CANNOT BE RIGHT....SHOULD BE A ELSE
           if(ex) then
             open(51,file=trim(xplasma_filename),status="unknown",iostat=istat)             
             close(51,status="delete")
           endif
        endif

        file_seq=file_seq+1
        istep = istep + 1
        ierr = 0
      
        call wrgeqdsk
        CALL put_plasma_state_eq
        call wrxpls
        CALL put_plasma_state_pa

        if(ierr .ne. 0) then
          ineg=70
          return
        endif
!     
        if(first_call) then
          if(.not. tr_init) then
            trfile=trim(suffix2)//"TR.DAT"
            write(6,*) "launching TRANSP background"
            cmd="fsp_runtr"//" "//trim(suffix2)//" "//trim(transp_machine)//        &
            " "//trim(transp_year)//" "//trim(trfile)//" "//"&"
            !cj ierr = jsystem(trim(cmd))
            call execute_command_line(trim(cmd),exitstat=ierr)
            if( ierr.ne.0) then
              write(6,*) "error launching TRANSP"
              stop
            endif
          else
            cmd="fsp_rstr"//" "//trim(suffix2)
            !cj ierr = jsystem(trim(cmd))
            call execute_command_line(trim(cmd),exitstat=ierr)
            if( ierr.ne.0) then
              write(6,*) "error launching TRANSP in restart"
              write(6,*) "try launching again from start"
              trfile=trim(suffix2)//"TR.DAT"
              write(6,*) "launching TRANSP background"
              cmd="fsp_runtr"//" "//trim(suffix2)//" "//trim(transp_machine)//        &
                 " "//trim(transp_year)//" "//trim(trfile)//" "//"&"
              !cj ierr = jsystem(trim(cmd))
              call execute_command_line(trim(cmd),exitstat=ierr)
              tr_init=.false.
                if( ierr.ne.0) then
                  write(6,*) "error launching TRANSP"
                  stop
                endif
            endif
          endif
          first_call=.false.
        endif      ! on first call


!       init
        if(.not.tr_init) then
          open(51,file="tmpfile",status="unknown",iostat=istat)
          t = acoef(ipos)
          write(51,"(f20.14,2x,A)") t, trim(xplasma_filename)
          close(51)
          extension="mv tmpfile transp_init.dat"
          !cj istat=jsystem(trim(extension))
          call execute_command_line(trim(extension),exitstat=ierr)
          extension=""
          tr_init=.true.
          write(6,*) "transp_init.dat created"
        else

!         step
!         save transp more often for debugging
          if( mod(istep,nstep_save) .eq. 0 ) then
            open(51,file="transp_save.dat",status="unknown",iostat=istat)
            write(51,*) "save transp for restart"
            close(51)
!           save tsc restart as well
            call rstrt2
          endif
 
!         iseq=iseq+1 ! ??????
          iseq=1
          trtimes(iseq)=acoef(ipos)
          xplasma_filelist(iseq)=xplasma_filename
          trevent(iseq)=" "
          open(51,file="tmpfile",status="unknown",iostat=istat)
          do i=iseq, iseq
            if(saw_tr) then
              extension="sawtooth_end"
              write(51,"(f20.14,2x,A,2x,A)") sawtooth_time, trim(saw1_filename), &
     &                              trim(extension)
              saw_tr=.false.
              extension=""
            endif
            write(51,"(f20.14,2x,A,2x,A)")      trtimes(i),                    &
     &                             trim(xplasma_filelist(i)),            &
     &                             trim(trevent(i))
          enddo
          close(51)
          extension="mv tmpfile transp_step.dat"
          !cj istat=jsystem(trim(extension))
          call execute_command_line(trim(extension),exitstat=ierr)
          extension=""
        endif    ! end of "else" on .not. tr_init

        inquire(file="transp_ready.dat",exist=ex)
        if(ex) then
          open(51,file="transp_ready.dat",status="old",iostat=istat)
          read(51,*) extension
          write(6,*) trim(extension)
          !cj istat=jsystem("mv transp_ready.dat transp_ready.old")
          extension="mv transp_ready.dat transp_ready.old"
          call execute_command_line(trim(extension),exitstat=ierr)
        endif

        call flush(6)
        extension="ps"
        call wait_transp(extension)
        acoef(ipos)=acoef(ipos)+acoef(ipos+1)   ! what is ipos?
        tr_flag=.false.

      endif    !cj end of if(tr_flag) then  dec-07-2009

!     time sync
      if(times+dt*udst .ge. acoef(ipos) .and. use_transp) then
        tr_flag=.true.
        dt=(acoef(ipos)-times)*usdt
        dts=acoef(ipos)-times
      endif

!
!     archiving TRANSP and TSC at the end of run 
!     if(times .ge. acoef(29) .and. from_xplasma) then
      if(times .ge. acoef(29) .and. kcycle .ge. ncycle) then

        if(use_transp) then
          open(51,file="transp_save.dat",status="unknown",iostat=istat)
          write(51,*) "save transp for restart"
          close(51)
        endif

        if(use_plasma_state) then
          call wrgeqdsk
          CALL put_plasma_state_eq
          call wrxpls
          CALL put_plasma_state_pa
!>>debug
          write(6,*) "after call to put_plasma_state_pa in trtscm", kcycle, times
          if(ierr .ne. 0) then
            ineg=70
          endif
          call rstrt2
          wrtrst = .true.
        endif

!     shut down TRANSP
        if(use_transp) then
          !cj istat=jsystem("date > transp_kill.dat")
          extension="date > transp_kill.dat"
          call execute_command_line(trim(extension),exitstat=ierr)
        endif
      endif

      return
      end
 
      subroutine getcarg(index, arg, numargs)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: index, numargs
      character*(*), intent(out) :: arg
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: numchars, ier
!-----------------------------------------------
      integer iargc, getarg
      numargs = iargc()
      numchars = getarg(index, arg)
#ifdef HAVE_MPI
      numargs = numargs - 4
#endif 
      end subroutine getcarg

      subroutine wait_transp(tr_state)

      use trtsc
      use plasma_state_mod
      use sawtooth_mod
      use clinam, only : ineg, nterm

      implicit none
      integer :: luntmp, luntty, idelete, istry, iswarn, isfail, istat
      character*160 :: tr_filename
      character*160 :: tr_locater, extension
      character*10  :: tr_state
      integer :: i, do_twice=1 !, jsystem

      luntmp = 51
      luntty = nterm                                          !nterm from tsc
      idelete = 0
      istry = 1                                               !try once per istry second
      iswarn = 120                                            !warn once every iswarn seconds
      isfail = 36000                                          !timeout after isfail seconds
      tr_locater="transp_ready.dat"

      if( istry .lt. 0 ) return

      do while (2>1)
      call wait_for_file(luntmp,tr_locater,luntty,idelete,istry,        &
     &                   iswarn,isfail,istat)

      if(istat .eq. 1) then
!
!     transp fails to return after isfail seconds
!     issue kill command to kill transp job
!     issue stop to kill tsc job
      close(luntmp)
      open(luntmp,file="transp_kill.dat",status="unknown",iostat=istat)
      write(luntty,*) "transp error, istat=1, job killed"
      write(luntmp,*) "transp error, job killed"
      close(luntmp)
      ineg=70
      call put_ineg(ineg)
      return
      endif
!
      close(luntmp)
      open(unit=luntmp,file=trim(tr_locater),status="unknown",          &
     &     position="rewind")
      read(luntmp,"(a160)",iostat=istat) tr_filename

      if(istat .ne. 0) then
      write(luntty,*) "error reading transp_ready file, istat= ",istat
      write(luntty,*) "transp_ready.dat is saved as _sav"
      extension="mv "//trim(tr_locater)//" "                            &
     &               //trim(tr_locater)//"_sav"
      !cj istat=jsystem(trim(extension))
      call execute_command_line(trim(extension),exitstat=ierr)
      extension=""
      tr_filename="ERROR"
      endif

      if( trim(tr_filename) .eq. "ERROR") then
      write(luntty,*) "transp error, transp_ready=ERROR, job killed"
      close(luntmp)
      open(luntmp,file="transp_kill.dat",status="unknown",iostat=istat)
      write(luntmp,*) "transp error, job killed"
      close(luntmp)
      ineg=70
      call put_ineg(ineg)
      return
      endif

      write(6,*) "transp_ready.dat: ", trim(tr_filename)
      call flush(6)
      close(luntmp)
      if(trim(tr_filename) .eq.                                          &
     &            trim(suffix2)//" (TRANSP) updated the state.") exit
      !cj istat=jsystem("mv transp_ready.dat transp_ready.old") 
      extension="mv transp_ready.dat transp_ready.old"
      call execute_command_line(trim(extension),exitstat=ierr)
      enddo

      if(trim(tr_state) .eq. "ps") then
      from_xplasma = .true.
      from_xplasma_p = .true.
      from_xplasma_j = .true.
      read_xplasma_p = .true.
      read_xplasma_j = .true.
      write(luntty,*) "xplasma file updated: ", trim(xplasma_filename)
      close(luntmp)
      !cj istat=jsystem("mv transp_ready.dat transp_ready.old") 
      extension="mv transp_ready.dat transp_ready.old"
      call execute_command_line(trim(extension),exitstat=ierr)

      call ps_get_plasma_state(ierr,trim(xplasma_filename)) 
      if(ierr .ne. 0) then
      write(6,*) "attempt to read plasma state failed",ierr
      ineg=70
      call put_ineg(ineg)
      !cj istat=jsystem("date > transp_kill.dat")
      extension="date > transp_kill.dat"
      call execute_command_line(trim(extension),exitstat=ierr)
      endif
      endif

      if(trim(tr_state) .eq. "saw1") then
      from_xplasma = .true.
      from_xplasma_p = .true.
      from_xplasma_j = .true.
      read_xplasma_p = .true.
      read_xplasma_j = .true.
      write(luntty,*) "saw1 file updated: ", trim(saw1_filename)
      close(luntmp)
      !cj istat=jsystem("mv transp_ready.dat transp_ready.old") 
      extension="mv transp_ready.dat transp_ready.old"
      call execute_command_line(trim(extension),exitstat=ierr)

      call ps_get_plasma_state(ierr,trim(saw1_filename),state=saw1)
      if(ierr .ne. 0) then
      write(6,*) "attempt to read saw1 state failed",ierr
      ineg=70
      call put_ineg(ineg)
      !cj istat=jsystem("date > transp_kill.dat")
      extension="date > transp_kill.dat"
      call execute_command_line(trim(extension),exitstat=ierr)
      endif
      endif


      return
      end

      subroutine put_ineg(neg)
      use clinam, ONLY : ineg
      integer :: neg

      ineg=neg

      return
      end

