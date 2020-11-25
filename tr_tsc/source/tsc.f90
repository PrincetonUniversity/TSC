      MODULE Input_Mod
      USE kind_spec_mod

      IMPLICIT NONE
      INTEGER :: nt_mod
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ncard_mod
      REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: time_mod
      CHARACTER*8, ALLOCATABLE, DIMENSION(:,:) :: card_mod*80
      LOGICAL*1, ALLOCATABLE, DIMENSION(:) :: iread_mod

      CONTAINS

      SUBROUTINE read_input_mod

      IMPLICIT NONE
      LOGICAL :: ex
      INTEGER :: istat, i, ic
      INTEGER :: jsystem, ierr
      CHARACTER*256 :: cmd

      nt_mod = 0
      INQUIRE(file="input_mod.dat",exist=ex)
      if(.NOT.ex) return

      CLOSE(51)
      OPEN(51,file="input_mod.dat",form="formatted",                    &
     &                             status="old",iostat=istat)
      if(istat .ne. 0) then
      WRITE(6,*) "input_mod.dat open error"
      stop
      endif
      READ(51,*) nt_mod
      ALLOCATE(ncard_mod(nt_mod))
      ALLOCATE(iread_mod(nt_mod))
      ALLOCATE(time_mod(nt_mod))
      ALLOCATE(card_mod(50,nt_mod))

      iread_mod = .false.

      DO i=1, nt_mod
      READ(51,*) time_mod(i)
      ic = 0
      DO WHILE (2>1)
      ic = ic + 1
      if(ic .gt. 50) then
      write(6,*) "number of entries in input_mod > limit"
      !cj istat=jsystem("date > transp_kill.dat")
      cmd= "date > transp_kill.dat"
      call execute_command_line(trim(cmd),exitstat=ierr)
      stop
      endif
      READ(51,("(A80)")) card_mod(ic,i)(1:80)
      if(card_mod(ic,i)(1:2) .eq. "99") exit
      ENDDO
      ncard_mod(i)=ic
      ENDDO
      CLOSE(51)

      WRITE(6,*) "*** INPUT MODIFICATION CARDS READ ***"
      DO i=1, nt_mod
      WRITE(6,*) time_mod(i)
      DO ic=1, ncard_mod(i)
      WRITE(6,("(A80)"))  card_mod(ic,i)(1:80)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE read_input_mod


      SUBROUTINE param_time_mod(nin,time)
      INTEGER :: nin, i, j, k
      REAL*8 :: time

      IF(nt_mod .le. 0) RETURN
      IF(time .lt. time_mod(1)) RETURN

      DO i=1, nt_mod
      j=nt_mod-i+1
      if(time .ge. time_mod(j)) EXIT
      ENDDO

      IF(.NOT.iread_mod(j)) THEN
      iread_mod(j)=.true.
      OPEN(51,file="input_mod.tmp",status="unknown")
      DO k=1, ncard_mod(j)
      write(51,("(A80)")) card_mod(k,j)(1:80)
!     write(6,("(A80)")) card_mod(k,j)(1:80)
      ENDDO
      CLOSE(51)
      CALL flush(6)

      CLOSE(nin)
      OPEN(nin,file="input_mod.tmp",status="unknown")
      CALL inpt
      CLOSE(nin,status="delete")
      ENDIF

      RETURN
      END SUBROUTINE param_time_mod


      END MODULE input_mod

!
!
!.....princeton tokamak simulation code version unx10.4 : 01 Apr 00
!
!
!
!                      * * * contributers * * *
!
!                      c. bathke
!                      j. delucia
!                      s.c.jardin
!                      c. kessel
!                      k.m.ling
!                      f.marcus
!                      b.merrill
!                      r.pillsbury
!                      n.pomphrey
!                      m.reusch
!                      r.sayer
!                      w.schneider
!                      y.sun
!                      r.weiner
!
!
!                            for information regarding availability,
!                            contact:     s.c.jardin
!                                         pppl
!                                         p.o.box 451
!                                         princeton,n.j. 08543
!                                         609-243-2635
!                                         jardin@pppl.gov
!
!
!
!              *** index ***
!
!......section 1 main program
!
!......1.   program tsc
!......1.10    parmchk1
!......1.20    abexit
!
!
!......section 2 . . physics
!
!......2.10    stepon
!......2.20    advanc
!......2.30    advanp
!......2.35    advan23
!......2.40    hyper
!
!
!......section 3 . . auxilliary calculations
!
!......3.10    efield
!......3.11    cirequ
!......3.11a   feedeg
!......3.12    rscale
!......3.20    tstep
!......3.30    icalc
!......3.40    globe
!......3.41    xptcalc
!......        binary
!......3.42    svd1
!......3.50    resist
!......3.60    enerd2
!......3.70    boundc
!......3.71    divplat
!......3.72     divdis
!......3.80    auxdef
!......     entry ajdef
!......     entry gsdef
!......     entry prdef
!......3.90    limpsi
!......3.93    rsmod
!......3.94    shapenp
!......3.95    moments
!......3.96    appvolt2
!......        appvolt2
!......3.97    thyris
!......3.98    vprep
!
!......section 4 . . elliptic solver
!
!......4.10    elliptic
!......4.11    facr
!......4.60    bounda
!......4.61    pfluxb
!......4.70    grnfnc
!......4.71    grnfncv
!......4.80    gradgf
!......4.90    gvect
!
!
!......section 5 . . initial conditions
!
!......5.10    init
!......5.20    define
!......5.30    inpt
!......5.31    readcard
!......5.32 entry writcard
!......5.40    rstrt1
!......     entry rstrt2
!......5.50    defolt
!......5.60    setup
!......5.70    growth
!......5.71    fieldg
!......5.72    multipol
!......        fieldc
!......        multip2
!......5.80    zero
!......5.90    curinit
!......5.91    feednorm
!
!
!......section 6 . . output routines
!
!......6.10    eqwrite
!......     entry eqread
!......6.21    vplot
!......6.31    arrow
!......6.41    outpr1
!......6.70    outpts
!......6.75    wrpest
!......        groupcur
!......6.80    outpl
!......6.90    cplot
!......        xptplt
!......6.904.1 dplplt
!......6.905   sumplot
!......6.906   coildr
!......6.91    saveit
!.....      entry plotit2
!......        dif1
!......        disk
!......6.96    frscj
!.....      entry pltindx
!......6.97    timer
!......6.98    limits
!......6.99    fcwire
!
!......section 7 . . initial equilibrium
!
!......7.00    initeq
!......7.07    newj2
!......7.10    ff
!......7.20    pinterp
!
!......section 8 . . surface averaging
!
!......8.10    qcalc
!......8.20    flxvol
!......8.22    grap
!......8.23    magaxis
!......         axm2d2
!......8.30    geval
!......8.40    peval
!......8.50    eeval
!......8.60    reval
!......8.70    defcoef
!
!......section 9. . . transport
!
!......9.10    advsfa
!......9.20    trcdef
!......9.21    auxheat
!......9.22    ripple
!......9.23    ripplot
!......9.24    ripdef
!
!......section 10 . . . ideal mhd stability
!
!......10.1 balloon
!......10.3 intgrte
!......10.4 der
!......10.5 alfa
!......10.6 betaf
!......10.8 gama
!......     de
!......     step
!......     intrp
!
!......section 11 . . . for comparison with actual data
!
!......11.1 spdpbx
!......11.11 spdtftr
!......11.11 spdpbxm
!......11.2 spdd3d
!......11.3 chop
!......11.4 appvolt
!
!......section 12 . . . coil temperature,stresses,and inductance
!
!......12.1 temprise
!......12.2 rhovst
!......12.3 cucpvst
!......12.4 sscpvst
!......12.5 fndt
!......12.6 sfilfx
!
!......section 13 . . . pellet subroutines from w.houlburg
!
!......13.1    inicon
!......13.2    pelabl
!......13.3    pellet
!......13.4    pelqe
!......13.5    pelqf
!......13.6    pelrat
!......13.7    pelrk4
!......13.8    zincon
!
!
!......section 9 . . . wall impurity transport subroutines
!                wall impurity transport subroutines                    *
!
!......9.30    metricw -- metric parameters for impurity mesh
!......9.31    sprop   -- surface average impurity properties
!......9.32    advimp  -- impurity twodee transport
!......9.33    impdcf  -- impurity parrallel diffusion coeficients
!......9.34    ssolve  -- steady state impurity populations
!......9.35    radionz -- impurity radiation/ionization coeficients
!......9.36    lookup  -- looks up values from radiation tables
!......9.37    walld   -- driver for impurity/wall routines
!
!
!     list of subroutines in wall interaction part
!
!********************************************************************
!     part 1 - - wall subroutines                                   *
!********************************************************************
!
!     1a - wallin -- reads in all initial data
!     1b - lsfit  -- generates a least squares fit for poly functions
!     1c - ploads -- calculates limiter particle heat loads
!
!     1d - rloads -- calculates radiant wall heat loads
!
!********************************************************************
!     part 4 - - wall temperature subroutines                       *
!********************************************************************
!
!     4a1 - wallst  -- wall steady state temperature subroutine
!     4a2 - walltt  -- wall transient temperature driver subroutine
!     4b  - updtemp -- updates wall temperatures
!     4c  - nodegen -- generates nodes
!     4d  - evapkt  -- calculates vaporization parameters
!     4e  - updtcd  -- updates thermal conductivity
!     4f  - updrcp  -- updates heat capacity
!     4g  - d1deg1  -- single variable interpolation routine
!     4h  - tridag  -- tri-diagonal matrix solver
!     4i  - invert  -- periodic tri-diagonal matrix solver
!     4j  - print   -- threedee array printer
!     4k  - psat    -- saturation pressure function  (beryllium)
!     4l  - hfgs    -- saturation enthalpy function
!     4m  - erf     -- error function approximation
!......end
!
      program tsc
!
!.....main program:  princeton tokamak simulation code
!             ref :  "dynamic modeling of transport and positional
!                     control of tokamaks"  pppl-2258  1985
!
!            s.c.jardin,n.pomphrey,j.delucia,j.comput.phys,66,481(1986)
!
      USE trtsc, ONLY : wrtrst, suffix2
      USE input_mod
      USE CLINAM
      USE COMWOY
      USE NONCOR
      USE SAPROP
      USE SCR15
      USE tsc_fmcfm_transport_mod

      IMPLICIT NONE

#ifdef HAVE_MPI
!...  JCdebug april-03-2009
      include 'mpif.h'
      INTEGER :: myPE, totPEs, ierr
#endif

      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ::  n_newton, i_newton, i_tmp
!
!.......................................................................
!.....1.0 initialization
!......................................................................
!tss  call dropfile(0)
!tss  call ioscale
!.....set logical unit no for plot output
!     s100 is now a string variable
!......special output needed by IPP, Max PLank
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iwoy,ios9,ios1,ios12,nxpest,nzpest,ios13,ios14,ios19
      INTEGER ios20,iprnti,iplti,iquit
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dtnext,vtrans,secleft,tcpu,x1, preskpa
!============

#ifdef HAVE_MPI
!...  JCdebug april-03-2009
      CALL MPI_INIT( ierr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, totPEs, ierr )
!     write(*,*) "JCdebug tsc: myPE is ", myPE, " of total ", totPEs
#endif

      call getcarg(1, suffix, numargs)
#ifdef I_HAVE_MPI_DEBUG
      write(*,*) myPE, "JCdebug tsc: numargs ", numargs, " suffix=", suffix
#endif

      n35 = 96
!......special input needed for LSC when ilhcd=1, ifk=2
      nlhcdin = 4
!.....unit 59 is tty output
      nterm = 6
!     call msglink(nterm,1)
!.....remaining logical units
      nin = 95
      nout = 96
      nmiss = 33
      no167a = 7
      ndivlu = 8
      n13 = 13
      lujayp = 14
      luwire = 15
      luwirer= 16
      neqdsk = 17
      neqdska = 35
      neqdskg = 36
      neqplt = 18
      nsprsi = 21
      nsprso = 22
      nsc2 = 23
      nsc3 = 24
      nenin = 25
      neqin = 26
      neqou = 27
      nlhcd = 29
      nlhcd2 = 32
      nsc1 = 88
      nwoy = 28
      nwoyi = 31
      nwoyo = 30
      nlscgnu = 34
      iwoy = 0
      apl = 0.0d0
!
!
!.....append suffix info to file names
      call timedate(itime,idate,imach,isuffix)
!     ifilein(1:5) = 'input'
!     ifilein(6:6) = isuffix(1:1)
!     ifileou(1:6) = 'output'
!     ifileou(7:7) = isuffix(1:1)
!     ifiles1(1:5) = 'oust3'
!     ifiles1(6:6) = isuffix(1:1)
      if( numargs .lt. 1 ) then
         ifilein = 'input' // isuffix(1:1)
         ifileou = 'output' // isuffix(1:1)
         ifiles1 = 'oust3' // isuffix(1:1)
      else    
         ifilein = 'input' // '.' // trim(suffix)
         ifileou = 'output' // '.' // trim(suffix)
         ifiles1 = 'oust3' // '.' // trim(suffix)
      end if
      open(nout,file=trim(ifileou),status='unknown',iostat=ios9)
      write(nout,9876)
 9876 format(" hello world ")
!tss  call idtous("print",nout)
      open(nin,file=trim(ifilein),status='old',iostat=ios1)
      open(nsc1,file=trim(ifiles1),status='unknown',                     &  
          iostat=ios12)
      if(ios9 .ne. 0) then
          write(*,*) "JCdebug tsc: output cannot be opened ", ios9
      else if(ios1 .ne. 0) then
          write(*,*) "JCdebug tsc: input cannot be opened ", ios1
      else if(ios12 .ne. 0) then
          write(*,*) "JCdebug tsc: oust3 cannot be opened ", ios12
      else
!         write(*,*) "JCdebug tsc: all files are opened "
          write(*,*)
      end if
!
!tss  call gfsize(3,50000000b)
!
        call  tversion (nout)
        call  tversion (nterm)
!
!.....print dimension parameters
!
      call problem_size
      call allocate_arrays
      call parmchk1
      write (nout, 1002)
!
!.....initialize everything
!
!     for tsc-transp coupling
!     note:
!     acoef(4950)=time to run transp
!     acoef(4951)=integration interval in transp
!     acoef(4950) = 1000000.0
!     acoef(4951) = 0.0
!     acoef(4955)=hydrogen fraction used for outputing to the plasma state
!     note: acoef(113)= deuterium fraction in DT mix
!     acoef(4955) = 0.0
!     acoef(4956)=power in mw for each beam
!     acoef(4957)=number of newton iteration for glf
!     acoef(4957)=1
!     acoef(4958)=step size for calculating dchi/dT'
!     acoef(4958)=-1.0        ! negative ==> do not do linearization

!
      call init
      call read_input_mod
      if(acoef(4954) .gt. 1.0e-06) zimp=acoef(4954)
!
!
       write (nterm, 1002)   name
 1002  format (1x, 8a10)
!
!  create pest file
      if(ipest.eq.0) goto 6
      nxpest=nx
      nzpest=nz
      if(isym.eq.1) nzpest=2*nz-1
      if(ipest.eq.1) goto 617
      if(ipest.eq.2) go to 619
      if( numargs .lt. 1 ) then
         ieqplt = 'eqplt' // isuffix(1:1)
      else 
         ieqplt = 'eqplt' // '.' // trim(suffix)
      end if
!     ieqplt(1:5) = 'eqplt'
!     ieqplt(6:6) = isuffix(1:1)
      open(neqplt,file=trim(ieqplt),status='unknown',form='unformatted',  &  
     &    iostat=ios13)
617   continue
!.....ipest=1
      if( numargs .lt. 1 ) then
         ieqdsk = 'eqdsk' // isuffix(1:1)
      else
         ieqdsk = 'eqdsk' // '.' // trim(suffix)
      end if
!     ieqdsk(1:5) = 'eqdsk'
!     ieqdsk(6:6) = isuffix(1:1)
      open(neqdsk,file=trim(ieqdsk),status='unknown',form='unformatted',  &  
     &     iostat=ios14)
      if( numargs .lt. 1 ) then
         filename = 'eqdskasci'
      else
         filename = 'eqdskasci' // '.' // trim(suffix)
      end if
      open(neqdska,file=trim(filename),status='unknown',iostat=ios19)
      if( numargs .lt. 1 ) then
         filename = 'geqdsk'
      else
         filename = 'geqdsk' // '.' // trim(suffix)
      end if
      open(neqdskg,file=trim(filename),status='unknown',iostat=ios20)
       go to 6
  619 continue
!.....ipest=2
!     ieqplt(1:5) = 'eqplt'
!     ieqplt(6:6) = isuffix(1:1)
      if( numargs .lt. 1 ) then
         ieqplt = 'eqplt' // isuffix(1:1)
      else 
         ieqplt = 'eqplt' // '.' // trim(suffix)
      end if
      open(neqplt,file=trim(ieqplt),status='unknown',                    &  
     &    iostat=ios13)
6     continue
!
!
!     write(*,*) "JCdebug tsc: after 6 ineg=",ineg," ncycle=",ncycle
      if(ineg.ne.0 .and. ineg.ne.12) go to 40
      if(ncycle .lt. 0 ) go to 15
      if(ineg.ne.0) go to 16
!
      iprnti=0
      iplti =0
      if(irst1.eq.1 .and.iprnt.lt.nskipr) iprnti=1
      if(irst1.eq.1 .and. iplt .lt.nskipl) iplti=1
      if(irst1.eq.1) iwoy = 1
      if(irst1.eq.0) iskipsf = nskipsfi
!
!
!                              initialize for vv poloidal current calc.
        call  ivvpol
!----------------------------------------------
!@@@a
!     initialisierung des reglers mit technischer abtastzeit
      if ((acoef(296).eq.4._R8) .or. (acoef(290).eq.4._R8)) then
           call reginit (dtnext)
      endif
!@@@e-------------------------------------------
!
!.......................................................................
!.....2.0 start main time loop
!.......................................................................
!     write(*,*) "JCdebug tsc: start main time loop"
      call trtsci  !cj moved here to init trxpl server    !pTransp
!cj      if( acoef(4991) .gt. .0 .or.                                      &
!cj          acoef(4992) .gt. .0 .or.                                      &
!cj          acoef(4993) .gt. .0      )                                    &
!cj     &    call trxpl_init 
!cj moved to trtsci, before subroutine init_plasma_state is called nov-09-2010
    5 continue
!
!.....update print and plot counters
!
      iprnt = iprnt + 1
      iplt = iplt + 1
      iplt2 = iplt2 + 1
!
!.....update time and cycle counters
!
!cj   call trtsci                                         !pTransp
      times = times + dt*udst
      time = times*usdt
      kcycle = kcycle + 1
      dtsf = dt
      call param_time_mod(nin,times)

      if(lrswtch.ne.0) go to 309
!
!.....map flux surfaces every nskipsf cycles for isurf.ne.0
      iskipsf = iskipsf + 1
      if(iskipsf.le.nskipsf) go to 308
      iskipsf = 1
      call prdef
      if(isurf.eq.0 ) go to 309
!
      call flxvol(2)
      if(ineg.ne.0) go to 40
  308 continue
      if(isurf.eq.0) go to 309
!.....cj may 18, 2010 
!.....replace adp, ade, adn in defcoef with profiles
!.....imported from transp database with trxpl
      if((kcycle .gt. 1) .and.                                          &
         (acoef(4991) .gt. .0 .or.                                      &
          acoef(4992) .gt. .0 .or.                                      &
          acoef(4993) .gt. .0 .or.                                      &
          acoef(4994) .gt. .0     ))                                    &
     &    call import_profile_from_transp
!
      call defcoef
      call sprop
!
!.....check for jet source
      if(acoef(816).gt.0._R8) call jetsource
!
!.....check for impurity line radiation calculation
      if(iimp.ge.1) then
      if(kcycle.gt.0) go to 310
      call impinit
      call sprop
  310 continue
      call advimp1d
      call rloads1d
      endif
      if(irunaway.gt.0) call runaway_sub
!
      n_newton = 1
      if(int(acoef(4957)) .gt. 1) n_newton=int(acoef(4957))
      do i_newton=1, n_newton
!cj...july...19...2001      call trcdef(i_newton)
!
! Note: if itrmod > 20, call tsc_fmcfm_transport rather than trcdef
! to interface transport models through the FMCFM module
! rather than computing transport model directly
! Bateman, Jan 2009
!
      IF ( itrmod > 20 ) THEN
        CALL tsc_transport ( itrmod, i_tmp )
        WRITE (*,*) &
        'JCdebug tsc: call tsc_transport at kcycle=', kcycle, &
        'i_newton =', i_newton, &
        'iskipsf = ', iskipsf
!         STOP
      ELSE
        CALL trcdef(i_newton)
!       WRITE (*,*) &
!       'JCdebug tsc: call trcdef at kcycle=', kcycle, &
!       'i_newton =', i_newton, &
!       'iskipsf = ', iskipsf
      ENDIF
!
      if(ineg.ne.0) go to 40
      if(i_newton .gt. 1) call sprop
      call advsfa(i_newton)
      if(ineg.ne.0) go to 40
      enddo
! end of Newton iteration

      call curdrive
      if(iwall.eq.0 .or. kcycle.le.0) go to 309
      call walld
  309 continue
!
!.....calculate eta,psimin and psilim
      call resist
      call efield
!
!.....initial plots for restart run
      if(iprnti.eq.1) then
      call savprof
      iprnti = 0
            endif
      if(iplti .eq.1) then
      call outpl
      iplti = 0
            endif
!@@@a
!
!     regleraufruf mit technischer abtastzeit
      if ((acoef(296).eq.4._R8) .or. (acoef(290).eq.4._R8)) then
!
      if (iwoy .eq. -1) iwoy = 1
         if(times.gt.dtnext) then
!     ********************************
           if (iwoy .eq. 0) iwoy = -1
!     **** a. krause 27.09.90 (raw measurement calculation(fp))
           call rawmeas
!     **** a. krause 27.09.90 (controll algorithm (fb))
           call regler
           dtnext = times + dtabt
         endif
      endif
!@@@e-------------------------------------------
!
!
!.....advance 2d variables 1 step in time
!
!     write(6,*) "before stepon :", ineg
      call stepon
!
      call icalc
      call tstep
      dts = dt*udst
      if(dt.gt.dtmin) go to 8
      ineg=4
      write(nout,4777) kcycle,times,dts,dtmins,amach,acoef(21),          &  
     &              ekin,acoef(22)
 4777 format(" kcycle,times,dts,dtmins,amach,acoef(21),ekin,acoef(22) ",  &  
     &      /,i3,1p7e12.4)
      iprnt = nskipr + 1
    8 continue
!
!.....write pest file
      if(itpest.gt.nopest) go to 184
      if(times-tpest(itpest) .lt. 0._R8) go to 184
      if(ipest.gt.0) call wrpest
      if(ipest.gt.0 .and. acoef(3001).eq.2._R8) call ufwr2d
      iprnt = nskipr+1
      iplt  = nskipl+1
!     iplt2 = nskip2+1
      itpest = itpest+1
  184 continue
!
!.....check if time for saving plot data
!
      hfluxtot = hfluxtot + hfluxp
      itot = itot + 1
      call globea
!                             calculate diamag and vv poloidal current
        if (mod(itot,4) .eq. 0)   call  ivvpol
      if(iplt2.le.nskip2) go to 14
      hfluxav = hfluxtot/itot
      if(kcycle.le.nskip2+1) hfluxav = beamp(istart)
      hfluxtot = 0._R8
      itot = 0._R8
      iplt2 = 0
!     **** a. krause 27.9.90  (vertical stability test)
      if((acoef(296).eq.4 .or. acoef(290).eq.4) .and.                    &  
     &         noplot(24).le.0) call vforce
      if(ineg.ne.0) go to 13
      call globe
!     *******************************
!...  suppress plots before controller is called the first time
      if((acoef(296).eq.4 .or. acoef(290).eq.4).and.iwoy.le.0) go to 13
      call saveit
   13 continue
      call outpts
   14 continue
!
!.....check if time for output
!
      if(iprnt.le.nskipr) goto 10
!
      if(acoef(296).eq.2._R8) then
      vtrans = vtran*udsv
      write(nterm,67) kcycle,vtrans
      write(nout,67)  kcycle,vtrans
   67 format(1x,"kcycle=",i7,2x,"vt=",1pe16.9," volts")
      endif
!
      iprnt=1
!     call outpr1
      call savprof
        if (iwayne.gt.0)   call  fwout
   10 continue
!
!.....                     is it time to output a movie frame ?
!
        if (dtmovie.gt.1.0E-8_R8.and. times.ge.acoef(33))   then
                         acoef(33) = acoef(33)+dtmovie
!
!
!...write time and plasma boundary (R,Z) into file 'boundary'
!     open(unit=46,file='boundary',status='unknown',position='append')
!     if(isym .eq. 1) then
!     write(46,461) times,2*kmax-4
!     go to 459
!     endif
!     write(46,461) times,nthe
!459  do k=2,nthe+1
!     write(46,460) xw(k,npsit),zw(k,npsit)
!     enddo
!     if(isym .eq. 1) then
!     do k=1,kmax-3
!     write(46,460) xplot(1,kmax-k),-zplot(1,kmax-k)
!     enddo
!     endif
!460  format(2f20.4)
!461  format(1f20.4,i10)
!     close(46)
!
!

      if(imovie.eq.8) call cplot(g,8)
!cj feb07_2010      if(imovie.le.8) call cplot(psi,1)
!
      if(imovie.ge.10) call moviedata
                  mframe = mframe + 1
         if (mod(mframe,10).eq.0)   write (nterm,400)   mframe, times
                  write (nout,400)   mframe, times
  400             format (' movie frame no',i6,'  time='1pe10.3,' sec')
                                                  endif
!
!.....check if time for plotting
!
      if(iplt.le.nskipl) go to 15
      iplt = 1
   16   if ((imovie.eq.0 .or. imovie.ge.10))   call outpl
!
!.....write restart file
!
      if(irst2.ne.0 ) call rstrt2
      wrtrst=.true.                          !pTransp

   15 continue
!
! check if necessary to write ajphi to file jayph for disruption calc
!                        and forces to file wiref for vv load analysis
!
      if(iwayne.gt.0.and.iwayne.lt.ptwir                                 &  
     &              .and.(times-tjphi).ge.0._R8)   call  fcwire
!
!
!     handshake with transp, waiting for return
      call trtscm                           !pTransp
!

      if(ineg.ne.0) go to 40
!
!......................test for plasma current and shape quantities
!                      outside range specified by acoef(102)
      call tquit(iquit)
      if(iquit.ne.0) go to 18
!
!.....check if time-limit is almost exceeded
      call tremain(secleft)
      if(secleft .le. 120) then
      write(nterm,1018)
      write(nout ,1018)
 1018 format(" * * program stopped because cpu limit exceeded * *")
      go to 18
                          endif
!
!.....check if max time range is exceeded
!
      if(time.lt.(tpro(ntpts)-dt).and.times.lt.acoef(29)) go to 17
   18 if(kcycle .gt. ncycle) go to 30
      ncycle = kcycle
      iprnt = nskipr
      iplt = nskipl
      iplt2 = nskip2
   17 continue
!
!.....end of time loop
!
      if(kcycle.le.ncycle .and. ncycle .gt.0) go to 5
!     write(*,*) "JCdebug tsc: after main time loop"
!.......................................................................
!.....3.0  termination
!.......................................................................
      go to 30
   40 call abexit
      if((imovie.eq.0 .or. imovie.ge.10) .and. iplt.gt.1) call outpl
   30 continue
      if(kcycle.lt.0 .or. ncycle.lt.0) go to 34
      call outpts
!
!.....write restart file
!
!     restart is now written in trtscm
!     if(irst2.ne.0 .and. iplt.gt.2) call rstrt2
!     if(irst2.ne.0 .and. iplt.gt.2) call rstrt2

!@@@  w. woyke restart fuer regler schreiben
      if (irst2.ne.0.and.(acoef(290).eq.4.or.acoef(296).eq.4))           &  
     &   call restcon(dtnext)
      if(ipest.ne.0 .and. nopest.eq.0) call wrpest
      if(lrswtch.eq.0) call eqwrite
      if(ilhcd .eq. 1 .and. ifk.ne.2 .and. ineg.eq.0                     &  
     &      .and. isurf.eq.1) call lhwrite(nlhcd,0)
#ifndef NCAR_DUMMY
      if(((imovie.eq.0 .or. imovie.ge.10)).and.(noplot(10).le.0))        &  
     &         call plotit1
!... plot for vertical stability test...................................
!     **** a. krause 27.9 90
      if((acoef(296).eq.4 .or. acoef(290).eq.4) .and.                    &  
     &         noplot(24).le.0) call vforcepl
!     *********************************
!... plot for function parametrization..................................
!     **** a. krause  27.9.90
      if((acoef(296).eq.4 .or. acoef(290).eq.4) .and.                    &  
     &         noplot(25).le.0) call fpplot
!     *********************************
#endif
      if(ncycle.le.0) go to 34
!
      if(acoef(3001).eq.2._R8) call ufwr1d
!
!
#ifndef NCAR_DUMMY
!.....time history plots
      if((imovie.eq.0 .or. imovie.ge.10)) call plotit2
#endif
!
!
#ifndef NCAR_DUMMY
!.....special summary plot
      if(((imovie.eq.0 .or. imovie.ge.10)).and.(noplot(12).le.0)) then
      call sumplot
      if(idiv.ne.0) call sumplot2
      endif
#endif
!
   34 continue
#ifndef NCAR_DUMMY
      if(ncycle.eq.0) call sumplot
#endif
      if(((imovie.gt.0 .and. imovie.lt.10)).or.(noplot(23).gt.0))        &  
     &    go to 35
#ifndef NCAR_DUMMY
      call map(.285_R8,1._R8,.285_R8,1._R8,0._R8,1.0_R8,0._R8,1.0_R8)
      call setld(1._R8,40._R8,1,0,1,0)
#endif

!     call timer(s100)
      call frscj(1)
   35 continue
!
#ifndef NCAR_DUMMY
!.....produce graphics index
      call pltindx
#endif
!
!......destroy scratch files
#ifndef HAVE_MPI
      close(nsc1,status='delete')
      if(ifrst(4).eq.0) close(nsc2,status='delete')
      if(ifrst(6).ne.0) close(nsc3,status='delete')
#endif
!     call timer(nout)
!     call timer(nterm)
      call second(tcpu)
      x1 = tcpu/60._R8
      write(nout,7777) x1,global(125)
      write(nterm,7777) x1,global(125)
 7777 format("CPU time (min)",1pe12.4,"(segment)",1pe12.4,               &  
     &          "(cumulative)")
      if(imovie.ge.10) call movieclose
#ifndef NCAR_DUMMY
      call plote
#endif

#ifdef HAVE_MPI
!...  JCdebug april-03-2009
      CALL MPI_FINALIZE(ierr)
#endif
      if( acoef(4991) .gt. .0 .or.                                      &
          acoef(4992) .gt. .0 .or.                                      &
          acoef(4993) .gt. .0 .or.                                      &
          acoef(4994) .gt. .0      )                                    &
     &    call trxpl_kill 

      stop
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
