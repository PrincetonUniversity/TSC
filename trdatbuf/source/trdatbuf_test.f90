program trdatbuf_test

  !  test program for trdatbuf module (dmc May 2005 -- rough draft)
  !  zfile passed as first (and only) argument (awb June 2005)

  use trdatbuf_module
  implicit NONE

  type (trdatbuf) :: d    ! trdatbuf buffer object
  type (profget) :: t1    ! object to aid profile interpolation
  type (pwrget) :: zpwr   ! specification of desired power parameters

  !  File containing trdatbuf buffer object (from an NSTX TRANSP run)

  character*20 :: zfile = '112989Z11PH.DAT'
  character*3 :: ext

  integer :: tty = 6
  integer :: lunfile = 99, icdf
  integer :: errct = 0

  integer :: ierr,iwarn,idum,iacur,iavsf,iarbz,iaddr,iz
  integer :: narg, id, iargc, getarg, nlen ! command line arguments

  real*8 :: ztimea,ztimeb     ! start and stop times
  integer :: intime1,intime2  ! time grid sizes

  real*8 :: zta               ! target time for interpolation
  real*8 :: zitem,zrbza

  !  interpolation info:
  integer :: it1,it2
  real*8 :: zf1,zf2

  real*8 :: zr1a,zr2a         ! midplane bdy intercept radii at zta

  logical :: ipresent

  !----------------------------------------------------
  !  a coarse profile grid...
  
  integer, parameter :: nzones = 10

  real*8 :: xi_bdy(nzones+1)
  real*8 :: xi_cen(nzones+1)
  real*8 :: rmajmp(2*nzones+1)
  real*8 :: bmidp(2*nzones+1)

  real*8 :: shif,rcen,rmin,rmin1,shif1

  !----------------------------------------------------
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE  = 1.0d0
  real*8, parameter :: HALF = 0.5d0
  !----------------------------------------------------

  ! moments info

  !  boundary
  integer :: nmom, im ! total number and index of moments
  integer :: icirc ! flag for circular flux surfaces
  integer :: it  ! time index
  real*8 :: zt  ! time
  real*8 :: zf  ! time scale (offset) factor
  real*8, dimension(:,:), allocatable :: rmcb, ymcb

  ! equilibrium (boundary and interior surfaces)

  real*8, dimension(:,:,:), allocatable :: zrmc2,zymc2 ! asymmetric Fourier moments set at xi bdys
  real*8 :: zdrshaf(nzones+1) ! "Shafranov" shift of interior surfaces
                              ! relative to the boundary surface
  integer :: mj            ! flux surface array dimension for zrmc2,zymc2
  integer :: lcentr        ! index to magnetic axis
  integer :: nzp1          ! no. of surfaces including mag. axis
  integer :: j ! index

  !----------------------------------------------------
  !
  real*8 :: zRi, zYi, rpld ! location for TF ripple and value of ripple

  !----------------------------------------------------
  !
  integer :: ntsaw ! number of sawtooth times in trdatbuf
  integer :: inum ! index for sawtooth time
  real*8, dimension(:), allocatable :: ztimes !  sawtooth times in trdatbuf
  real*8 ::  ztime_prev, ztime ! preceding sawtooth time and next time
  real*8 :: ztime_pre, ztime_post ! sawtooth start and end time

  !----------------------------------------------------
  ! power-channel data
  !
  integer ::  nbeam, ib ! number of beams (or antenna) and index
  real*8, dimension(:), allocatable :: tonarr,toffarr ! on and off time arrays for data channel
  real*8 :: ton, toff !  min on and max off time for power-channel data
  real*8 :: zthresh ! on/off power threshhold
  real*8 :: zdtfix !  time to search from threshhold

  !----------------------------------------------------
  ! look for command line arguments - physics data file
  !
  narg = iargc()
  if(narg.lt.1)then
     write(tty,*) 'usage: trdatbuf_test <run-id>PH.DAT'
     stop
  endif

  id=getarg(1,zfile)

  !----------------------------------------------------
  !  (re)init trdatbuf object
  
  write(tty,*) ' initializing trdatbuf object... '

  call trdatbuf_init(d)   ! initialize new object (allocate buffers)

  call trdatbuf_reInit(d) ! re-initialize object (deallocate buffers then
  !  execute a call to trdatbuf_init(d)

  !---------------------------------------------------
  !  calls test presence of data -- should be .FALSE., nothing read yet.

  if(.not.tdb_present1(d,'CUR',iaddr)) then
     write(tty,*) ' CUR 1d data not in trdatbuf object (file not yet read)'
  else
     write(tty,*) ' ?? tdb_present1 unexpectedly returned TRUE'
     errct = errct + 1
  endif

  if(.not.tdb_present1(d,'VSF',iaddr)) then
     write(tty,*) ' VSF 1d data not in trdatbuf object (file not yet read)'
  else
     write(tty,*) ' ?? tdb_present1 unexpectedly returned TRUE'
     errct = errct + 1
  endif

  if(.not.tdb_present2(d,'TER',iaddr)) then
     write(tty,*) ' TER 2d data not in trdatbuf object (file not yet read)'
  else
     write(tty,*) ' ?? tdb_present2 unexpectedly returned TRUE'
     errct = errct + 1
  endif

  !---------------------------------------------------
  !  read data

  nlen = len_trim(zfile)
  ext = zfile(nlen-2:nlen)

  if(ext.eq."DAT" .or. ext.eq."dat")then
    open(unit=lunfile,file=zfile,status='old')
    call  trdatbuf_read_tsc(lunfile,d,ierr)
    close(lunfile)
  elseif(ext.eq."CDF" .or. ext.eq."cdf")then
    call trdatbuf_open_tsc(zfile,"r",icdf)
    call trdatbuf_rd_cdf_tsc(icdf,d,ierr)
    call trdatbuf_close_tsc(icdf)
  else
     write(tty,*)  ext, ' extension not recognized'
     write(tty,*) 'Please supply a .DAT (ascii) of .cdf (NetCDF) extension'
      stop
  endif

  if(ierr.ne.0) then
     write(tty,*) ' ?? error reading data: ', trim(zfile)
     stop
  else
     write(tty,*) ' '
     write(tty,*) ' file read successful: ', trim(zfile)
     write(tty,*) ' '
  endif

  !---------------------------------------------------
  !  calls test presence of data -- should be .TRUE., data was read.

  if(.not.tdb_present1(d,'CUR',iaddr)) then
     write(tty,*) ' ?? tdb_present1 unexpectedly returned FALSE'
     errct = errct + 1
  else
     write(tty,*) ' CUR 1d data detected in trdatbuf object.'
  endif

  if(.not.tdb_present1(d,'VSF',iaddr)) then
     write(tty,*) ' ?? tdb_present1 unexpectedly returned FALSE'
     errct = errct + 1
  else
     write(tty,*) ' VSF 1d data detected in trdatbuf object.'
  endif

  if(.not.tdb_present1(d,'RBZ',iaddr)) then
     write(tty,*) ' ?? tdb_present1 unexpectedly returned FALSE'
     errct = errct + 1
  else
     write(tty,*) ' RBZ 1d data detected in trdatbuf object.'
  endif

  write(tty,*) ' '

  if(.not.tdb_present2(d,'TER',iaddr)) then
     write(tty,*) ' ?? tdb_present2 unexpectedly returned FALSE'
     errct = errct + 1
  else
     write(tty,*) ' TER 2d data detected in trdatbuf object.'
  endif

  if(.not.tdb_present2(d,'NER',iaddr)) then
     write(tty,*) ' ?? tdb_present2 unexpectedly returned FALSE'
     errct = errct + 1
  else
     write(tty,*) ' NER 2d data detected in trdatbuf object.'
  endif

  !---------------------------------------------------
  !  summary timebase information

  call tdb_tlims(d,ztimea,ztimeb)
  write(tty,*) ' time range of trdatbuf data: ',ztimea,ztimeb

  call tdb_ntimes(d,intime1,intime2)
  write(tty,*) '   no. of pts in f(t) timebase:   ',intime1
  write(tty,*) '   no. of pts in f(x,t) timebase: ',intime2

  zta = (ztimea + ztimeb)*HALF
  write(tty,*) ' time for profile testing:'
  write(tty,*) '   zta = ',zta,' seconds.'

  !---------------------------------------------------
  !  1d data lookup at time zta
  !  note: for efficiency, key lookup should only be done once
  !    per trdatbuf object read, and, one time zone lookup should
  !    be amortized over several 1d time interpolation calls...

  ipresent = tdb_present1(d,'CUR',iacur)   ! iacur returned
  ipresent = tdb_present1(d,'VSF',iavsf)   ! iavsf returned
  ipresent = tdb_present1(d,'RBZ',iarbz)   ! iarbz returned

  call tdb_lookup1(d,zta,it1,it2,zf1,zf2)  ! lookup up time zta
  !  interpolation indices & factors returned: it1,it2,zf1,zf2

  write(tty,*) ' at time: ',zta

  zitem = tdb_fintrp1(d,iacur,it1,it2,zf1,zf2)
  write(tty,*) ' plasma current (A):  ',zitem

  zitem = tdb_fintrp1(d,iavsf,it1,it2,zf1,zf2)
  write(tty,*) ' surface voltage (V): ',zitem

  zitem = tdb_fintrp1(d,iarbz,it1,it2,zf1,zf2)
  write(tty,*) ' R*B_vacuum (cm*T):   ',zitem

  ! save this item for future reference...
  zrbza = zitem 

  !---------------------------------------------------
  call tdb_rmp_bdy(d,zta,zr1a,zr2a)

  write(tty,*) ' at time ',zta,' (s):'
  write(tty,*) '    inner midplane boundary intercept radius (cm): ',zr1a
  write(tty,*) '    outer midplane boundary intercept radius (cm): ',zr2a
  write(tty,*) ' '

  !---------------------------------------------------
  !  profile testing...

  !  INITIALIZATION call -- to support data mapping, this call must
  !  give an upper limit on the number of zones that will ever occur in
  !  a run; it is called once at the start or restart of a simulation...

  call tdb_symini(d,nzones)  ! nzones does not vary in this test program...

  !  sqrt(Phi/Philim) [(normalized toroidal flux)**0.5] radial grid...

  xi_bdy(1)=ZERO
  do iz=1,nzones
     xi_bdy(iz+1)=xi_bdy(iz)+ONE/nzones  ! .1, .2, ... , 1.0
     xi_cen(iz)=HALF*(xi_bdy(iz)+xi_bdy(iz+1))
  enddo
  xi_cen(nzones+1)=xi_bdy(nzones+1)+HALF/nzones  ! half zone beyond boundary...

  !  SIMPLIFICATIONS:
  !  instead of using a real equilibrium, we'll just a shafranof shift
  !  on the midplane; instead of using a real field we'll just use the
  !  vacuum field

  !  at time zta...
  shif=0.15
  rcen = HALF*(zr1a+zr2a)
  rmin = HALF*(zr2a-zr1a)

  !  major radius grid on midplane -- note indices:
  !    1 --> inner plasma midplane boundary intercept
  !    nzones+1 --> magnetic axis
  !    2*nzones+1 --> outer plasma midplane boundary intercept

  rmajmp(nzones+1)=rcen + shif*rmin
  do iz=1,nzones
     rmin1=xi_bdy(iz+1)*rmin
     shif1=(ONE-xi_bdy(iz+1)**2)*shif
     rmajmp(nzones+1-iz) = rcen + shif1*rmin - rmin1
     rmajmp(nzones+1+iz) = rcen + shif1*rmin + rmin1
  enddo

  !  vacuum field:

  bmidp = zrbza/rmajmp

  !  ** profget structure for time zta **

  call tdb_profget_init(t1,nzones,.FALSE.)
  t1%time = zta
  t1%xibdys = xi_bdy
  t1%rmajmp = rmajmp
  t1%bmidp = bmidp

  t1%ibdy = .FALSE.   ! interpolate to zone centers; bdy values set to
                      ! half the sum of values of the adjacent zones

  t1%item = 'NER'     ! identifier for electron density...

  call tdb_profin(d,t1,ierr)
  if(ierr.ne.0) then
     write(tty,*) ' ?? tdb_profin "NER" call failed.'
     errct = errct + 1
  else
     write(tty,*) ' '
     write(tty,*) ' "NER" (electron density) at zone centers & bdys (cm**-3):'
     write(tty,*) ' ...(data mapped to zone centers; bdys interpolated)'
     do iz=1,nzones+1
        write(tty,8888) iz,xi_bdy(iz),t1%data_zb(iz),xi_cen(iz),t1%data_zc(iz)
     enddo
     ! t1%data_zb(1:nzones+1) is also available -- boundary centered data
  endif

  ! now, get the electron temeperature, with a different interpolation option:

  t1%ibdy = .TRUE.    ! interpolate to zone bdys; center values set to 
                      ! half the sum of values of the adjacent bdys

  t1%item = 'TER'

  call tdb_profin(d,t1,ierr)
  if(ierr.ne.0) then
     write(tty,*) ' ?? tdb_profin "NER" call failed.'
     errct = errct + 1
  else
     write(tty,*) ' '
     write(tty,*) ' "TER" (electron temperature) at zone centers & bdys (eV):'
     write(tty,*) ' ...(data mapped to zone bdys; centers interpolated)'
     do iz=1,nzones+1
        write(tty,8888) iz,xi_bdy(iz),t1%data_zb(iz),xi_cen(iz),t1%data_zc(iz)
     enddo
     ! t1%data_zc(1:nzones+1) is also available -- zone centered data
  endif

8888 format(1x,i4,2(2x,2(1x,1pe13.6)))
  write(tty,*) ' '

  !  cleanup... avoid memory leaks!

  call tdb_profget_free(t1)

  !---------------------------------------------------
  !  more queries as to file contents

  ipresent = tdb_logchk_special(d,'RPL',iwarn)
  write(tty,*) ' -> TF field ripple vs. (R,Z) in data file: ',ipresent
  errct = errct + iwarn
  if(ipresent)then
    zRi = zr2a
    zYi = 0.
    rpld=tdb_get_rpld(d,zRi,zYi,ierr)
    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_get_rpld call failed.'
       errct = errct + 1
    else
      write(tty,*) " TF field ripple ln(B~/B) = ", rpld, " at R,Z = ", zRi,zYi
    endif
  endif

  ipresent = tdb_logchk_special(d,'SAW',iwarn)
  write(tty,*) ' -> sawtooth event times in data file: ',ipresent
  errct = errct + iwarn
  if(ipresent)then
    call tdb_num_sawtimes(d,ntsaw)
    write(tty,*) 'number of sawtooth times in trdat buffer = ', ntsaw

    inum=ntsaw
    allocate(ztimes(ntsaw))
    call tdb_sawtimes(d,ztimes,inum)
    write(tty,*) 'sawtooth times = ', ztimes(1:inum)

    ztime_prev = ztimes(1)
    call tdb_find_next_sawtime(d,ztime_prev,ztime)
    write(tty,*) 'sawtooth time following ', ztime_prev, ' is ', ztime

    ztime_pre = ztimes(1)
    call tdb_post_sawtime(d,ztime_pre,ztime_post,iwarn)
    write(tty,*) ' for sawtooth start time = ', ztime_pre, ' end time = ', ztime_post

    call tdb_post_peltime(d,ztime_pre,ztime_post,iwarn)
    write(tty,*) ' for pellet start time = ', ztime_pre, ' end time = ', ztime_post
  endif

  !--------------------------------------------
  !
  !   ECH/ECCD antenna powers

  ipresent = tdb_logchk_special(d,'ECP',iwarn)
  write(tty,*) ' -> ECH/ECCD antenna powers in data file: ',ipresent
  errct = errct + iwarn
  if(ipresent)then
    nbeam = tdb_nchan_find(d,"EC")
    write(tty,*) 'Number of ECH/ECCD antennae: ', nbeam

    zthresh = 0.10
    zdtfix = 10.
    allocate(tonarr(nbeam), toffarr(nbeam))
    call tdb_onoff_atime(d,"EC",zthresh,zdtfix,ton   ,toff   ,ierr)
    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_onoff_atime call failed.'
       errct = errct + 1
    else
       write(tty,*) "EC min on/max off time = ", ton, toff
    endif
    call tdb_onoff_times(d,"EC",zthresh,zdtfix,tonarr,toffarr,ierr)
    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_onoff_times call failed.'
       errct = errct + 1
    endif

    call tdb_pwrget_init(zpwr)
    zpwr%item = 'ECP'
    zpwr%ztime1 = ton
    zpwr%ztime2 = toff
    zpwr%pweight = .TRUE.
    zpwr%nbeam = nbeam
    zpwr%tbon(1:nbeam) = tonarr(1:nbeam)
    zpwr%tboff(1:nbeam) = toffarr(1:nbeam)

    call tdb_pwrdata_avg(d,zpwr,ierr)

    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_pwrdata_avg "ECP" call failed.'
       errct = errct + 1
    else
       write(tty,*) 'Average Power over time interval ', ton, toff
       do ib=1,nbeam
         write(tty,*) 'Antenna ', ib, ' : ', zpwr%zparam(ib)
       enddo
    endif
  endif

  !--------------------------------------------
  !
  !  ICRF antenna powers

  ipresent = tdb_logchk_special(d,'RFP',iwarn)
  write(tty,*) ' -> ICRF antenna powers in data file: ',ipresent
  errct = errct + iwarn
  if(ipresent)then
    nbeam = tdb_nchan_find(d,"RF")
    write(tty,*) 'Number of ICRF antennae: ', nbeam

    zthresh = 0.10
    zdtfix = 10.
    allocate(tonarr(nbeam), toffarr(nbeam))
    call tdb_onoff_atime(d,"RF",zthresh,zdtfix,ton   ,toff   ,ierr)
    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_onoff_atime call failed.'
       errct = errct + 1
    else
       write(tty,*) "RF min on/max off time = ", ton, toff
    endif
    call tdb_onoff_times(d,"RF",zthresh,zdtfix,tonarr,toffarr,ierr)
    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_onoff_times call failed.'
       errct = errct + 1
    endif

    call tdb_pwrget_init(zpwr)
    zpwr%item = 'RFP'
    zpwr%ztime1 = ton
    zpwr%ztime2 = toff
    zpwr%pweight = .TRUE.
    zpwr%nbeam = nbeam
    zpwr%tbon(1:nbeam) = tonarr(1:nbeam)
    zpwr%tboff(1:nbeam) = toffarr(1:nbeam)

    call tdb_pwrdata_avg(d,zpwr,ierr)

    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_pwrdata_avg "RFP" call failed.'
       errct = errct + 1
    else
       write(tty,*) 'Average Power over time interval ', ton, toff
       do ib=1,nbeam
         write(tty,*) 'Antenna ', ib, ' : ', zpwr%zparam(ib)
       enddo
    endif
    !
    !  ICRF antenna frequencies
    !
     ipresent = tdb_logchk_special(d,'RFF',iwarn)
     write(tty,*) '    -> ICRF frequencies vs time: ',ipresent
     errct = errct + iwarn
     if(ipresent)then
       zpwr%item = 'RFF'

       call tdb_pwrdata_avg(d,zpwr,ierr)

       if(ierr.ne.0) then
          write(tty,*) ' ?? tdb_pwrdata_avg "RFF" call failed.'
          errct = errct + 1
       else
          write(tty,*) 'Average frequencies over time interval ', ton, toff
          do ib=1,nbeam
            write(tty,*) 'Antenna ', ib, ' : ', zpwr%zparam(ib)
          enddo
       endif
     endif
  endif

  !--------------------------------------------
  !
  !  Lower Hybrid antenna powers

  ipresent = tdb_logchk_special(d,'LHP',iwarn)
  write(tty,*) ' -> Lower Hybrid antenna powers in data file: ',ipresent
  errct = errct + iwarn
  if(ipresent)then
    nbeam = tdb_nchan_find(d,"LH")
    write(tty,*) 'Number of Lower Hybrid antennae: ', nbeam

    zthresh = 0.10
    zdtfix = 10.
    allocate(tonarr(nbeam), toffarr(nbeam))
    call tdb_onoff_atime(d,"LH",zthresh,zdtfix,ton   ,toff   ,ierr)
    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_onoff_atime call failed.'
       errct = errct + 1
    else
       write(tty,*) "LH min on/max off time = ", ton, toff
    endif
    call tdb_onoff_times(d,"LH",zthresh,zdtfix,tonarr,toffarr,ierr)
    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_onoff_times call failed.'
       errct = errct + 1
    endif

    call tdb_pwrget_init(zpwr)
    zpwr%item = 'LHP'
    zpwr%ztime1 = ton
    zpwr%ztime2 = toff
    zpwr%pweight = .TRUE.
    zpwr%nbeam = nbeam
    zpwr%tbon(1:nbeam) = tonarr(1:nbeam)
    zpwr%tboff(1:nbeam) = toffarr(1:nbeam)

    call tdb_pwrdata_avg(d,zpwr,ierr)

    if(ierr.ne.0) then
       write(tty,*) ' ?? tdb_pwrdata_avg "LHP" call failed.'
       errct = errct + 1
    else
       write(tty,*) 'Average Power over time interval ', ton, toff
       do ib=1,nbeam
         write(tty,*) 'Antenna ', ib, ' : ', zpwr%zparam(ib)
       enddo
    endif
  endif

  !--------------------------------------------
  !
  !  NB data

  ipresent = tdb_logchk_special(d,'NB2',iwarn)
  write(tty,*) ' -> neutral beam injection powers in data file: ',ipresent
  errct = errct + iwarn

  if(ipresent) then
     nbeam = tdb_nchan_find(d,"NB")
     write(tty,*) 'Number of Neutral Beams: ', nbeam

     zthresh = 0.01
     zdtfix = 10.
     allocate(tonarr(nbeam), toffarr(nbeam))

     call tdb_onoff_atime(d,"NB",zthresh,zdtfix,ton   ,toff   ,ierr)
     if(ierr.ne.0) then
        write(tty,*) ' ?? tdb_onoff_atime call failed.'
        errct = errct + 1
     else
        write(tty,*) "NB min on/max off time = ", ton, toff
     endif

     call tdb_onoff_times(d,"NB",zthresh,zdtfix,tonarr,toffarr,ierr)
     if(ierr.ne.0) then
        write(tty,*) ' ?? tdb_onoff_times call failed.'
        errct = errct + 1
     endif

     call tdb_pwrget_init(zpwr)
     zpwr%ztime1 = ton
     zpwr%ztime2 = toff
     zpwr%pweight = .TRUE.
     zpwr%nbeam = nbeam
     zpwr%tbon(1:nbeam) = tonarr(1:nbeam)
     zpwr%tboff(1:nbeam) = toffarr(1:nbeam)

    !--------------------------------------------
    !
    !  NB Power data

     ipresent = tdb_logchk_nbi(d,'PWR',iwarn)
     write(tty,*) '    NBI powers vs. time: ',ipresent
     errct = errct + iwarn
     if(ipresent) then
       zpwr%item = 'PWR'

       call tdb_pwrdata_avg(d,zpwr,ierr)

       if(ierr.ne.0) then
          write(tty,*) ' PWR tdb_pwrdata_avg call failed.'
          errct = errct + 1
       else
          write(tty,*) 'Average NB Power over time interval ', ton, toff
          do ib=1,nbeam
          write(tty,*) 'Neutral Beam ', ib, ' : ', zpwr%zparam(ib)
          enddo
       endif
     endif
    !--------------------------------------------
    !
    !  NB Voltage data

  
     ipresent = tdb_logchk_nbi(d,'VLT',iwarn)
     write(tty,*) '    NBI voltages vs. time: ',ipresent
     errct = errct + iwarn
     if(ipresent) then
       zpwr%item = 'VLT'

       call tdb_pwrdata_avg(d,zpwr,ierr)

       if(ierr.ne.0) then
          write(tty,*) ' ?? tdb_pwrdata_avg "VLT" call failed.'
          errct = errct + 1
       else
          write(tty,*) 'Average NB Voltages over time interval ', ton, toff
          do ib=1,nbeam
          write(tty,*) 'Neutral Beam ', ib, ' : ', zpwr%zparam(ib)
          enddo
       endif
     endif
  
    !--------------------------------------------
    !
    !  NB Full Energy data

     ipresent = tdb_logchk_nbi(d,'FUL',iwarn)
     write(tty,*) '    NBI full energy fractions vs. time: ',ipresent
     errct = errct + iwarn
     if(ipresent) then
       zpwr%item = 'FUL'

       call tdb_pwrdata_avg(d,zpwr,ierr)

       if(ierr.ne.0) then
          write(tty,*) ' ?? tdb_pwrdata_avg "FUL" call failed.'
          errct = errct + 1
       else
          write(tty,*) 'Average NB full energy fraction over time interval ', ton, toff
          do ib=1,nbeam
          write(tty,*) 'Neutral Beam ', ib, ' : ', zpwr%zparam(ib)
          enddo
       endif
     endif
  
    !--------------------------------------------
    !
    !  NB Half Energy data

     ipresent = tdb_logchk_nbi(d,'HLF',iwarn)
     write(tty,*) '    NBI half energy fractions vs. time: ',ipresent
     errct = errct + iwarn
     if(ipresent) then
       zpwr%item = 'HLF'

       call tdb_pwrdata_avg(d,zpwr,ierr)

       if(ierr.ne.0) then
          write(tty,*) ' ?? tdb_pwrdata_avg "HLF" call failed.'
          errct = errct + 1
       else
          write(tty,*) 'Average NB half energy fraction over time interval ', ton, toff
          do ib=1,nbeam
          write(tty,*) 'Neutral Beam ', ib, ' : ', zpwr%zparam(ib)
          enddo
       endif
     endif
  endif

  !--------------------------------------------
  !
  !  plasma boundary position and shape...
  !
  ! First find number of moments
  !

  call tdb_find_nmoms(d,nmom,icirc)
  if(icirc.eq.1) then
     write(tty,*) ' Circular Flux Surfaces being used'
  endif

  ! for time zt get time index, it and scale (offset) factor, zf

      zt = zta
  call tdb_bdy_timefac(d,zt,it,zf)
    !
  allocate(rmcb(0:nmom,2), ymcb(0:nmom,2) )
  !  do im=0,nmom
  !     rmcb(im,1)=tdb_bdy_getmom(d,it,zf,im,TDB_MOMS_RCOS) ! Rcos mom.
  !     rmcb(im,2)=tdb_bdy_getmom(d,it,zf,im,TDB_MOMS_RSIN)
  !     ymcb(im,1)=tdb_bdy_getmom(d,it,zf,im,TDB_MOMS_ZCOS)
  !     ymcb(im,2)=tdb_bdy_getmom(d,it,zf,im,TDB_MOMS_ZSIN) ! Zsin mom.
  !  enddo

  call tdb_bdy_getmoms(d,it,zf,rmcb,ymcb,nmom)

  write(tty,*) ' Plasma Boundary moments '
  write(tty,*) ' m            rmc            rms            zmc            zms'
  do im=0,nmom
     write(tty,"(1x,i4,1p4e15.6)") im, rmcb(im,1), rmcb(im,2), ymcb(im,1), ymcb(im,2)
  enddo

  !
  !--------------------------------------------
  !
  !  Equilibrium Geometry
  ! using the MMX data, fetch the current equilibrium geometry.
  ! this routine MUST NOT be called unless the MMX data exists.
  !

  ipresent = tdb_logchk_special(d,'MMX',iwarn)
  write(tty,*) ' -> MHD equilibrium vs. time in data file: ',ipresent
  write(tty,*) '    (if F, FALSE, equilibrium bdy vs. time is available).'
  errct = errct + iwarn
  
  if(ipresent) then
     nzp1 = nzones + 1
     lcentr = 1  ! axis index (can differ from 1 if needed)
     mj = lcentr + nzones
     allocate(zrmc2(mj,0:nmom,2), zymc2(mj,0:nmom,2) )

     call tdb_getmmx(d,zt,xi_bdy,zrmc2,zymc2,mj,nmom,lcentr,nzp1,zdrshaf)

     write(tty,*) ' Equilibrium Geometry moments '
     write(tty,*)
     do j=lcentr,lcentr+nzones
       write(tty,*)
       write(tty,*) ' Surface index ', j,  '  Surfaces from axis ', j-lcentr
       write(tty,*) '  xi = ', xi_bdy(j-lcentr+1)
       write(tty,*) '  Shafranov shift = ', zdrshaf(j)
       write(tty,*)
       write(tty,*) ' m          rmc            rms            zmc            zms'
       do im=0,nmom
          write(tty,"(1x,i4,1p4e15.6)") im, zrmc2(j,im,1), zrmc2(j,im,2), zymc2(j,im,1), zymc2(j,im,2)
       enddo
     enddo

  endif

  write(tty,*) '-----------------------'
  write(tty,*) ' trdatbuf_test exit, errct = ',errct

end program trdatbuf_test
