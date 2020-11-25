subroutine ne_from_trdatbuf(times,xi_bdy,rmajmp,bmidp,ner,nzones,ierr)

  use trdatbuf_obj_tsc
  use trdatbuf_aux
  use clinam, ONLY : npsi
  implicit NONE

  type (trdatbuf) :: d    ! trdatbuf buffer object
  type (profget) :: t1    ! object to aid profile interpolation

  !  File containing trdatbuf buffer object (from an NSTX TRANSP run)

  character*20 :: zfile = 'TRDATBUF_PH'
  character*3 :: ext

  integer :: tty = 6
  integer :: lunfile = 99, icdf
  integer :: errct = 0

  integer :: ierr,iaddr,iz

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
  
  integer, intent(in) :: nzones 
  real*8, intent(in) :: times

  real*8 :: xi_bdy(nzones+1)
  real*8 :: rmajmp(2*nzones+1)
  real*8 :: bmidp(2*nzones+1)
  real*8 :: ner(nzones+1)

  !----------------------------------------------------
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE  = 1.0d0
  real*8, parameter :: HALF = 0.5d0
  !----------------------------------------------------

  integer :: j ! index
  real*8 :: fac, t0=0.0, ner0=0.05e+13

  integer :: first_time = 0
  integer :: istat
  character*120 :: filename
  save ztimea, ztimeb, intime1, intime2, d

  !----------------------------------------------------
 if (first_time .eq. 0) then 
  write(tty,*) ' initializing trdatbuf object... '

  call trdatbuf_init(d)   ! initialize new object (allocate buffers)

  call trdatbuf_reInit(d) ! re-initialize object (deallocate buffers then

  !---------------------------------------------------
  !  read data

    ext = "DAT"
    filename = trim(zfile) // "." // trim(ext)
    open(unit=lunfile,file=filename,status='old',iostat=istat)
    if (istat .eq. 0) then 
    call  trdatbuf_read_tsc(lunfile,d,ierr)
    close(lunfile)
    else
      close(lunfile)
      ext = "dat"
      filename = trim(zfile) // "." // trim(ext)
      open(unit=lunfile,file=filename,status='old',iostat=istat)
      if (istat .eq. 0) then
      call  trdatbuf_read_tsc(lunfile,d,ierr)
      close(lunfile)
      else
        close(lunfile)
        ext = "CDF"
        filename = trim(zfile) // "." // trim(ext)
        call trdatbuf_open_tsc(filename,"r",icdf)
        if(icdf .ne. 0) then
        call trdatbuf_rd_cdf_tsc(icdf,d,ierr)
        call trdatbuf_close_tsc(icdf)
        else
          ext = "cdf"
          filename = trim(zfile) // "." // trim(ext)
          call trdatbuf_open_tsc(filename,"r",icdf)
          if(icdf .ne. 0) then
          call trdatbuf_rd_cdf_tsc(icdf,d,ierr)
          call trdatbuf_close_tsc(icdf)
          else
          write(tty,*)  ext, ' extension not recognized'
          write(tty,*) 'Please supply a .DAT (ascii) of .cdf (NetCDF) extension'
          stop
          endif
        endif
      endif
    endif

  if(ierr.ne.0) then
     write(tty,*) ' ?? error reading data: ', trim(filename)
     stop
  else
     write(tty,*) ' '
     write(tty,*) ' file read successful: ', trim(filename)
     write(tty,*) ' '
  endif

  !---------------------------------------------------
  !  calls test presence of data -- should be .TRUE., data was read.

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

  !  INITIALIZATION call -- to support data mapping, this call must
  !  give an upper limit on the number of zones that will ever occur in
  !  a run; it is called once at the start or restart of a simulation...

  call tdb_symini(d,npsi)  ! nzones does not vary in this test program...

 first_time = 1
 endif

  zta = max(times,ztimea)
  fac = (times-t0)/(zta-t0)

  !---------------------------------------------------
  !  note: for efficiency, key lookup should only be done once
  !    per trdatbuf object read, and, one time zone lookup should
  !    be amortized over several 1d time interpolation calls...

  call tdb_lookup1(d,zta,it1,it2,zf1,zf2)  ! lookup up time zta
  !  interpolation indices & factors returned: it1,it2,zf1,zf2

! write(tty,*) ' at time: ',zta

  !---------------------------------------------------
  !  profile testing...

  !  sqrt(Phi/Philim) [(normalized toroidal flux)**0.5] radial grid...

  !  ** profget structure for time zta **

  call tdb_profget_init(t1,nzones,.FALSE.)
  t1%time = zta
  t1%xibdys(:nzones+1) = xi_bdy(:nzones+1)
  t1%rmajmp(:2*nzones+1) = rmajmp(:2*nzones+1)
  t1%bmidp(:2*nzones+1) = bmidp(:2*nzones+1)

  t1%ibdy = .FALSE.   ! interpolate to zone centers; bdy values set to
                      ! half the sum of values of the adjacent zones

  t1%item = 'NER'     ! identifier for electron density...

  call tdb_profin(d,t1,ierr)
  if(ierr.ne.0) then
     write(tty,*) ' ?? tdb_profin "NER" call failed.'
     errct = errct + 1
  else
     ner = t1%data_zb
!    write(tty,*) ' '
!    write(tty,*) ' "NER" (electron density) at zone centers & bdys (cm**-3):'
!    write(tty,*) ' ...(data mapped to zone centers; bdys interpolated)'
!    do iz=1,nzones+1
!       write(tty,8888) iz,xi_bdy(iz),t1%data_zb(iz),xi_cen(iz),t1%data_zc(iz)
!    enddo
     ! t1%data_zb(1:nzones+1) is also available -- boundary centered data
  endif

8888 format(1x,i4,2(2x,2(1x,1pe13.6)))

!    temporary, for testing
     do j = 1, nzones+1
     ner(j) = (1.0-fac)*ner0+fac*ner(j)
     enddo

  !  cleanup... avoid memory leaks!
  
  call tdb_profget_free(t1)

  !---------------------------------------------------

  return
end subroutine ne_from_trdatbuf

