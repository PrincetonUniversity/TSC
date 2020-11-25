      subroutine spdnstx
!			Rewrite to use new NSTX shot data files.
!			ROS.  20 Nov 2009 
!
!                       Generalize to also work for DIII-D
!                       SCJ   Feb 2010
!!!!!!!      USE ezcdf
!!!!!!!#define TSC_USE_NETCDF
      USE CLINAM, ONLY : numargs, filename, suffix, idata
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 stat
!============
      save
 
      REAL(KIND=r8), dimension(:),   allocatable :: time
 
      character*50,   dimension(:),   allocatable :: psupply_name
      REAL(KIND=r8), dimension(:,:), allocatable :: psupply_current
 
      character*50,   dimension(:),   allocatable :: halosensor_name
      REAL(KIND=r8), dimension(:,:),  allocatable :: halosensor_current
 
      character*50,   dimension(:),   allocatable :: flux_name
      integer,        dimension(:),   allocatable :: flux_active
      REAL(KIND=r8), dimension(:),   allocatable :: flux_r
      REAL(KIND=r8), dimension(:),   allocatable :: flux_z
      REAL(KIND=r8), dimension(:,:), allocatable :: flux_data
 
      character*50,   dimension(:),   allocatable :: mirnov_name
      integer,        dimension(:),   allocatable :: mirnov_active
      REAL(KIND=r8), dimension(:),   allocatable :: mirnov_r
      REAL(KIND=r8), dimension(:),   allocatable :: mirnov_z
      REAL(KIND=r8), dimension(:),   allocatable :: mirnov_tangle
      REAL(KIND=r8), dimension(:),   allocatable :: mirnov_pangle
      REAL(KIND=r8), dimension(:,:), allocatable :: mirnov_data
      REAL(KIND=r8), dimension(:), allocatable :: midplane_nedl
      REAL(KIND=r8), dimension(:), allocatable :: nbi_1_power
      REAL(KIND=r8), dimension(:), allocatable :: diamagnetic_flux
 
      character*50,   dimension(:),   allocatable :: current_name
      integer,        dimension(:),   allocatable :: current_active
      REAL(KIND=r8), dimension(:,:), allocatable :: current_data
      REAL(KIND=r8), dimension(:,:), allocatable :: mpts_Te
      REAL(KIND=r8), dimension(:), allocatable :: mpts_radius
      REAL(KIND=r8), dimension(:,:), allocatable :: mpts_ne
      REAL(KIND=r8), dimension(:), allocatable :: maxte
       REAL(KIND=r8)                              :: amaxtemp, tffactor
 
      integer                 :: ntest = 11
      integer                 :: numnumerrec_shotfile
      integer                 :: nhalosensors
!
      integer                 :: ncid, varid, shotnumber, ifirst
      integer                 :: ialloc, ier, i, debug ,l, n, jj
      integer                 :: ntimes,ncurrent,numflux,nmirnov
      integer                 :: numcurrent
      integer                 :: ntfwindings
      integer                 :: nnedl
      integer                 :: ndifl
      integer                 :: nbi1p
      integer                 :: nradp
      integer                 :: imachine !0 for NSTX, 1 for DIII
      integer, dimension(3)   :: dimlens = (/1,1,1/)
      character*4             :: type
      character*80            :: coilcur_filename, halocur_filename
      character*80            :: fluxlooploc_filename,  blooploc_filename
      character*80            :: fluxloopdata_filename, bloopdata_filename
      data ifirst/1/
!============      
      if(ifirst.eq.0) go to 1000
      ifirst = 0
 
!---------------------------------------------------------------------
      debug = 0
	if(debug .gt. 0) print *, 'here is start of spdnstx:'
	open (ntest, file='testout', status='replace', form='formatted')
!---------------------------------------------------------------------
!-- open the file
	if (numargs .lt. 1 )	then
         			filename = 'nstx_shot.dat'
	else
	filename = 'nstx_shot.dat' // '.' // trim(suffix)
				endif
      ncid = 71
      open(ncid,file=trim(filename),status='old',err=20)
      go to 50
   20	print *, 'Error in reading  nstx_shot.dat'
        return
   50 continue
!	       6        = number of following numeric records to be read
!	  123456	= shotnumber (i8)
!	     411        = ntimes = number of time points (i8)
!	      12        = ncurrent = number of PF coils + 1 (for the TF) (i8)
!	       4        = nhalosensors = number of halo sensors
!	      44        = numflux = number of flux loops
!	      43        = nmirnov = number of B loops
!	coilcur_filename        (a80)
!	halocur_filename        (a80)
!	fluxlooploc_filename    (a80)
!	blooploc_filename       (a80)
!	fluxloopdata_filename   (a80)
!	bloopdata_filename      (a80)

!---read above data:
        ntfwindings = 36
        imachine = 0
	read (ncid, '(i8)')	numnumerrec_shotfile 
	if (numnumerrec_shotfile < 1)	then
		print *, 'Error, # records =',numnumerrec_shotfile, '< 1'
					return
					endif
	read (ncid, '(i8)')	shotnumber 
	read (ncid, '(i8)')	ntimes 
	read (ncid, '(i8)')	ncurrent 
	read (ncid, '(i8)')	nhalosensors 
	read (ncid, '(i8)')	numflux 
	read (ncid, '(i8)')	nmirnov 
        if(numnumerrec_shotfile.ge.7) read (ncid,'(i8)') ntfwindings
        if(numnumerrec_shotfile.ge.8) read (ncid,'(i8)') imachine ! 0 for NSTX, 1 for DIII
	write (ntest,*)  'shotnumber   ntimes  ncurrent nhalosensors numflux  nmirnov'
	write (ntest,'(6i10)') shotnumber,ntimes,ncurrent,nhalosensors,numflux,nmirnov
	if (debug > 0)	then
 	print*,' shotnumber  ntimes ncurrent nhalosensors numflux  nmirnov ntfwindings'
 	print*, shotnumber,ntimes,ncurrent,nhalosensors,numflux,nmirnov,ntfwindings
			endif
        tffactor = 2.E-7*ntfwindings
	if (ntimes < 1 .or. ncurrent < 1 ) then
     		print *, 'Error: one of the above is < 1'
     						return
     						endif
!---------------------------------------------------------------------
!-- read the filenames
	read (ncid, '(a80)') coilcur_filename
        write (ntest,'(a80)') coilcur_filename
!
        if(nhalosensors .gt. 0) then
            read (ncid,  '(a80)') halocur_filename 
            write(ntest, '(a80)') halocur_filename 
        endif
!
        if(numflux .gt. 0) then      
            read (ncid, '(a80)') fluxlooploc_filename
            write(ntest,'(a80)') fluxlooploc_filename
        endif
!
        if(nmirnov .gt. 0) then
            read (ncid, '(a80)') blooploc_filename
	    write(ntest,'(a80)') blooploc_filename
        endif
!
        if(numflux .gt. 0) then
            read (ncid, '(a80)') fluxloopdata_filename 
	    write (ntest,'(a80)')fluxloopdata_filename
        endif
!
        if(nmirnov .gt. 0) then
            read (ncid, '(a80)') bloopdata_filename
	    write(ntest,'(a80)') bloopdata_filename
        endif
      close(ncid)
!---------------------------------------------------------------------
!						allocate arrays to hold data
      ALLOCATE (time (ntimes), STAT=ialloc )
      ALLOCATE (psupply_name (ncurrent), STAT=ialloc )
      ALLOCATE (psupply_current (ntimes,ncurrent), STAT=ialloc )
!									halo sensors
      if(nhalosensors.gt.0) then
        ALLOCATE (halosensor_name (nhalosensors), STAT=ialloc )
        ALLOCATE (halosensor_current (ntimes,nhalosensors), STAT=ialloc )
      endif
!									flux loops
      if(numflux.gt.0) then
        ALLOCATE (flux_active (numflux), STAT=ialloc )
        ALLOCATE (flux_name (numflux), STAT=ialloc )
        ALLOCATE (flux_r (numflux), STAT=ialloc )
        ALLOCATE (flux_z (numflux), STAT=ialloc )
        ALLOCATE (flux_data (ntimes,numflux), STAT=ialloc )
      endif
!									B loops
      if(nmirnov.gt.0) then
        ALLOCATE (mirnov_active (nmirnov), STAT=ialloc )
        ALLOCATE (mirnov_name (nmirnov), STAT=ialloc )
        ALLOCATE (mirnov_r (nmirnov), STAT=ialloc )
        ALLOCATE (mirnov_z (nmirnov), STAT=ialloc )
        ALLOCATE (mirnov_tangle(nmirnov), STAT=ialloc )
        ALLOCATE (mirnov_pangle(nmirnov), STAT=ialloc )
        ALLOCATE (mirnov_data (ntimes,nmirnov), STAT=ialloc )
      endif
!---------------------------------------------------------------------
!									dummy for current_data
	numcurrent = ntimes
	ALLOCATE (current_active (numcurrent), STAT=ialloc )
	ALLOCATE (current_data (ntimes,numcurrent), STAT=ialloc )
	current_data = 1.d0      
!                                              line integrated density trace - dummy for now
	nnedl = ntimes
	ALLOCATE (midplane_nedl(nnedl), STAT=ialloc )
	midplane_nedl = 1.d0
!                                              diamagnetic flux density trace - dummy for now
	ndifl = ntimes
	ALLOCATE (diamagnetic_flux(ndifl), STAT=ialloc )
	diamagnetic_flux = 1.d0

!                                              neutral beam power trace - dummy for now
	nbi1p = ntimes
	ALLOCATE (nbi_1_power(nbi1p), STAT=ialloc )
	nbi_1_power = 0.d0

!				calculate max temperature as function of time (set to 1. for now)
!
	ALLOCATE (maxte (ntimes),STAT=ialloc)
		do n=1,ntimes
		maxte(n) = 1.d0
		enddo

!---------------------------------------------------------------------
!cj#ifdef TSC_USE_NETCDF
!cj Sept 17, 2010 read coil current & flux loop from netcdf
      if(idata.eq.11) then
      call get_cc_loop_data(coilcur_filename, &
                            ntimes, time, &
                            ncurrent, psupply_name, psupply_current, &
                            numflux, flux_name, flux_data, flux_r, flux_z)
      endif
!cj#endif

!---------------------------------------------------------------------
!				Read coil current names, current and time values
!
!cj#ifndef TSC_USE_NETCDF
      if(idata.eq.9) then
	call getshot (ncid, coilcur_filename, ntimes, time, ncurrent, psupply_name, psupply_current)
      endif
!cj#endif
	if(nhalosensors.gt.0) call getshot (ncid, halocur_filename, ntimes, time, nhalosensors, halosensor_name, halosensor_current)
!cj#ifndef TSC_USE_NETCDF
      if(idata.eq.9) then
        if(imachine.eq.0) then
!.........NSTX
	  if(numflux.gt.0) call getshot (ncid, fluxloopdata_filename, ntimes, time, numflux, flux_name, flux_data)
        else
!.........DIII-D
           if(numflux.gt.0) call D3fluxvals(ncid, fluxloopdata_filename, ntimes, time, numflux, flux_name, flux_data)
        endif
      endif
!cj#endif
!
	if(nmirnov.gt.0) call getshot (ncid, bloopdata_filename,    ntimes, time, nmirnov, mirnov_name, mirnov_data)
!---------------------------------------------------------------------
!				Read flux and B loop locations and angles
        if(debug.gt.0) print *, imachine,numflux, &
                       'imachine and numflux before call to getlocang'
!
!cj#ifndef TSC_USE_NETCDF
      if(idata.eq.9) then
        if(imachine.eq.0) then
!.........NSTX
	  if(numflux.gt.0) call getloc (ncid, fluxlooploc_filename, numflux, flux_r, flux_z)
        else
!.........DIII-D
          if(numflux.gt.0) call d3fluxloops( ncid, fluxlooploc_filename,numflux, flux_r, flux_z)
        endif
      endif
!cj#endif
	if(nmirnov.gt.0) call getlocang (ncid, blooploc_filename,    nmirnov, mirnov_r, mirnov_z, mirnov_tangle, mirnov_pangle)
        if(debug.gt.0) print *, 'after last call to getlocang:'

!cj 2010-10-06 compare netcdf ascii data
!cj      write(93,9971) ((psupply_current(n,l),n=1,ntimes), l=1,ncurrent)
!cj      write(94,9971) ((flux_data(n,l),n=1,ntimes), l=1,numflux)
!cj 9971 format(1p6e12.4)

!cj      print *
!cj      print *, "from parent call ..."
!cj      print *, time(1), time(ntimes)
!cj      print *, trim(psupply_name(1)), '  ', trim(psupply_name(2))
!cj      print *, trim(psupply_name(ncurrent-1)), '  ', trim(psupply_name(ncurrent))
!cj      if(numflux.gt.0) print *, trim(flux_name(1)), '  ', trim(flux_name(numflux)) 
!cj      print *, "psupply_current ntimes=1"
!cj      write(*,'(4E12.4)') psupply_current(1,1), psupply_current(1,2), &
!cj                          psupply_current(1,ncurrent-1), psupply_current(1,ncurrent)
!cj      print *, "psupply_current ntimes"
!cj      write(*,'(4E12.4)') psupply_current(ntimes,1), psupply_current(ntimes,2), &
!cj                          psupply_current(ntimes,ncurrent-1), psupply_current(ntimes,ncurrent)
!cj      if(numflux.gt.0) print *, flux_data(1,1), flux_data(1,numflux),  &
!cj               flux_data(ntimes,1), flux_data(ntimes,numflux)
!cj      if(numflux.gt.0) print *, flux_r(1), flux_r(numflux), flux_z(1), flux_z(numflux)

!!!!!!!	stop
!
!--------------------------------------------------------------------- 
!
!      ALLOCATE (maxte (dimlens(1)),STAT=ialloc)
!      do n=1,dimlens(1)
!        amaxtemp = 0._R8
!          do l=1,dimlens(2)
!          amaxtemp = max(amaxtemp,mpts_Te(n,l))
!          end do
!        maxte(n) = amaxtemp
!      end do
!      if(debug.eq.1) then
!      print*,'maxTE data'
!      do n=1,dimlens(1)
!      print*,maxte(n)
!      enddo
!        endif
 
!---------------------------------------------------------------------
 1000 continue
!
	if (debug > 0)	then
	print *, 'here is call to spdnstxo:'
	print *, 'shotnumber, ntimes, ncurrent, numflux, nmirnov, numcurrent, nhalosensors, nnedl, ndifl, nbi1p'
	print *, shotnumber, ntimes, ncurrent, numflux, nmirnov, numcurrent, nhalosensors, nnedl, ndifl, nbi1p
			endif
!
!				note added quantities: halosensor_current, nhalosensors
!.....
      call spdnstxo(shotnumber,time,ntimes,                              &  
     &     psupply_current, ncurrent,                                    &  
     &     flux_active,flux_r,flux_z,flux_data,numflux,                  &  
     &     mirnov_active,mirnov_r,mirnov_z,mirnov_tangle,mirnov_pangle,  &  
     &                   mirnov_data,nmirnov,                            &  
     &     current_active,current_data,numcurrent ,                      &  
     &     halosensor_current, nhalosensors,                             &  
     &     midplane_nedl, nnedl, diamagnetic_flux, ndifl, nbi_1_power,   &
     &     nbi1p, maxte, tffactor,imachine)
!
!
!	print *, 'here is call to testxo:'
!      call testxo (shotnumber, time, ntimes, psupply_current, ncurrent,   &
!     &             flux_active, flux_r, flux_z, flux_data, numflux,       &
!     &     mirnov_active,mirnov_r,mirnov_z,mirnov_tangle,mirnov_pangle,   &  
!     &                   mirnov_data,nmirnov,                             &
!     &    current_active, current_data, numcurrent,                       &  
!     &     halosensor_current, nhalosensors,                              &  
!     &     midplane_nedl, nnedl, diamagnetic_flux, ndifl, nbi_1_power,    &
!     &     nbi1p, maxte)
!      return
 
 
      end
!---------------------------------------------------------------------
!---------------------------------------------------------------------
	subroutine testxo (shotnumber, stime, ntimes, si, ncurrent,	&
     &			   iflux_active, rfl, zfl, sfl, numfluxloops,   &
     &     mirnov_active,bmirnov_r,bmirnov_z,bmirnov_tangle,            &  
     &     bmirnov_pangle, bmirnov_data,nmirnov,                        &
     &     icurrent_active, scurrent, numcurrent,                       &  
     &     halosensor_current, nhalosensors,                            &  
     &     expnedl, nnedl, expdifl, ndifl, nbipower, nbi1p, maxte)
!
      USE CLINAM
      USE POLCUR
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

	integer  shotnumber, ntimes, ncurrent, numfluxloops, nmirnov
	integer  mirnov_active (nmirnov), iflux_active (numfluxloops)
	integer  numcurrent, nnedl, ndifl, nbi1p, nhalosensors
	integer  icurrent_active (numcurrent)
!
	real*8  stime (ntimes), si (ncurrent)
	real*8  rfl (numfluxloops), zfl (numfluxloops)
	real*8  sfl (ntimes, numfluxloops)
	real*8  bmirnov_r (nmirnov), bmirnov_z(nmirnov)
	real*8  bmirnov_tangle (nmirnov), bmirnov_pangle(nmirnov)
	real*8  bmirnov_data (ntimes, nmirnov)
	real*8  scurrent (ntimes, numcurrent)
	real*8  halosensor_current (ntimes, nhalosensors)
	real*8  expnedl (nnedl), expdifl (ndifl), nbipower(nbi1p), maxte (ntimes)

!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pprsav,rnormsv,fbchisv,pcursv,vlpsav,xmagzsv,zmagzsv
      REAL*8 tevv0sv,ffac0sv,zeffvsv,beampsv,alphrsv,betarsv
      REAL*8 frparsv,thalosv,whalosv,heactsv
      REAL*8 azervsv, ezervsv, dzervsv, rzervsv
      REAL*8 vescure,totcure,vescure2,chicur,rline,difle,rdis,pval
      REAL*8 rlinae,abeamp,fac,pone,r1,z1,graddum,dpxarg,dpzarg
      REAL*8 gsdum,psarg,pxzarg,pxxarg,pzzarg,p1,csum,atot,anumfl
!============
!
      dimension pprsav(ptpts),rnormsv(ptpts), fbchisv(ptpts),         &
     & pcursv(ptpts),vlpsav(ptpts),                                    &
     & xmagzsv(ptpts),zmagzsv(ptpts),tevv0sv(ptpts),                  &
     & ffac0sv(ptpts),zeffvsv(ptpts),beampsv(ptpts),                  &
     & alphrsv(ptpts),betarsv(ptpts),frparsv(ptpts),                  &
     & thalosv(ptpts), whalosv(ptpts), heactsv(ptpts),                &
     & azervsv(ptpts), ezervsv(ptpts), dzervsv(ptpts), rzervsv(ptpts)
!
!
!
!============
!From the coil currents data file:
! tinterp   Ip1 IOH  IPF1AU  IPF2U IPF3U IPF4 IPF5 IPF3L  IPF2L IPF1b  IPF1AL  ITF
!            1   2     3       4     5    6     7    8      9     10     11     12
!
	integer :: icnstxplas = 1, icnstxoh   = 2, icnstxpf1au = 3, icnstxpf2u = 4
	integer :: icnstxpf3u = 5, icnstxpf4  = 6, icnstxpf5   = 7, icnstxpf3l = 8
	integer :: icnstxpf2l = 9, icnstxpf1b =10, icnstxpf1al =11, icnstxtf  = 12
!
	REAL*8 totcureo
	integer                 :: ialloc
      INTEGER icount,l,i,ios22,n,indx,iz,ip,lz,lp,k,ii

	print *, 'here is start of testxo:'
!print *, 'shotnumber, ntimes, ncurrent, numfluxloops, nmirnov, numcurrent, nhalosensors, nnedl, ndifl, nbi1p'
	print *, shotnumber, ntimes, ncurrent, numfluxloops
	print *, 'nmirnov, numcurrent, nhalosensors, nnedl, ndifl, nbi1p'

	ALLOCATE (ifrst (50), STAT=ialloc )
	print *, 'ifrst(9) = ', ifrst(9)
!							Initialize
!
!...........................................................
!
!...save trajectories from input file
!
!...........................................................
      do 5 l=1,ntpts
      pprsav(l) = ppres(l)
      vlpsav(l) = vloopv(l)
      rnormsv(l) = rnorm(l)
      fbchisv(l) = fbchia(l)
      xmagzsv(l) = xmagz(l)
      zmagzsv(l) = zmagz(l)
      tevv0sv(l) = tevv0(l)
      ffac0sv(l) = ffac0(l)
      zeffvsv(l) = zeffv(l)
      frparsv(l) = frcparv(l)
      betarsv(l) = betarv(l)
      alphrsv(l) = alpharv(l)
      thalosv(l) = thalov(l)
      whalosv(l) = whalov(l)
      heactsv(l) = heactv(l)
      azervsv(l) = azerv(l)
      ezervsv(l) = ezerv(l)
      dzervsv(l) = dzerv(l)
      rzervsv(l) = rzerv(l)
    5 continue
!
!.....open the file to write the simulation and experimental data
      if( numargs .lt. 1 ) then
         filename = 'tscnstxdata_out'
      else
         filename = 'tscnstxdata_out' // '.' // trim(suffix)
      end if
      open (no167a,file=trim(filename),status='unknown',iostat=ios22)
!
!
!.....end of initialization section
!
	return
	end
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      subroutine getshot (ncid, filename, ntimes, time, numcols, column_name, column_value)
!
!			Read NSTX shot data files
!
	integer :: debug = 0
	integer :: ntest = 11
	integer :: ncid, ntimes, numcols
	integer :: numcolsout
	character*80  filename
	character*50  column_name(numcols)    
	real*8        time (ntimes)  
	real*8        column_value(ntimes, numcols)  
!
	column_value = 0.d0					! initialize array
!
      ncid = 71
	open (ncid, file=filename, status='old', form='formatted',err=100)
	read (ncid, '(13x,80a13)')  column_name
		do i=1,numcols
		column_name(i) = trim (adjustl (column_name(i)))
		if (debug > 0)  print *, column_name(i)
		enddo
	write (ntest,'(12x,80a9)') (column_name (i), i=1,numcols)
	numcolsout = numcols
	if (numcolsout > 22)  numcolsout = 22
	do i=1,ntimes
	      read (ncid, '(80f13.1)')  time(i), (column_value (i,j), j=1,numcols)
!>>>>>debug
!     read(ncid,'(80f13.1)')  time(i), column_value(i,1),(column_value (i,j), j=numcols-2,numcols-1), &
!                             (column_value (i,j), j=2,numcols-3),column_value(i,numcols)
	write (ntest,'(f8.3,1x,1p40e9.1)')  time(i), (column_value (i,j), j=1,numcolsout)
!	print *, time(i), (column_value (i,j), j=1,5)
	enddo		
	print *, 'getshot: data read from ', trim(filename)
	print *, 'getshot: ntimes, numcols, tmin, tmax =', ntimes, numcols, time(1), time(ntimes)
	print *, 'getshot: ', column_value(1,1), column_value(1,numcols),  &
                              column_value(ntimes,1), column_value(ntimes,numcols)
      close(ncid)
	return
  100	print *, 'Error opening file:', filename
	end
!---------------------------------------------------------------------
      subroutine D3fluxvals(ncid, filename, ntimes, time, numcols, column_name, column_value)
!
!			Read DIII shot data files for fluxes
!
	integer :: debug = 0
	integer :: ntest = 11
	integer :: ncid, ntimes, numcols
	integer :: numcolsout
	character*80  filename
	character*50  column_name(numcols)    
	real*8        time (ntimes)  
	real*8        column_value(ntimes, numcols)  
!
	column_value = 0.d0					! initialize array
!
      ncid = 71
      debug = 1
	open (ncid, file=filename, status='old', form='formatted',err=100)
	read (ncid, '(13x,80a13)')  column_name
		do i=1,numcols
		column_name(i) = trim (adjustl (column_name(i)))
		if (debug > 0)  print *, column_name(i)
		enddo
	write (ntest,'(12x,80a9)') (column_name (i), i=1,numcols)
	numcolsout = numcols
	if (numcolsout > 22)  numcolsout = 22
	do i=1,ntimes
	      read (ncid, '(80f13.1)')  time(i), (column_value (i,j), j=1,numcols)
	write (ntest,'(f8.3,1x,1p40e9.1)')  time(i), (column_value (i,j), j=1,numcolsout)
!	print *, time(i), (column_value (i,j), j=1,5)
	enddo		
	print *, 'D3fluxvals: data read from : ', filename
	print *, 'D3fluxvals: ntimes, numcols, tmin, tmax =', ntimes, numcols, time(1), time(ntimes)
      close(ncid)
	return
  100	print *, 'Error opening file:', filename
	end
!---------------------------------------------------------------------
      subroutine getlocang (ncid, filename, numloops, rloop, zloop, tangle, pangle)
!
!			Read NSTX shot flux & B loop location files
!
	integer :: ntest = 11
	integer :: ncid, ntimes, numloops
	character*80  filename
	character*80  column_name 
	character*15  loop_name 
	real*8        rloop(numloops), zloop(numloops), tangle(numloops), pangle(numloops)  
!
!
	rloop  = 0.d0 ; 	zloop  = 0.d0				! initialize arrays
	tangle = 0.d0 ; 	pangle = 0.d0				! initialize arrays
! 
        ncid = 71
	open (ncid, file=filename, status='old', form='formatted',err=100)
	read (ncid, '(a80)')  column_name
	write (ntest,'(a80)') column_name
        do i=1,numloops
	  read (ncid, '(a15, 8f15.1)',end=200)  loop_name, pangle(i), rloop(i), zloop(i), tangle(i)
	  write (ntest,'(a15,8f15.3)' )  loop_name, pangle(i), rloop(i), zloop(i), tangle(i)
        enddo
        close(ncid)
!		
	print *, 'getlocang: loop location & angle data read from : ', filename
	print *, 'getlocang: numloops =', numloops
	return
  100	print *, 'Error opening file:', filename
  	return
  200	print *, 'End of file:', filename
  	return
	end
      subroutine getloc (ncid, filename, numloops, rloop, zloop)
!
!			Read NSTX shot flux loop location files
!
	integer :: ntest = 11
	integer :: ncid, ntimes, numloops
	character*80  filename
	character*80  column_name 
	character*15  loop_name 
	real*8        rloop(numloops), zloop(numloops), pangle  
!
!
	rloop  = 0.d0 ; 	zloop  = 0.d0				! initialize array
! 
        ncid = 71
	open (ncid, file=filename, status='old', form='formatted',err=100)
	read (ncid, '(a80)')  column_name
	write (ntest,'(a80)') column_name
        do i=1,numloops
	  read (ncid, '(a15, 8f15.1)',end=200)  loop_name, pangle, rloop(i), zloop(i)
	  write (ntest,'(a15,8f15.3)' )  loop_name, pangle, rloop(i), zloop(i) 
        enddo
        close(ncid)
!		
	print *, 'getloc: loop location & angle data read from : ', filename
	print *, 'getloc: numloops =', numloops
	return
  100	print *, 'Error opening file:', filename
  	return
  200	print *, 'End of file:', filename
  	return
	end
        subroutine d3fluxloops( ncid, filename,numloops,rloop,zloop)
!
!			Read DIII shot flux loop location files
!
	integer :: ntest = 11
	integer :: ncid, numloops
	character*80  filename
	character*80  column_name 
	character*15  loop_name 
	real*8        rloop(numloops), zloop(numloops), tangle(numloops), pangle(numloops)  
!
!
	rloop  = 0.d0 ; 	zloop  = 0.d0				! initialize arrays
	tangle = 0.d0 ; 	pangle = 0.d0				! initialize arrays
!
      ncid = 71
	open (ncid, file=filename, status='old', form='formatted',err=100)
	read (ncid, '(a80)')  column_name
	write (ntest,'(a80)') column_name
		do i=1,numloops
	read (ncid, '(a10,2f9.0)',end=200)  loop_name, rloop(i), zloop(i)
	write (ntest,'(a10,2f9.4)' )  loop_name, rloop(i), zloop(i)
		enddo		
	print *, 'd3fluxloops: loop location & angle data read from : ', filename
	print *, 'd3fluxloops: numloops =', numloops
      close(ncid)
	return
  100	print *, 'Error opening file:', filename
  	return
      stop 10
  200	print *, 'End of file:', filename
  	return
      stop 11
	end


      subroutine get_cc_loop_data(coilcur_filename, &
                            ntimes, time, &
                            ncurrent, psupply_name, psupply_current, &
                            numflux, flux_name, flux_data, flux_r, flux_z)
!			Jin Chen
!			17 Sept 2010 
!
      use netcdf
      USE CLINAM, ONLY : ineg

!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
      REAL(KIND=r8) time(ntimes)
 
      character*13 psupply_name(ncurrent)
      REAL(KIND=r8) psupply_current(ntimes,ncurrent)
      integer, dimension(:), allocatable :: psupply_index
      character*13,   dimension(:),   allocatable :: dummy_cc_name
      REAL(KIND=r8), dimension(:,:), allocatable :: dummy_cc
      REAL(KIND=r8), dimension(:), allocatable :: dummy_ip, dummy_itf
 
      character*13 flux_name(numflux)
      REAL(KIND=r8) flux_r(numflux)
      REAL(KIND=r8) flux_z(numflux)
      REAL(KIND=r8) flux_data(ntimes,numflux)
      integer loop_index(numflux)
      REAL(KIND=r8), dimension(:,:), allocatable :: dummy_flux
 
      integer                 :: ncid, varid, shotnumber, ialloc, ier, i, j
      integer                 :: ntimes,ncurrent,numflux
      character*80            :: coilcur_filename
      integer                 :: ntim,ncur,nflux,ntim_loop

      character*20            :: tDimName, tDimName_loop
      integer                 :: tDimId, tVarId, tDimId_loop

      character*20            :: ccDimName
      integer                 :: ccDimId, iccVarId, ccVarId, ccNameId
      integer                 :: ipVarId, itfVarId 

      character*20            :: loopDimName
      integer                 :: loopDimId, iloopVarId, loopVarId, rloopVarId, zloopVarId, loopNameId

      integer                 :: pfip_name_len
      character*50            :: pfip_name
      integer                 :: ohbt_name_len
      character*50            :: ohbt_name
      integer                 :: loop_name_len
      character*50            :: loop_name
      integer                 :: pos_name_len
      character*50            :: pos_name
      integer                 :: shot_len
      integer                 :: shot
      integer                 :: date_len
      character*50            :: date
      integer                 :: sdate_len
      character*50            :: sdate
      integer                 :: machine_len
      character*50            :: machine

!============      

!coilcur_filename = "/p/swim/jchen/IPS/ips/components/epa/tsc/NETCDF_sept08/out.cdf"
      ncid = 71

      ier = nf90_open(trim(coilcur_filename), nf90_nowrite, ncid)
      if (ier /= nf90_noerr) then
         ineg=37
         write(*,*)  "open error", trim(nf90_inq_libvers())
         return
      else 
         write(*,*)  "open ok", trim(coilcur_filename), trim(nf90_inq_libvers())
      endif
 37   format (1x,a,'(',i4,')',a,i4,'(',')') 

!time
      ier = nf90_inq_dimid(ncid, "time_cc", tDimId)
      if (ier /= nf90_noerr) write(*,*)  "inquire tdimid error"
      ier = nf90_inquire_dimension(ncid, tDimId, tDimName, ntim)
      if (ier /= nf90_noerr) write(*,*)  "inquire tdim error"
      if(ntim .ne. ntimes) then
         ineg=37
         write(*,37) "ntimes: nstx_shot.dat",ntimes,                    &
                     "not match netcdf",ntim
         return
      else
         print *, tDimId, trim(tDimName), ntim, ntimes
      endif

      ier = nf90_inq_dimid(ncid, "time_loop", tDimId_loop)
      if (ier /= nf90_noerr) write(*,*)  "inquire tdimid_loop error"
      ier = nf90_inquire_dimension(ncid, tDimId_loop, tDimName_loop, ntim_loop)
      if (ier /= nf90_noerr) write(*,*)  "inquire tdim_loop error"
      if(ntim_loop .ne. ntimes) then
         ineg=37
         write(*,37) "ntimes: nstx_shot.dat",ntimes,                    &
                     "not match netcdf",ntim_loop
         return
      else
         print *, tDimId_loop, trim(tDimName_loop), ntim_loop, ntimes
      endif

!coil current
      ier = nf90_inq_dimid(ncid, "index_cc", ccDimId)
      if (ier /= nf90_noerr) write(*,*)  "inquire ccdimid error"
      ier = nf90_inquire_dimension(ncid, ccDimId, ccDimName, ncur)
      if (ier /= nf90_noerr) write(*,*)  "inquire ccdim error"
      if(ncur .ne. (ncurrent-2) ) then
         ineg=37
         write(*,37) "ncurrent: nstx_shot.dat",ncurrent,                &
                     "not match netcdf ncur",ncur
         write(*,*)  "ncurrent: nstx_shot.dat ncurrent ?= ncur + 2"
         stop
      else
         print *, ccDimId, trim(ccDimName), ncur, ncurrent
      endif

!flux loop
      if(numflux.gt.0) then
      ier = nf90_inq_dimid(ncid, "index_loop", loopDimId)
      if (ier /= nf90_noerr) write(*,*)  "inquire loopdimid error"
      ier = nf90_inquire_dimension(ncid, loopDimId, loopDimName, nflux)
      if (ier /= nf90_noerr) write(*,*)  "inquire loopdim error"
      if(nflux .ne. numflux) then
         ineg=37
         write(*,*) "numflux: nstx_shot.dat doesn't match netcdf"
         write(*,37) "numflux: nstx_shot.dat",numflux,                  &
                     "not match netcdf nflux",nflux
         return
      else
         print *, loopDimId, trim(loopDimName), nflux, numflux
      endif
      endif

      ALLOCATE (psupply_index (ncur), STAT=ialloc )
      ALLOCATE (dummy_cc_name (ncur), STAT=ialloc )
      ALLOCATE (dummy_cc (ncur,ntim), STAT=ialloc )
      ALLOCATE (dummy_ip (ntim), STAT=ialloc )
      ALLOCATE (dummy_itf (ntim), STAT=ialloc )
      if(nflux.gt.0) then
        ALLOCATE (dummy_flux (nflux,ntim), STAT=ialloc )
      endif

!time
      ier = nf90_inq_varid(ncid, "time_cc", tVarId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire tvarid error"
      ier = nf90_get_var(ncid, tVarId, time, start=(/1/), count=(/ntim/))
      if(ier /= nf90_NoErr) write(*,*)  "get tvar error"
      print *, time(1), time(ntim)
!  -1.24600000000000        8.16500000000000     

!coil current
      ier = nf90_inq_varid(ncid, "index_cc", iccVarId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire iccvarid error"
      ier = nf90_get_var(ncid, iccVarId, psupply_index, start=(/1/), count=(/ncur/))
      if(ier /= nf90_NoErr) write(*,*)  "get iccvar error"
      print *, psupply_index(1), psupply_index(ncur)
!           1          20 

      ier = nf90_inq_varid(ncid, "name_cc", ccNameId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire ccnameid error"
      ier = nf90_get_var(ncid, ccNameId, psupply_name, start=(/1,1/), count=(/13,ncur/))
      if(ier /= nf90_NoErr) write(*,*)  "get ccnamevar error", trim(nf90_strerror(ier)) 
         do j=1, ncur
            psupply_name(j+1) = psupply_name(j)
         enddo
         print *, trim(psupply_name(1+1)), trim(psupply_name(ncur+1))
!f1aecoilb

      ier = nf90_inq_varid(ncid, "cc", ccVarId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire ccvarid error"
      ier = nf90_get_var(ncid, ccVarId, dummy_cc, start=(/1,1/), count=(/ncur,ntim/))
      if(ier /= nf90_NoErr) write(*,*)  "get ccvar error"
      print *, dummy_cc(2,2), dummy_cc(ncur-1,2),  &
               dummy_cc(2,ntim), dummy_cc(ncur-1,ntim)
!  0.000000000000000E+000   184.501760000001       -44.3896600000000     807.195400000000     

      ier = nf90_inq_varid(ncid, "ip", ipVarId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire ipvarid error"
      ier = nf90_get_var(ncid, ipVarId, dummy_ip, start=(/1/), count=(/ntim/))
      if(ier /= nf90_NoErr) write(*,*)  "get ipvar error"
      psupply_name(1) = 'ip'
      print *, dummy_ip(1), dummy_ip(ntim)

      ier = nf90_inq_varid(ncid, "itf", itfVarId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire itfvarid error"
      ier = nf90_get_var(ncid, itfVarId, dummy_itf, start=(/1/), count=(/ntim/))
      if(ier /= nf90_NoErr) write(*,*)  "get itfvar error"
      psupply_name(ncur+2) = 'itf'
      print *, dummy_itf(1), dummy_itf(ntim)

      if(nflux.gt.0) then
!flux loop
      ier = nf90_inq_varid(ncid, "index_loop", iloopVarId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire iloopvarid error"
      ier = nf90_get_var(ncid, iloopVarId, loop_index, start=(/1/), count=(/nflux/))
      if(ier /= nf90_NoErr) write(*,*)  "get iloopvar error"
      print *, loop_index(1), loop_index(nflux)
!           1          44 

      ier = nf90_inq_varid(ncid, "name_loop", loopNameId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire loopnameid error"
      ier = nf90_get_var(ncid, loopNameId, flux_name, start=(/1,1/), count=(/13,nflux/))
      if(ier /= nf90_NoErr) write(*,*)  "get loopnamevar error", trim(nf90_strerror(ier)) 
      print *, trim(flux_name(1)), trim(flux_name(nflux))
!PSF1APSI3L

      if(ntim_loop .eq. ntim) then
      ier = nf90_inq_varid(ncid, "flux_loop", loopVarId)
      else
      ier = nf90_inq_varid(ncid, "flux_intloop", loopVarId)
      endif
      if(ier /= nf90_NoErr) write(*,*)  "inquire loopvarid error"
      ier = nf90_get_var(ncid, loopVarId, dummy_flux, start=(/1,1/), count=(/nflux,ntim/))
      if(ier /= nf90_NoErr) write(*,*)  "get loopvar error"
      print *, dummy_flux(1,1), dummy_flux(nflux,1),  &
               dummy_flux(1,ntim), dummy_flux(nflux,ntim)
!  0.000000000000000E+000  0.000000000000000E+000 -3.557545804110389E-002 0.000000000000000E+000

      ier = nf90_inq_varid(ncid, "r_loop", rloopVarId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire rloopvarid error"
      ier = nf90_get_var(ncid, rloopVarId, flux_r, start=(/1/), count=(/nflux/))
      if(ier /= nf90_NoErr) write(*,*)  "get rloopvar error"
      ier = nf90_inq_varid(ncid, "z_loop", zloopVarId)
      if(ier /= nf90_NoErr) write(*,*)  "inquire zloopvarid error"
      ier = nf90_get_var(ncid, zloopVarId, flux_z, start=(/1/), count=(/nflux/))
      if(ier /= nf90_NoErr) write(*,*)  "get zloopvar error"
      print *, flux_r(1), flux_r(nflux), flux_z(1), flux_z(nflux)
!  0.892900000000000        1.77590000000000       0.168300000000000     -1.30530000000000     
      endif

!global attribute
      ier = nf90_inquire_attribute(ncid,nf90_global,"pfip_name", len = pfip_name_len)
      if(ier /= nf90_NoErr) write(*,*)  "inquire global att error"
      ier = nf90_get_att(ncid, nf90_global, "pfip_name", pfip_name)
      print *, 'ga pfip_name:', pfip_name_len, pfip_name
      ier = nf90_inquire_attribute(ncid,nf90_global,"ohbt_name", len = ohbt_name_len)
      if(ier /= nf90_NoErr) write(*,*)  "inquire global att error"
      ier = nf90_get_att(ncid, nf90_global, "ohbt_name", ohbt_name)
      print *, 'ga ohbt_name:', ohbt_name_len, ohbt_name
      ier = nf90_inquire_attribute(ncid,nf90_global,"loop_name", len = loop_name_len)
      if(ier /= nf90_NoErr) write(*,*)  "inquire global att error"
      ier = nf90_get_att(ncid, nf90_global, "loop_name", loop_name)
      print *, 'ga loop_name:', loop_name_len, loop_name
      ier = nf90_inquire_attribute(ncid,nf90_global,"pos_name", len = pos_name_len)
      if(ier /= nf90_NoErr) write(*,*)  "inquire global att error"
      ier = nf90_get_att(ncid, nf90_global, "pos_name", pos_name)
      print *, 'ga pos_name:', pos_name_len, pos_name
      ier = nf90_inquire_attribute(ncid,nf90_global,"shot", len = shot_len)
      if(ier /= nf90_NoErr) write(*,*)  "inquire global att error"
      ier = nf90_get_att(ncid, nf90_global, "shot", shot)
      print *, 'ga shot:', shot_len, shot
      ier = nf90_inquire_attribute(ncid,nf90_global,"date", len = date_len)
      if(ier /= nf90_NoErr) write(*,*)  "inquire global att error"
      ier = nf90_get_att(ncid, nf90_global, "date", date)
      print *, 'ga date:', date_len, date
      ier = nf90_inquire_attribute(ncid,nf90_global,"sdate", len = sdate_len)
      if(ier /= nf90_NoErr) write(*,*)  "inquire global att error"
      ier = nf90_get_att(ncid, nf90_global, "sdate", sdate)
      print *, 'ga sdate:', sdate_len, sdate
      ier = nf90_inquire_attribute(ncid,nf90_global,"machine", len = machine_len)
      if(ier /= nf90_NoErr) write(*,*)  "inquire global att error"
      ier = nf90_get_att(ncid, nf90_global, "machine", machine)
      print *, 'ga machine:', machine_len, machine

      ier = nf90_close(ncid)
      if (ier /= nf90_noerr) then
         write(*,*)  "close error" 
      else
         write(*,*)  "OK close ", trim(coilcur_filename)
      endif

!swap data: flip the index
      do i=1, ntim
         psupply_current(i,1) = dummy_ip(i)
         do j=1, ncur
            psupply_current(i,j+1) = dummy_cc(j,i)
         enddo
         psupply_current(i,ncur+2) = dummy_itf(i)

         if(nflux.gt.0) then
         do j=1, nflux
         flux_data(i,j) = dummy_flux(j,i)
         enddo
         endif
      enddo
         i=1
         write(*,'(4E12.4)') psupply_current(i,1), psupply_current(i,2), &
                             psupply_current(i,ncur+1), psupply_current(i,ncur+2)
         i=ntim
         write(*,'(4E12.4)') psupply_current(i,1), psupply_current(i,2), &
                             psupply_current(i,ncur+1), psupply_current(i,ncur+2)

      DEALLOCATE ( psupply_index, STAT=ialloc)
      DEALLOCATE ( dummy_cc_name, STAT=ialloc)
      DEALLOCATE ( dummy_cc, STAT=ialloc)
      DEALLOCATE ( dummy_ip, STAT=ialloc)
      DEALLOCATE ( dummy_itf, STAT=ialloc)
      DEALLOCATE ( dummy_flux, STAT=ialloc)

      return
      end
