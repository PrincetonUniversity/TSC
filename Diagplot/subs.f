      subroutine readmag(icolor,flux_r,flux_z,flux_active,nloopsp,nc,nf)
!			Rewrite to use new NSTX shot data files.
!			ROS.  20 Nov 2009 
!
!                       Generalize to also work for DIII-D
!                       SCJ   Feb 2010
!!!!!!!      USE ezcdf
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
      integer,        dimension(nloopsp)   :: flux_active, icolor
      REAL(KIND=r8), dimension(nloopsp)  :: flux_r
      REAL(KIND=r8), dimension(nloopsp)  :: flux_z
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
       REAL(KIND=r8)                              :: amaxtemp, tffactor
 
      integer  :: icolor
      integer            :: ntest = 11
      integer            :: numnumerrec_shotfile
      integer            :: nhalosensors, nf
!
      integer            :: ncid, varid, shotnumber, ifirst
      integer            :: ialloc, ier, i, debug ,l, n, jj
      integer            :: ntimes,nc,nloopsp,nmirnov
      integer            :: numcurrent
      integer            :: ntfwindings
      integer            :: nradp
      integer            :: imachine !0 for NSTX, 1 for DIII
      integer, dimension(3)   :: dimlens = (/1,1,1/)
      character*4        :: type
      character*80       :: coilcur_filename, halocur_filename,filename
      character*80       :: fluxlooploc_filename,  blooploc_filename
      character*80       :: fluxloopdata_filename, bloopdata_filename
      data ifirst/1/
!============      
      if(ifirst.eq.0) go to 1000
      ifirst = 0
      do i=1,nloopsp
        flux_active(i) = 1
      enddo
 
!---------------------------------------------------------------------
      debug = 0
	if(debug .gt. 0) print *, 'here is start of spdnstx:'
	open (ntest, file='testout', status='replace', form='formatted')
!---------------------------------------------------------------------
!-- open the file
         			filename = 'nstx_shot.dat'
      ncid = 71
      open(ncid,file=trim(filename),status='old',err=20)
      go to 50
   20	print *, '>>> Cannot find file nstx_shot.dat <<<'
        return
   50 continue
!	       6        = number of following numeric records to be read
!	  123456	= shotnumber (i8)
!	     411        = ntimes = number of time points (i8)
!	      12        = nc = number of PF coils + 2 (for the IP and TF) (i8)
!	       4        = nhalosensors = number of halo sensors
!	      44        = nf = number of flux loops
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
	read (ncid, '(i8)')	nc
	read (ncid, '(i8)')	nhalosensors 
	read (ncid, '(i8)')	nf
	read (ncid, '(i8)')	nmirnov 
        if(numnumerrec_shotfile.ge.7) read (ncid,'(i8)') ntfwindings
        if(numnumerrec_shotfile.ge.8) read (ncid,'(i8)') imachine ! 0 for NSTX, 1 for DIII
	write (ntest,*) 'shotnumber   ntimes  nc nhalosensors ',
     &  '    nf  nmirnov'
	write (ntest,'(6i10)') shotnumber,ntimes,nc,nhalosensors,
     &       nf,nmirnov
	if (debug > 0)	then
 	print*,' shotnumber  ntimes nc nhalosensors nf  ',
     &     ' nmirnov ntfwindings'
 	print*, shotnumber,ntimes,nc,nhalosensors,nf,
     &       nmirnov,ntfwindings
			endif
        tffactor = 2.E-7*ntfwindings
	if (ntimes < 1 .or. nc < 1 ) then
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
        if(nf .gt. 0) then
            read (ncid, '(a80)') fluxlooploc_filename
            write(ntest,'(a80)') fluxlooploc_filename
        endif
!
        if(nmirnov .gt. 0) then
            read (ncid, '(a80)') blooploc_filename
	    write(ntest,'(a80)') blooploc_filename
        endif
!
        if(nf .gt. 0) then
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
!	allocate arrays to hold data
      ALLOCATE (time (ntimes), STAT=ialloc )
      ALLOCATE (psupply_name (nc), STAT=ialloc )
      ALLOCATE (psupply_current (ntimes,nc), STAT=ialloc )
!									halo sensors
      if(nhalosensors.gt.0) then
        ALLOCATE (halosensor_name (nhalosensors), STAT=ialloc )
        ALLOCATE (halosensor_current(ntimes,nhalosensors),STAT=ialloc)
      endif
!									flux loops
      if(nf.gt.0) then
        ALLOCATE (flux_name (nf), STAT=ialloc )
        ALLOCATE (flux_data (ntimes,nf), STAT=ialloc )
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
	numcurrent = 1
!                                              line integrated density trace - dummy for now
	midplane_nedl(1) = 1.d0
!                                              diamagnetic flux density trace - dummy for now
	diamagnetic_flux(1) = 1.d0

!                                              neutral beam power trace - dummy for now

!				calculate max temperature as function of time (set to 1. for now)
!
		do n=1,ntimes
		enddo
!---------------------------------------------------------------------
!				Read coil current names, current and time values
!
	call getshot (ncid, coilcur_filename, ntimes, time, nc,
     &       psupply_name, psupply_current)
	if(nhalosensors.gt.0) call getshot (ncid, halocur_filename, 
     &       ntimes, time, nhalosensors, halosensor_name, 
     &  halosensor_current)
        if(imachine.eq.0) then
!.........NSTX
	  if(nf.gt.0) call getshot (ncid, fluxloopdata_filename,
     &    ntimes, time, nf, flux_name, flux_data)
        else
!.........DIII-D
           if(nf.gt.0) call D3fluxvals(ncid,fluxloopdata_filename,
     &            ntimes, time, nf, flux_name, flux_data)
        endif
!
	if(nmirnov.gt.0) call getshot (ncid, bloopdata_filename,    
     &           ntimes, time, nmirnov, mirnov_name, mirnov_data)
!---------------------------------------------------------------------
!				Read flux and B loop locations and angles
!
        if(imachine.eq.0) then
!.........NSTX
	  if(nf.gt.0) call getloc (ncid, fluxlooploc_filename,
     &     nf, flux_r, flux_z)
        else
!.........DIII-D
          if(nf.gt.0) call d3fluxloops( ncid,fluxlooploc_filename,
     &            nf, flux_r, flux_z)
        endif
	if(nmirnov.gt.0) call getlocang (ncid, blooploc_filename,    
     &     nmirnov, mirnov_r, mirnov_z, mirnov_tangle, mirnov_pangle)
        if(debug.gt.0) print *, 'after last call to getlocang:'
!!!!!!!	stop
!
!--------------------------------------------------------------------- 
!
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
	print *, 'shotnumber, ntimes, nc, nf, nmirnov,
     &       numcurrent, nhalosensors'
	print *, shotnumber, ntimes, nc, nf, nmirnov,
     &        numcurrent, nhalosensors
			endif
!
!				note added quantities: halosensor_current, nhalosensors
!.....
!
!
!      return
 
 
      end
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      subroutine getshot (ncid, filename, ntimes, time, numcols, 
     &      column_name, column_value)
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
	open (ncid,file=filename,status='old',form='formatted',err=100)
	read (ncid, '(13x,80a13)')  column_name
		do i=1,numcols
		column_name(i) = trim (adjustl (column_name(i)))
		if (debug > 0)  print *, column_name(i)
		enddo
	write (ntest,'(12x,80a9)') (column_name (i), i=1,numcols)
	numcolsout = numcols
	if (numcolsout > 22)  numcolsout = 22
	do i=1,ntimes
	read (ncid, '(80f13.1)')  time(i), (column_value (i,j), 
     &        j=1,numcols)
!>>>>>debug
!     read(ncid,'(80f13.1)')  time(i), column_value(i,1),(column_value (i,j), j=numcols-2,numcols-1), &
!                             (column_value (i,j), j=2,numcols-3),column_value(i,numcols)
	write (ntest,'(f8.3,1x,1p40e9.1)')  time(i), 
     &       (column_value (i,j), j=1,numcolsout)
!	print *, time(i), (column_value (i,j), j=1,5)
	enddo		
	print *, 'getshot: data read from : ', filename
	print *, 'getshot: ntimes, numcols, tmin, tmax =', 
     &     ntimes, numcols, time(1), time(ntimes)
      close(ncid)
	return
  100	print *, '>>> Cannot find file:', filename
	end
!---------------------------------------------------------------------
      subroutine D3fluxvals(ncid, filename, ntimes, time, numcols, 
     &       column_name, column_value)
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
      debug = 0
	open (ncid, file=filename, status='old', form='formatted',
     &        err=100)
	read (ncid, '(13x,80a13)')  column_name
		do i=1,numcols
		column_name(i) = trim (adjustl (column_name(i)))
		if (debug > 0)  print *, column_name(i)
		enddo
	write (ntest,'(12x,80a9)') (column_name (i), i=1,numcols)
	numcolsout = numcols
	if (numcolsout > 22)  numcolsout = 22
	do i=1,ntimes
	      read (ncid, '(80f13.1)')  time(i), (column_value (i,j), 
     &             j=1,numcols)
	write (ntest,'(f8.3,1x,1p40e9.1)')  time(i), 
     &       (column_value (i,j), j=1,numcolsout)
!	print *, time(i), (column_value (i,j), j=1,5)
	enddo		
	print *, 'getshot: data read from : ', filename
	print *, 'getshot: ntimes, numcols, tmin, tmax =', ntimes, 
     &      numcols, time(1), time(ntimes)
      close(ncid)
	return
  100	print *, '>>> Cannot find file:', filename
	end
!---------------------------------------------------------------------
      subroutine getlocang (ncid, filename, numloops, rloop, zloop, 
     &         tangle, pangle)
!
!			Read NSTX shot flux & B loop location files
!
	integer :: ntest = 11
	integer :: ncid, ntimes, numloops
	character*80  filename
	character*80  column_name 
	character*15  loop_name 
	real*8        rloop(numloops), zloop(numloops), 
     &   tangle(numloops), pangle(numloops)  
!
!
	rloop  = 0.d0 ; 	zloop  = 0.d0				! initialize arrays
	tangle = 0.d0 ; 	pangle = 0.d0				! initialize arrays
! 
        ncid = 71
	open (ncid,file=filename,status='old',form='formatted',err=100)
	read (ncid, '(a80)')  column_name
	write (ntest,'(a80)') column_name
        do i=1,numloops
	  read (ncid, '(a15, 8f15.1)',end=200)  loop_name, pangle(i), 
     &    rloop(i), zloop(i), tangle(i)
	  write (ntest,'(a15,8f15.3)' )  loop_name, pangle(i), 
     &    rloop(i), zloop(i), tangle(i)
        enddo
        close(ncid)
!		
	print *, 'loop location & angle data read from : ', filename
	print *, 'numloops =', numloops
	return
  100	print *, '>>> Cannot find file:', filename
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
	open (ncid,file=filename,status='old',form='formatted',err=100)
	read (ncid, '(a80)')  column_name
	write (ntest,'(a80)') column_name
        do i=1,numloops
	  read (ncid, '(a15, 8f15.1)',end=200)  loop_name, pangle, 
     &       rloop(i), zloop(i)
	  write (ntest,'(a15,8f15.3)' )  loop_name, pangle, rloop(i), 
     &       zloop(i) 
        enddo
        close(ncid)
!		
	print *, 'loop location & angle data read from : ', filename
	print *, 'numloops =', numloops
	return
  100	print *, '>>> Cannot find file:', filename
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
	real*8        rloop(numloops), zloop(numloops), 
     &       tangle(numloops), pangle(numloops)  
!
!
	rloop  = 0.d0 ; 	zloop  = 0.d0				! initialize arrays
	tangle = 0.d0 ; 	pangle = 0.d0				! initialize arrays
!
      ncid = 71
	open (ncid, file=filename, status='old', form='formatted',
     &        err=100)
	read (ncid, '(a80)')  column_name
	write (ntest,'(a80)') column_name
		do i=1,numloops
	read (ncid, '(a10,2f9.0)',end=200)  loop_name, rloop(i), 
     &         zloop(i)
	write (ntest,'(a10,2f9.4)' )  loop_name, rloop(i), zloop(i)
		enddo		
	print *, 'loop location & angle data read from : ', filename
	print *, 'numloops =', numloops
      close(ncid)
	return
  100	print *, '>>> Cannot find file:', filename
  	return
      stop 10
  200	print *, 'End of file:', filename
  	return
      stop 11
	end


