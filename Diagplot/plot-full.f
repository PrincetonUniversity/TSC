      Program test
c
c     implicit real*8 (a-h,o-z)
c     include 'comfile.i'
      parameter(numfluxloops=44)
      dimension xplot(10000),yplot(10000),zplot(10000),
     1 eplot(10000), fplot(10000),gplot(10000),hplot(1000)
      dimension pflux(numfluxloops),eflux(numfluxloops),
     1          error(numfluxloops),icolor(numfluxloops),
     2          flux_r(numfluxloops),flux_z(numfluxloops),
     3          iflux_active(numfluxloops)
      dimension pfluxplt(10000,numfluxloops),
     1          efluxplt(10000,numfluxloops)
      dimension csuma(10),ppcur(10),fbcur(10),div(2)
      character*4 istring4
      character*14 label,istring14
      character*2 inum(75)
      data (inum(i),i=1,75)/
     1 " 1"," 2"," 3"," 4"," 5"," 6"," 7"," 8"," 9",
     1 "10","11","12","13","14","15","16","17","18","19",
     1 "20","21","22","23","24","25","26","27","28","29",
     1 "30","31","32","33","34","35","36","37","38","39",
     1 "40","41","42","43","44","45","46","47","48","49",
     1 "50","51","52","53","54","55","56","57","58","59",
     1 "60","61","62","63","64","65","66","67","68","69",
     1 "70","71","72","73","74","75"/
c
c
c
c....1.0  initialize constants and open files
      pi = 3.1415926535
c
      no167a = 9
      open(7,file='plotout',status='unknown',iostat=ios)
      open(no167a,file='tscnstxdata_out',status='old',iostat=ios22)
      call ncarcgm(1,"out.cgm")
      call dders(-1)
c
c.....initialize NCAR GKS graphics
      iwkid = 2
      iwtype = 8
      lunit = 2
      ierrf = 6
c     call gopks(ierrf,1000)
      call gopwk(iwkid,lunit,iwtype)
      call gacwk(iwkid)
c
c.....initialize program by ploting flux loop positions
c
c.....first pass through data to determine error catagories
      do i=1,numfluxloops
      error(i) = 0.
      enddo
c
      do irec=1,10000
      read(no167a,7010,end=102) times,plcurs,vespplass,vescurs,plcure,
     1                          vescure, vespplase,chicur
      read(no167a,7010) (csuma(k),k=1,10)
      read(no167a,7010) (ppcur(k),k=1,10)
      read(no167a,7010) (fbcur(k),k=1,10)
      read(no167a,7010) (pflux(i),i=1,numfluxloops)
      read(no167a,7010) (eflux(i),i=1,numfluxloops)
      read(no167a,7011)
      do i=1,numfluxloops
      error(i) = error(i) + (pflux(i)-eflux(i))**2
      enddo
      enddo
 102  continue
      do i=1,numfluxloops
      error(i) = sqrt(error(i)/irec)
      icolor(i) = 2
      if(error(i).gt. 0.090) icolor(i) = 3
      if(error(i).lt. 0.025) icolor(i) = 1
      enddo
      write(6,1005) ((i,error(i)),i=1,numfluxloops)
 1005 format(3(i5,1pe10.2))
c    
      call readmag(icolor,flux_r,flux_z,iflux_active,numfluxloops)
 1000 continue
      write(6,1002)
      read(5,*)nn
      call frame(0)
c
      if(nn.lt.0) go to 2000
c

c
      if(nn.gt.100 .and. nn.lt.200) then
      label = "current no"
      label(11:12) = inum(nn-100)
      endif        
c
c
      if(nn.gt.0 .and. nn.le.numfluxloops) then
      label = "flux loop "
      label(11:12) = inum(nn)
      endif
c
      if(nn.eq.0) then
      label = "vessel current"
      endif
      if(nn.eq.201) then
      label = "CHI current"
      endif
c
      write(6,1010) nn,label
 1010 format("nn = ",i3,3x,a14)
c
      rewind(no167a)
c
      ymin = 0.
      ymax = 0.
      do irec=1,10000
      read(no167a,7010,end=101) times,plcurs,vespplass,vescurs,plcure,
     1       vescure, vespplase, chicure,vescure2,chicurs
      read(no167a,7010) (csuma(k),k=1,10)
      read(no167a,7010) (ppcur(k),k=1,10)
      read(no167a,7010) (fbcur(k),k=1,10)
      read(no167a,7010) (pflux(i),i=1,numfluxloops)
      read(no167a,7010) (eflux(i),i=1,numfluxloops)
      read(no167a,7011)
c    
c
      do i=1,numfluxloops
      pfluxplt(irec,i) = pflux(i)
      efluxplt(irec,i) = eflux(i)
      enddo
c
      xplot(irec) = times

      if(nn.eq.200) then
      yplot(irec) = csuma(1)/100000.
      zplot(irec) = csuma(2)/1000.
      eplot(irec) = csuma(3)/1000.
      fplot(irec) = csuma(5)/1000.
                    endif
c
      if(nn.gt.100 .and. nn.lt.200) then
      yplot(irec) = csuma(nn-100)
      zplot(irec) = ppcur(nn-100)
      fplot(irec) = fbcur(nn-100)
      endif
c
c
c
      if(nn.gt.0 .and. nn.le.numfluxloops) then
      yplot(irec) = pflux(nn)
      zplot(irec) = eflux(nn)
      endif
      if(nn.eq.201) then
      yplot(irec) = chicurs
      zplot(irec) = chicure
      endif
c
      if(nn.eq.0) then
c     yplot(irec) = vespplass
c     zplot(irec) = vespplase - ppcur(7)
c     zplot(irec) = vespplase
      eplot(irec) = plcure
      fplot(irec) = plcurs
      gplot(irec) = vescure
      hplot(irec) = vescurs
      endif
c
      ymin = amin1(ymin,yplot(irec),zplot(irec))
      ymax = amax1(ymax,yplot(irec),zplot(irec))
      if(nn.eq.0) then
      ymin = amin1(ymin,eplot(irec),fplot(irec))
      ymax = amax1(ymax,eplot(irec),fplot(irec))
      ymin = amin1(ymin,gplot(irec),hplot(irec))
      ymax = amax1(ymax,gplot(irec),hplot(irec))
      endif
      if(nn.eq.200) then
c
      ymin = -380.
      ymax =  180.
      endif
c
      imax = irec
      enddo
      close(no167a)
 101  continue
      write(6,3333) imax,ymin,ymax
 3333 format(" imax,ymin,ymax =",i5,1p2e12.4)
      if(nn.eq.300) go to 300
c
      bmax = .89
      bmin = .20
      if(nn.eq.200) then
      bmax = .49
      bmin = .20
      endif
      call colora("blue")
      call maps(xplot(1),xplot(imax),ymin,ymax,.20,.99,bmin,bmax)
      if(nn.eq.0) go to 198
      if(nn.eq.200) go to 200
      call tracec(1hS,xplot,yplot,imax,-1,-1,0.,0.)
      call colora("red")
      call tracec(1hE,xplot,zplot,imax,-1,-1,0.,0.)
      call colora("blue")
      call tracec(1hS,xplot,yplot,imax,-1,-1,0.,0.)
c
c      if(nn.gt.100.and. nn.lt.200) then
c      call tracec(1hF,xplot,fplot,imax,-1,-1,0.,0.)
c      call tracec(1hF,xplot,fplot,imax,-1,-1,0.,0.)
c      endif
c
  198 continue
      if(nn.eq.0) then
      call colora("red")
      call tracec(1hE,xplot,eplot,imax,-1,-1,0.,0.)
      call tracec(1hS,xplot,fplot,imax,-1,-1,0.,0.)
      call colora("green")
      call tracec(1hS,xplot,hplot,imax,-1,-1,0.,0.)
      call tracec(1hE,xplot,gplot,imax,-1,-1,0.,0.)
      call tracec(1hE,xplot,gplot,imax,-1,-1,0.,0.)
      call colora("blue")
      endif
c
c     bottom label
      x = xplot(1) + 0.3*(xplot(imax)-xplot(1))
      y = ymin     - .15*(ymax       -ymin    )
      icase = 1
      isize = 3
      iorient = 0
      istring4 = "time"
      nchar = 4
      call setlch(x,y,icase,isize,iorient,-1)
      call gtext(istring4,nchar,0)
c
c     left label
      x = xplot(1) - .15*(xplot(imax)-xplot(1))
      y = ymin     + .20*(ymax       -ymin    )
      icase = 1
      isize = 3
      iorient = 1
      istring14 = label
      nchar = 14
      call setlch(x,y,icase,isize,iorient,-1)
      call gtext(istring14,nchar,0)
      go to 201
 200  continue
c
      call trace(xplot,yplot,imax,-1,-1,0.,0.)
      call colora("red")
      call trace(xplot,zplot,imax,-1,-1,0.,0.)
      call colora("green")
      call trace(xplot,eplot,imax,-1,-1,0.,0.)
      call colora("blue")
      call trace(xplot,fplot,imax,-1,-1,0.,0.)
      call trace(xplot,fplot,imax,-1,-1,0.,0.)
      call colora("blue")
      go to 201
 300  continue
c
      rmax = 0.
      zmax = 0.
      rmin = 0.
      zmin = 0.
      do i=1,numfluxloops
      rmax = amax1(rmax,1.3*flux_r(i))
      rmin = amin1(rmin,    flux_r(i))
      zmax = amax1(zmax,    flux_z(i))
      zmin = amin1(zmin,1.4*flux_z(i))
      enddo
      rdis = 0.3*(rmax-rmin)
      zdis = 0.2*(zmax-zmin)
c
      iflux_active(7) = 0
      iflux_active(5) = 0
      iflux_active(10) = 0
      iflux_active(12) = 0
      iflux_active(14) = 0
      iflux_active(19) = 0
      iflux_active(25) = 0
      iflux_active(28) = 0
      iflux_active(29) = 0
      iflux_active(31) = 0
      iflux_active(21) = 0
      iflux_active(23) = 0
      iflux_active(24) = 0
      iflux_active(26) = 0
      iflux_active(30) = 0
      iflux_active(31) = 0
      iflux_active(32) = 0
      iflux_active(35) = 0
      iflux_active(37) = 0

      fluxmin = 0.
      fluxmax = 0.
      do 303 j=1,numfluxloops
      if(iflux_active(j).eq.0) go to 303
      do 302 i=1,imax
      fluxmin = amin1(fluxmin,pfluxplt(i,j),efluxplt(i,j))
      fluxmax = amax1(fluxmax,pfluxplt(i,j),efluxplt(i,j))
 302  continue
 303  continue
c
      do 301 i=2,numfluxloops

      r1norm = (flux_r(i)     -rmin)/(rmax-rmin)
      r2norm = (flux_r(i)+rdis-rmin)/(rmax-rmin)
      z1norm = (flux_z(i)-zdis-zmin)/(zmax-zmin)
      z2norm = (flux_z(i)     -zmin)/(zmax-zmin)
c 
      if(iflux_active(i).eq.0) go to 301
      call map(xplot(1),xplot(imax),fluxmin,fluxmax,
     1          r1norm,r2norm,z1norm,z2norm)
      call colora("red")
      call trace(xplot,pfluxplt(1,i),imax,-1,-1,0.,0.)
      call colora("green")
      call trace(xplot,efluxplt(1,i),imax,-1,-1,0.,0.)
      call colora("blue")
      call setlch(xplot(1),fluxmax,0,1,0,-1)
      call gtext(inum(i),2,0)
      write(6,1001) i,flux_r(i),flux_z(i),error(i)
 1001 format(i3,1p4e12.4)
c
 301  continue
      call colora("blue")
      call maps(xplot(1),xplot(imax),fluxmin,fluxmax,
     1         .7,1.,.8,1.)
 201  continue
c  
c
c.....empty buffers
      call guwk(2,0)

c
      go to 1000
 2000 continue
      call gclrwk(iwkid,1)
      call gdawk(iwkid)
      call gclwk(iwkid)
      call plote

      stop
c
 7010 format(1p10e12.4)
 7011 format(1x) 
 1002 format(/,"type -1:quit",
     .       /,"      0: I(vv)",
     .       /,"      #: flux loop #",
     .       /,"  100+n: cur n",
     1       /,"    200:special plot of currents",
     2       /,"    201:chi currents",
     3       /,"    300:summary of flux loops" )
      end

      subroutine readmag(icolor,flux_r,flux_z,flux_active,numflux)
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
      integer,        dimension(:),   allocatable :: flux_active
      REAL(KIND=r8), dimension(:)  :: flux_r
      REAL(KIND=r8), dimension(:)  :: flux_z
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
      integer                 :: ntest = 11
      integer                 :: numnumerrec_shotfile
      integer                 :: nhalosensors
!
      integer                 :: ncid, varid, shotnumber, ifirst
      integer                 :: ialloc, ier, i, debug ,l, n, jj
      integer                 :: ntimes,ncurrent,numflux,nmirnov
      integer                 :: numcurrent
      integer                 :: ntfwindings
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
         			filename = 'nstx_shot.dat'
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
	call getshot (ncid, coilcur_filename, ntimes, time, ncurrent, psupply_name, psupply_current)
	if(nhalosensors.gt.0) call getshot (ncid, halocur_filename, ntimes, time, nhalosensors, halosensor_name, halosensor_current)
        if(imachine.eq.0) then
!.........NSTX
	  if(numflux.gt.0) call getshot (ncid, fluxloopdata_filename, ntimes, time, numflux, flux_name, flux_data)
        else
!.........DIII-D
           if(numflux.gt.0) call D3fluxvals(ncid, fluxloopdata_filename, ntimes, time, numflux, flux_name, flux_data)
        endif
!
	if(nmirnov.gt.0) call getshot (ncid, bloopdata_filename,    ntimes, time, nmirnov, mirnov_name, mirnov_data)
!---------------------------------------------------------------------
!				Read flux and B loop locations and angles
        if(debug.gt.0) print *, imachine,numflux, 'imachine and numflux before call to getlocang'
!
        if(imachine.eq.0) then
!.........NSTX
	  if(numflux.gt.0) call getloc (ncid, fluxlooploc_filename, numflux, flux_r, flux_z)
        else
!.........DIII-D
          if(numflux.gt.0) call d3fluxloops( ncid, fluxlooploc_filename,numflux, flux_r, flux_z)
        endif
	if(nmirnov.gt.0) call getlocang (ncid, blooploc_filename,    nmirnov, mirnov_r, mirnov_z, mirnov_tangle, mirnov_pangle)
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
	print *, 'here is call to spdnstxo:'
	print *, 'shotnumber, ntimes, ncurrent, numflux, nmirnov, numcurrent, nhalosensors, nnedl, ndifl, nbi1p'
	print *, shotnumber, ntimes, ncurrent, numflux, nmirnov, numcurrent, nhalosensors, nnedl, ndifl, nbi1p
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
	print *, 'getshot: data read from : ', filename
	print *, 'getshot: ntimes, numcols, tmin, tmax =', ntimes, numcols, time(1), time(ntimes)
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
	print *, 'getshot: data read from : ', filename
	print *, 'getshot: ntimes, numcols, tmin, tmax =', ntimes, numcols, time(1), time(ntimes)
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
	print *, 'loop location & angle data read from : ', filename
	print *, 'numloops =', numloops
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
	print *, 'loop location & angle data read from : ', filename
	print *, 'numloops =', numloops
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
	print *, 'loop location & angle data read from : ', filename
	print *, 'numloops =', numloops
      close(ncid)
	return
  100	print *, 'Error opening file:', filename
  	return
      stop 10
  200	print *, 'End of file:', filename
  	return
      stop 11
	end


