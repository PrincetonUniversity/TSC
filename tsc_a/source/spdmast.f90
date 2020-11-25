      subroutine spdmast
 
      USE ezcdf
      USE CLINAM, ONLY : numargs, filename, suffix
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
 
      character*50,   dimension(:),   allocatable :: current_name
      integer,        dimension(:),   allocatable :: current_active
      REAL(KIND=r8), dimension(:,:), allocatable :: current_data
 
      integer                 :: ncid, varid, shotnumber, ifirst
      integer                 :: ialloc, ier, i, debug
      integer                 :: ntimes,ncurrent,numflux,nmirnov
      integer                 :: numcurrent
      integer                 :: nnedl
      integer, dimension(3)   :: dimlens = (/1,1,1/)
      character*4             :: type
      data ifirst/1/
      if(ifirst.eq.0) go to 1000
      ifirst = 0
 
!---------------------------------------------------------------------
      debug = 0
 
!---------------------------------------------------------------------
!-- open the file
      if(numargs .lt. 1) then
        filename = 'mast_magnetics.nc'
      else
        filename = 'mast_magnetics.nc' // '.' // trim(suffix)
      end if
      call cdfopn(ncid,trim(filename),'r')
!
!---read the shot number
      call cdfGetVar(ncid,'shotnumber',shotnumber)
!
!---------------------------------------------------------------------
!-- read the time vector
      call cdfInqVar(ncid,'time',dimlens,type)
      ALLOCATE(time(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'time',time)
      ntimes = dimlens(1)
      if (debug .eq. 1) then
         print*,'TIME ---->'
         do i=1,dimlens(1)
            print*,time(i)
         end do
      endif
 
!---------------------------------------------------------------------
!-- read the power supply names
      call cdfInqVar(ncid,'psupply_name',dimlens,type,ier)
      ALLOCATE(psupply_name(dimlens(2)), STAT=ialloc )
      call cdfGetVar(ncid,'psupply_name',psupply_name)
      ncurrent = dimlens(2)
      if (debug .eq. 1) then
         print*,'POWER SUPPLY NAMES ---->'
         do i=1,dimlens(2)
            print*,i,trim(psupply_name(i))
         end do
      endif
 
!-- read the power supply currents
      call cdfInqVar(ncid,'psupply_current',dimlens,type,ier)
      ALLOCATE(psupply_current(dimlens(1),dimlens(2)), STAT=ialloc )
      call cdfGetVar(ncid,'psupply_current',psupply_current)
      if (debug .eq. 1) then
         print*,'POWER SUPPLY 1 CURRENT ---->'
         do i=1,dimlens(1)
            print*,psupply_current(i,16)
         end do
      endif
!---------------------------------------------------------------------
!-- read the flux loop names
      call cdfInqVar(ncid,'flux_name',dimlens,type,ier)
      ALLOCATE(flux_name(dimlens(2)), STAT=ialloc )
      call cdfGetVar(ncid,'flux_name',flux_name)
      numflux = dimlens(2)
      if (debug .eq. 1) then
         print*,'FLUX LOOP NAMES ---->'
         do i=1,dimlens(2)
            print*,i,flux_name(i)
         end do
      endif
 
!-- read the flux loop active status
      call cdfInqVar(ncid,'flux_active',dimlens,type,ier)
      ALLOCATE(flux_active(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'flux_active',flux_active)
      if (debug .eq. 1) then
         print*,'FLUX LOOP ACTIVE STATUS ---->'
         do i=1,dimlens(1)
            print*,i,flux_active(i)
         end do
      endif
 
!-- read the flux loop R position
      call cdfInqVar(ncid,'flux_r',dimlens,type,ier)
      ALLOCATE(flux_r(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'flux_r',flux_r)
      if (debug .eq. 1) then
         print*,'FLUX LOOP R ---->'
         do i=1,dimlens(1)
            print*,i,flux_r(i)
         end do
      endif
 
!-- read the flux loop Z position
      call cdfInqVar(ncid,'flux_z',dimlens,type,ier)
      ALLOCATE(flux_z(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'flux_z',flux_z)
      if (debug .eq. 1) then
         print*,'FLUX LOOP Z ---->'
         do i=1,dimlens(1)
            print*,i,flux_z(i)
         end do
      endif
 
!-- read the flux data
      call cdfInqVar(ncid,'flux_data',dimlens,type,ier)
      ALLOCATE(flux_data(dimlens(1),dimlens(2)), STAT=ialloc )
      call cdfGetVar(ncid,'flux_data',flux_data)
      if (debug .eq. 1) then
         print*,'FLUX LOOP 1 DATA ---->'
         do i=1,dimlens(1)
            print*,flux_data(i,1)
         end do
      endif
 
!---------------------------------------------------------------------
!-- read the mirnov loop names
      call cdfInqVar(ncid,'mirnov_name',dimlens,type,ier)
      ALLOCATE(mirnov_name(dimlens(2)), STAT=ialloc )
      call cdfGetVar(ncid,'mirnov_name',mirnov_name)
      nmirnov = dimlens(2)
      if (debug .eq. 1) then
         print*,'MIRNOV LOOP NAMES ---->'
         do i=1,dimlens(2)
            print*,i,mirnov_name(i)
         end do
      endif
 
!-- read the mirnov loop active status
      call cdfInqVar(ncid,'mirnov_active',dimlens,type,ier)
      ALLOCATE(mirnov_active(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'mirnov_active',mirnov_active)
      if (debug .eq. 1) then
         print*,'MIRNOV LOOP ACTIVE STATUS ---->'
         do i=1,dimlens(1)
            print*,i,mirnov_active(i)
         end do
      endif
 
!-- read the mirnov loop R position
      call cdfInqVar(ncid,'mirnov_r',dimlens,type,ier)
      ALLOCATE(mirnov_r(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'mirnov_r',mirnov_r)
      if (debug .eq. 1) then
         print*,'MIRNOV LOOP R ---->'
         do i=1,dimlens(1)
            print*,i,mirnov_r(i)
         end do
      endif
 
!-- read the mirnov loop Z position
      call cdfInqVar(ncid,'mirnov_z',dimlens,type,ier)
      ALLOCATE(mirnov_z(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'mirnov_z',mirnov_z)
      if (debug .eq. 1) then
         print*,'MIRNOV LOOP Z ---->'
         do i=1,dimlens(1)
            print*,i,mirnov_z(i)
         end do
      endif
 
!-- read the mirnov loop toroidal angle
      call cdfInqVar(ncid,'mirnov_tangle',dimlens,type,ier)
      ALLOCATE(mirnov_tangle(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'mirnov_tangle',mirnov_tangle)
      if (debug .eq. 1) then
         print*,'MIRNOV LOOP TOROIDAL ANGLE ---->'
         do i=1,dimlens(1)
            print*,i,mirnov_tangle(i)
         end do
      endif
 
!-- read the mirnov loop poloidal angle
      call cdfInqVar(ncid,'mirnov_pangle',dimlens,type,ier)
      ALLOCATE(mirnov_pangle(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'mirnov_pangle',mirnov_pangle)
      if (debug .eq. 1) then
         print*,'MIRNOV LOOP POLOIDAL ORIENTATION ---->'
         do i=1,dimlens(1)
            print*,i,mirnov_pangle(i)
         end do
      endif
 
!-- read the mirnov data
      call cdfInqVar(ncid,'mirnov_data',dimlens,type,ier)
      ALLOCATE(mirnov_data(dimlens(1),dimlens(2)), STAT=ialloc )
      call cdfGetVar(ncid,'mirnov_data',mirnov_data)
      if (debug .eq. 1) then
         print*,'MIRNOV LOOP 1 DATA ---->'
         do i=1,dimlens(1)
            print*,mirnov_data(i,1)
         end do
      endif
 
!---------------------------------------------------------------------
!-- read the current measurment names
      call cdfInqVar(ncid,'current_name',dimlens,type,ier)
      ALLOCATE(current_name(dimlens(2)), STAT=ialloc )
      call cdfGetVar(ncid,'current_name',current_name)
      numcurrent = dimlens(2)
      if (debug .eq. 1) then
         print*,'CURRENT MEASUREMENT NAMES ---->'
         do i=1,dimlens(2)
            print*,i,current_name(i)
         end do
      endif
 
!-- read the current measurment active status
      call cdfInqVar(ncid,'current_active',dimlens,type,ier)
      ALLOCATE(current_active(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'current_active',current_active)
      if (debug .eq. 1) then
         print*,'CURRENT MEASUREMENT ACTIVE STATUS ---->'
         do i=1,dimlens(1)
            print*,i,current_active(i)
         end do
      endif
 
!-- read the current data
      call cdfInqVar(ncid,'current_data',dimlens,type,ier)
      ALLOCATE(current_data(dimlens(1),dimlens(2)), STAT=ialloc )
      call cdfGetVar(ncid,'current_data',current_data)
        if(debug.eq.1) then
         print*,'CURRENT MEASUREMENT 10 DATA ---->'
            write(*,1001)(current_data(i,10),i=1,dimlens(1))
 1001       format(1x,1p10e12.4)
        endif
 
 
!-- read the line integrated density trace
      call cdfInqVar(ncid,'midplane_nedl',dimlens,type,ier)
      ALLOCATE(midplane_nedl(dimlens(1)), STAT=ialloc )
      call cdfGetVar(ncid,'midplane_nedl',midplane_nedl)
      nnedl = dimlens(1)
      if (debug .eq. 1) then
         print*,'LINE INTEGRATED DENSITY ---->'
         do i=1,dimlens(1)
            print*,midplane_nedl(i)
         end do
      endif
 
!---------------------------------------------------------------------
!-- close the file
      call cdfcls(ncid)
 1000 continue
!
!
!.....
      call spdmasto(shotnumber,time,ntimes,                              &  
     &     psupply_current,ncurrent,                                     &  
     &     flux_active,flux_r,flux_z,flux_data,numflux,                  &  
     &     mirnov_active,mirnov_r,mirnov_z,mirnov_tangle, mirnov_pangle,  &  
!    &                                                                   &  
     &                   mirnov_data,nmirnov,                            &  
     &     current_active,current_data,numcurrent ,                      &  
     &     midplane_nedl,nnedl )
      return
 
 
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
