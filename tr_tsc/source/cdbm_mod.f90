!#include "f77_dcomplx.h"
      subroutine cdbm_mod(chicdbm)
!
      USE trtsc, ONLY : use_user_rot,use_user_pres
      USE CLINAM
!      USE BALCLI
!      USE RUNAWAY
      USE SAPROP
!      USE SCR3
!      USE NEWPLOT

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,npsihm,k,i,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sum, AREAL, ptoti, felec, fion
!============      
      INTEGER :: istat = 0 
!============      
      logical :: first_read_rot=.true.
      logical :: first_read_rott=.true.
      logical :: first_read_rotp=.true.
      logical :: first_read_pres=.true.
      logical :: ex_rot=.false. 
      logical :: ex_rott=.false. 
      logical :: ex_rotp=.false. 
      logical :: ex_pres=.false. 
      real*8 :: tint, rint, wrot1, wrot2, frot,tfluxd,x1,x2,            &
     &          pres1,pres2,presint
      real*8, dimension(:),allocatable :: time_rot, rad_rot
      real*8, dimension(:),allocatable :: time_rotp, rad_rotp
      real*8, dimension(:),allocatable :: time_rott, rad_rott
      real*8, dimension(:),allocatable :: time_pres, rad_pres
      real*8, dimension(:,:),allocatable :: wrot,vrott,vrotp,prestot
      integer :: j1, j2, kk, l, lsav
      integer :: ntime_rot, nrad_rot
      integer :: ntime_rott, nrad_rott
      integer :: ntime_rotp, nrad_rotp
      integer :: ntime_pres, nrad_pres
!============  
!     arrays needed for the CDBM model
      integer jj,icdbmodel,jmin,jlim_cdbm,jmin_cdbm
      REAL*8 xchi,chiemksCT,chiimksCT,fs_min,chi_jlim,qmin_CT
      INTEGER isw,npsiCDBM,npsiCT
      REAL*8, DIMENSION(npsit) :: dpdr,dpidr,dptdr,dErdr,vrot_exp,      &
     &                            chicdbm,vrott_exp,vrotp_exp,Erad,     &
     &                            chicdbm0,angrotp_exp,pres_exp,        &
     &                            cdbm_fs,cdbm_fe,cdbm_wexb
!============  
!----------------------------------------------------------------------------
!     toroidal rotation profile
      do j=1,npsit 
         angrotp_exp(j) = acoef(4986)*float(npsit-j+1)/float(npsit)
         vrot_exp(j) = angrotp_exp(j)*rmajora(j)
         if (j .eq. 1) vrot_exp(j)=angrotp_exp(j)*xmag
         vrott_exp(j) = 0.0
         vrotp_exp(j) = 0.0
      end do
!        
!  read rotation data from file, if provided
      if(use_user_rot .and. first_read_rot) then
        inquire(file="user_rot_data", exist=ex_rot)
        if (ex_rot) then
          open(55,file="user_rot_data",form="formatted",status="old")
          read(55,*) ntime_rot, nrad_rot
          allocate(time_rot(ntime_rot), rad_rot(nrad_rot))
          allocate(wrot(nrad_rot,ntime_rot))
          read(55,*) time_rot(1:ntime_rot)
          read(55,*) rad_rot(1:nrad_rot)
          read(55,*) wrot(1:nrad_rot,1:ntime_rot)
          close (55)
        endif
        first_read_rot=.false.
      endif

      if(use_user_rot .and. first_read_rott) then
        inquire(file="user_rot_tor_data", exist=ex_rott)
        if (ex_rott) then
          open(55,file="user_rot_tor_data",form="formatted",status="old")
          read(55,*) ntime_rott, nrad_rott
          allocate(time_rott(ntime_rott), rad_rott(nrad_rott))
          allocate(vrott(nrad_rott,ntime_rott))
          read(55,*) time_rott(1:ntime_rott)
          read(55,*) rad_rott(1:nrad_rott)
          read(55,*) vrott(1:nrad_rott,1:ntime_rott)
          close (55)
        endif
        first_read_rott=.false.
      endif
      
      if(use_user_rot .and. first_read_rotp) then
        inquire(file="user_rot_pol_data", exist=ex_rotp)
        if (ex_rotp) then
          open(55,file="user_rot_pol_data",form="formatted",status="old")
          read(55,*) ntime_rotp, nrad_rotp
          allocate(time_rotp(ntime_rotp), rad_rotp(nrad_rotp))
          allocate(vrotp(nrad_rotp,ntime_rotp))
          read(55,*) time_rotp(1:ntime_rotp)
          read(55,*) rad_rotp(1:nrad_rotp)
          read(55,*) vrotp(1:nrad_rotp,1:ntime_rotp)
          close (55)
        endif
        first_read_rotp=.false.
      endif
!
      if(use_user_rot .and. ex_rott) then
        do  kk = 1,ntime_rott-1
          if(times .gt. time_rott(kk) .and. times .le. time_rott(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_rott(j1))/(time_rott(j2)-time_rott(j1))
            exit
          endif
        enddo
        if(times .le. time_rott(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif 
        if(times .gt. time_rott(ntime_rott)) then
          j1 = ntime_rott-1
          j2 = ntime_rott
          tint = 1.
        endif
        do i=1,npsit
           tfluxd = sqrt((float(i-1)*dpsi)/(float(npsit-1)*dpsi))
           do ii=1,nrad_rott-1
              x1 = (rad_rott(ii)-xmag)/rminor
              x2 = (rad_rott(ii+1)-xmag)/rminor
              if(tfluxd .gt. x1 .and. tfluxd .le. x2) then
                 rint = (tfluxd-x1)/(x2-x1)
                 wrot1 = vrott(ii,j1)+(vrott(ii+1,j1)-vrott(ii,j1))*rint
                 wrot2 = vrott(ii,j2)+(vrott(ii+1,j2)-vrott(ii,j2))*rint
                 frot = wrot1 + (wrot2-wrot1)*tint
                 exit
              endif
           enddo
           vrott_exp(i) = frot
        end do
      else
        do i=1,npsit
           vrott_exp(i) = 0.0_R8
        enddo
      endif
!----------------------------------------------------------------------------
      if(use_user_rot .and. ex_rotp) then
        do  kk = 1,ntime_rotp-1
          if(times .gt. time_rotp(kk) .and. times .le. time_rotp(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_rotp(j1))/(time_rotp(j2)-time_rotp(j1))
            exit
          endif
        enddo
        if(times .le. time_rotp(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif 
        if(times .gt. time_rotp(ntime_rotp)) then
          j1 = ntime_rotp-1
          j2 = ntime_rotp
          tint = 1.
        endif
        do i=1,npsit
           tfluxd = sqrt((float(i-1)*dpsi)/(float(npsit-1)*dpsi))
           do ii=1,nrad_rotp-1
              x1 = (rad_rotp(ii)-xmag)/rminor
              x2 = (rad_rotp(ii+1)-xmag)/rminor
              if(tfluxd .gt. x1 .and. tfluxd .le. x2) then
                 rint = (tfluxd-x1)/(x2-x1)
                 wrot1 = vrotp(ii,j1)+(vrotp(ii+1,j1)-vrotp(ii,j1))*rint
                 wrot2 = vrotp(ii,j2)+(vrotp(ii+1,j2)-vrotp(ii,j2))*rint
                 frot = wrot1 + (wrot2-wrot1)*tint
                 exit
              endif
           enddo
           vrotp_exp(i) = frot
        end do
      else
        do i=1,npsit
           vrotp_exp(i) = 0.0_R8
        enddo
      endif
!----------------------------------------------------------------------------
      if(use_user_rot .and. ex_rot) then
        do  kk = 1,ntime_rot-1
          if(times .gt. time_rot(kk) .and. times .le. time_rot(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_rot(j1))/(time_rot(j2)-time_rot(j1))
            exit
          endif
        enddo
        if(times .le. time_rot(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif 
        if(times .gt. time_rot(ntime_rot)) then
          j1 = ntime_rot-1
          j2 = ntime_rot
          tint = 1.
        endif
        do i=1,npsit
           tfluxd = sqrt((float(i-1)*dpsi)/(float(npsit-1)*dpsi))
           do ii=1,nrad_rot
              if(tfluxd .gt. rad_rot(ii) .and. tfluxd .le. rad_rot(ii+1)) then
                 rint = (tfluxd-rad_rot(ii))/(rad_rot(ii+1)-rad_rot(ii))
                 wrot1 = wrot(ii,j1) + (wrot(ii+1,j1)-wrot(ii,j1))*rint
                 wrot2 = wrot(ii,j2) + (wrot(ii+1,j2)-wrot(ii,j2))*rint
                 frot = wrot1 + (wrot2-wrot1)*tint
                 exit
              endif
           enddo
           angrotp_exp(i) = frot
           vrot_exp(i) = angrotp_exp(i)*rmajora(i)
        end do
      endif

!  read total pressure data (including fast ions) from file, if provided
      if(use_user_pres .and. first_read_pres) then
        inquire(file="user_pres_data", exist=ex_pres)
        if (ex_pres) then
          open(55,file="user_pres_data",form="formatted",status="old")
          read(55,*) ntime_pres, nrad_pres
          allocate(time_pres(ntime_pres), rad_pres(nrad_pres))
          allocate(prestot(nrad_pres,ntime_pres))
          read(55,*) time_pres(1:ntime_pres)
          read(55,*) rad_pres(1:nrad_pres)
          read(55,*) prestot(1:nrad_pres,1:ntime_pres)
          close (55)
        endif
        first_read_pres=.false.
      endif
!----------------------------------------------------------------------------
      if(use_user_pres .and. ex_pres) then
        do  kk = 1,ntime_pres-1
          if(times .gt. time_pres(kk) .and. times .le. time_pres(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_pres(j1))/(time_pres(j2)-time_pres(j1))
            exit
          endif
        enddo
        if(times .le. time_pres(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif
        if(times .gt. time_pres(ntime_pres)) then
          j1 = ntime_pres-1
          j2 = ntime_pres
          tint = 1.
        endif
        do i=1,npsit
           tfluxd = sqrt((float(i-1)*dpsi)/(float(npsit-1)*dpsi))
           if(tfluxd .le. rad_pres(1)) then
              rint = 0.0 
              pres1 = prestot(1,j1)
              pres2 = prestot(1,j2)
              presint= pres1 + (pres2-pres1)*tint
           endif     
           if(tfluxd .gt. rad_pres(nrad_pres)) then
              pres1 = prestot(nrad_pres,j1)
              pres2 = prestot(nrad_pres,j2)
              presint= pres1 + (pres2-pres1)*tint
           endif     
           do ii=1,nrad_pres
              if(tfluxd .gt. rad_pres(ii) .and. tfluxd .le. rad_pres(ii+1)) then
                 rint = (tfluxd-rad_pres(ii))/(rad_pres(ii+1)-rad_pres(ii))
                 pres1 = prestot(ii,j1) + (prestot(ii+1,j1)-prestot(ii,j1))*rint
                 pres2 = prestot(ii,j2) + (prestot(ii+1,j2)-prestot(ii,j2))*rint
                 presint= pres1 + (pres2-pres1)*tint
                 exit
              endif     
           enddo
           pres_exp(i) = presint*usdp
        end do
      else
        do i=1,npsit
           pres_exp(i) = 0.
        enddo
      endif



!----------------------------------------------------------------------------
!   calculate gradients used in the CDBM model if itrmod=16 
      icdbmodel = INT(acoef(3200))
      npsihm = int(pwidthc * AREAL(npsit))
      call cdbm_dpdr(pres_exp,dpdr,dpidr,dptdr)  ! calculate pressure gradient profiles
!        calculate radial electric field shear profile
      SELECT CASE(icdbmodel)
         CASE (0,1)
            do j=1,npsit
               Erad(j) = 0.0_R8
               dErdr(j) = 0.0_R8
               cdbm_wexb(j) = 0.0_R8
            end do
         CASE (2,3,4,5)
            call cdbm_dErdr(dpidr,vrot_exp,vrott_exp,vrotp_exp,Erad,cdbm_wexb)
      END SELECT
      if(use_user_pres .and. ex_pres) then
         call cdbm_chi(icdbmodel,dpdr,dptdr,cdbm_wexb,chicdbm0,cdbm_fs,cdbm_fe,chicdbm)
      else
         call cdbm_chi(icdbmodel,dpdr,dpdr,cdbm_wexb,chicdbm0,cdbm_fs,cdbm_fe,chicdbm)
      endif
      if (acoef(3206) .gt. 0.) then
         do l=1,100
            if(times .ge. tpest(l) .and. times .lt. tpest(l+1)) lsav=l
         end do
         if(times- tpest(lsav) .le. 1.0E-02) then 
            write(nout,1069) kcycle,times
            write(nout,1070)
            do j=1,npsit-2
              write(nout,2069) j*dpsi,rminora(j),dpidr(j),vrot_exp(j),vrott_exp(j),       &
     &                       vrotp_exp(j),Erad(j),cdbm_wexb(j),chicdbm0(j), &
     &                       cdbm_fs(j),cdbm_fe(j),chicdbm(j)
      
           end do
         endif
 1069 format(1h1," cycle",i7,"    time",1pe12.4)  
 1070 format("rho  dpi/dr   vrot   vtorNC    vpolNC   Erad   wExB",     &  
     & "  chi0   fs    fe    chicdbm")
 2069 format(1p12e11.3)
      endif
!-----------------------------------------------------------------------

      end 



