Subroutine w20solv ( zz, zvr, zvi, iu, eu, imax, emax, G_ave, kps, ierr )
  use w20data
  implicit none

  !Maximum number of iterations
  integer, parameter :: nitmax = 50

  !Eigenvectors resulting from solved system of equations
  real*8, intent(inout), dimension(neq,neq) :: zvr, zvi

  !Complex frequency of the mode 
  complex*16, intent(inout), dimension(neq) :: zz

  !Pointing indices (im-->unstable mode used, iu-->number of unstable modes)
  integer, intent(inout) :: iu, eu, imax, emax

  !External error checking flag
  integer, intent(inout) :: ierr

  !Geometric average needed to compute momentum pinches
  real*8 :: G_ave, G_ave_tmp, kps, kps_tmp

  !Internal error checking flag for eigenvalue solver
  integer :: ifail, kflag, it_err, it_conv
  
  !Internal counters for number of unstable modes and unstable mode used
  integer :: im, nu

  !Current number of iterations
  integer :: niter, niterm

  !Iteration indices
  integer :: i, j, nmod

  !Skip mode flag
  integer :: iskip, foundmode

  !Complex frequencies of fastest growing eigenmode
  complex*16, dimension(nitmax+1) :: wz, xbest, wzfnd

  !Input frequency for dispersion relation
  complex*16 :: wzin, dw1, dw

  !Temporary storage for growth rates and eigenmodes
  real*8                     :: wimax, wemax, wamax
  real*8, dimension(neq,neq) :: zvrtmp, zvitmp
  complex*16, dimension(neq) :: zztmp

!  complex*16, external :: cmplex

  !Mask array for active modes
!  integer, dimension(neq) :: mask_mod
  real*8,  dimension(neq) :: mask

  !Variables dealing with error control
  real*8 ::  zepsilon, ztol
  real*8 ::  delWZ, delZmin, delZ

  !Pointer to mode being used
  integer :: nmod_ptr

  !Initial count of unstable modes
  integer :: unst

  !Final test for largest growing mode in the other direction
  real :: womax, woomax, wiimax
  integer :: omax, ou, au, amax

  !Matrices and vectors for eigenvalue system
  real*8, dimension(neq) :: omega, gamma, zbeta, zalfi, zalfr, somega, sgamma
  real*8, dimension(neq,neq) :: zamr, zami, zbmr, zbmi
  
  real*8 :: zb

!  OPEN (UNIT=10,FILE='gamma.out')
!  OPEN (UNIT=11,FILE='gamma2.out')

  au = 0
  iu = 0
  eu = 0
  ou = 0
  unst = 0
  amax = 0
  imax = 0
  emax = 0
  omax = 0
  ierr = 0
  zepsilon = 1.0E-6
  ztol = 1.0E-3
  zz(1:neq)=(0.0,0.0)
  mask(1:neq)=0.0
!  mask_mod(1:neq)=0
  wamax = 0.0
  wimax = 0.0
  wemax = 0.0
  womax = 0.0
  wiimax = 0.0
  woomax = 0.0
  it_err = 0
  it_conv = 0
!  write (6,*) "Enter dispersion relation solver:"

  modes:  do nmod=1,neq

!     write (6,*) "Processing mode number",nmod
!     write (6,*) "======================"

     !Initialization of some parameters
     nmod_ptr = 0
     niter = 0
     ierr  = 0
     wz(1:nitmax+1)=(0.0,0.0)
     wzin=(0.0,0.0)
     wzfnd(1:nitmax+1)=(0.0,0.0)
     zztmp(1:neq)=(0.0,0.0)
     xbest(1:nitmax+1)=(0.0,0.0)
     delZ = 0.0
     delWZ = 0.0
     dw = (0.0,0.0)
     dw1 = (0.0,0.0)

     call w20wguess( wzin )

        if ( searchmode .eq. S_ELC .and. dreal(wzin).lt.0.0 ) then
           wzin=wzin+2*abs(dreal(wzin))
        end if

     wz(1)    = wzin
     xbest(1) = wzin     

     iskip = 0

     !Algorithm to find converged frequency begins here
     freq: Do
        !Initialize complex frequencies

        niter = niter + 1

        omega(1:neq)=0.0
        gamma(1:neq)=0.0


        !Maximum number of iterations 
        If ( niter .gt. nitmax ) Then
           it_err = 1
           If ( 0.1*abs(delWZ) .lt. ztol ) then
              exit freq
           !Hard error if tolerance not found within an order of magnitude
           Else If ( 0.1*abs(delWZ) .ge. ztol ) then
              ierr = 2
              it_err = 2
              Cycle modes
           End If
        End If

        !Obtain dispersion relation
        call w20disp(zamr, zami, zbmr, zbmi, wzin, G_ave_tmp, kps_tmp )

        !Diagonalize equation system
        call r8tomsqz(neq,neq,zamr,zami,zbmr,zbmi,zalfr,zalfi,zbeta,zvrtmp,zvitmp,ifail)

        !Check for errors
        If ( ifail .gt. 0 ) Then
           Write(6,*) "In w20solv: eror in eigenvalue solver"
           ierr = 1
           exit modes
        End If

        zztmp = w20omg ( zalfr, zalfi, zbeta, neq )
        omega = dreal(zztmp)
        gamma = dimag(zztmp)

        !Begin carry out a series out tests to improve convergence and reliability

        !Test 1: if all modes are unstable, then return from this call
        unst = w20nunst ( zztmp, neq, S_ANY )

        if ( unst .eq. 0 .and. niter .ge. 3) then 
           it_conv = 1
           cycle modes
        end if

        !Test 2: lock onto the mode being looked at
        !        It should be nmod, but sometimes the solver rotates eigenvalues
        !        Thus find the mode that is most like the mode last used
        !        In the first iteration, choose the mode to be nmod
        if ( niter .eq. 1 ) then
           nmod_ptr = nmod
        else
           !For niter > 1, lock onto the mode last used
!           nmod_ptr = w20wmatch ( wzin, zztmp, neq )
           nmod_ptr = w20wmatch ( wzfnd(niter-1), zztmp, neq )
        end if

        !Save the matched mode for the next iteration
        wzfnd(niter)=zztmp(nmod_ptr)

        !Test 3: if the mode is stable, go to the next mode
        If ( gamma(nmod_ptr) .lt. 0.0 .and. niter .ge. 3 ) then
           it_conv = 1
           cycle modes
        end if

        foundmode = w20s_fnd ( zztmp(nmod_ptr) )

        !Test 4: if the mode is in the wrong direction, go to the next mode
        if ( foundmode .ne. searchmode .and. niter .gt. 3 ) then
           it_conv = 1
           cycle modes
        end if
        !Choose the next trial frequency as the average of the old and the new mode
        if ( niter .le. 1 ) then
           wz(niter+1)   = zztmp(nmod_ptr)
        else
           wz(niter+1)   = 0.5*(wzin+(0.D0,1.D0)*gamma(nmod_ptr)+omega(nmod_ptr))
        end if


        wzin          = wz(niter+1)
        Xbest(niter+1)= wzin

!        write (6,*) "Next mode is", wz(niter+1)

!        If ( niter .gt. 1 ) CALL AITKEN(wz,Xbest,niter+1,nitmax+1)       
        
        !Compute estimate of error
        delWZ = w20delwz( Xbest(niter+1), Xbest(niter) )

!        write (6,*) "Computed error is", delWz

        !Exit inner loop (freq) when convergence is found
        if ( delWZ .lt. ztol .and. niter .ge. 2) then
!           write (6,*) "Mode converged, modes found are:"
!           do j=1,niter
!              write (6,*) j, wzfnd(j)
!           end do
           it_conv = 1
           Exit freq
        end if

     End Do freq

     !If the unstable mode is in the direction we are looking for
     ! then evaluate whether this mode is faster growing than the 
     ! previously stored mode. If so, carry out the replacement

     foundmode = w20s_fnd ( zztmp(nmod_ptr) )
     
     if ( foundmode .eq. searchmode ) then
        if ( w20gamma( zztmp, neq, nmod_ptr ) .gt. wamax ) then       
           wamax = w20gamma( zztmp, neq, nmod_ptr )
           amax  = nmod_ptr
           G_ave = G_ave_tmp
           kps   = kps_tmp

           do j=1,neq
              foundmode = w20s_fnd( zztmp(j) )
              if ( j .ne. amax .and. foundmode .eq. searchmode .and. w20gamma( zztmp, neq, j ) .gt. wamax ) then
!                 zztmp(j)=zztmp(j)-(0.0,1.0)*zztmp(j)
              end if
           end do

           imax  = w20wmunst ( zztmp, neq, S_ION ) !nmod_ptr
           wimax = w20gamma ( zztmp, neq, imax ) !wamax
           emax  = w20wmunst ( zztmp, neq, S_ELC )
           wemax = w20gamma ( zztmp, neq, emax )

           eu  = w20nunst( zztmp, neq, S_ELC )
           iu  = w20nunst( zztmp, neq, S_ION )

           zz  = zztmp
           zvr = zvrtmp
           zvi = zvitmp
        end if
              
     end if

  End Do modes

  if ( it_conv .ne. 1 ) then
     Write (6,*) "In w20solv: Excessive iterations in Weiland model"    

     if ( it_err .gt. 1 ) then
!        write (6,*)"Mode not converged, modes found are:"

!        do j=1,niter
!           write (6,*) j, wzfnd(j)
!        end do

        write (6,*) "mode=",nmod,"  ptr=",nmod_ptr
        write (6,*) "delWZ=",delWZ
        write (6,*) "GRDNE=",gne
        write (6,*) "GRDTE=",gte
        write (6,*) "GRDNI=",gni
        write (6,*) "GRDTI=",gti
        write (6,*) "GRDNZ=",gnz
        write (6,*) "GRDTZ=",gtz
        write (6,*) "NE=",ne
        write (6,*) "TE=",te
        write (6,*) "NH=",nh
        write (6,*) "TI=",th
        write (6,*) "NZ=",nz
        write (6,*) "NS=",ns_ne*ne
        write (6,*) "Q=",q
        write (6,*) "Shear=",shear
        write (6,*) "Elong=",kappa
        write (6,*) "Aimp=",aimp
        write (6,*) "Zimp=",zimp
        write (6,*) "Ahyd=",ahyd
        write (6,*) "Rmaj=",rmaj
        write (6,*) "Rmin=",rmin
        write (6,*) "btor=",btor
        !              write (20,'(5ES15.6)') (wz(j),wzfnd(j),w20delwz(wz(j),wzfnd(j) ),j=1,niter)
        !              stop
     end if
  end if

 ! if (searchmode .eq. S_ION) then
 !    write(15,*)betae,gne,beta,wimax,dreal(zz(imax))
 !    write(6,*)"Vector found is:"
 !    do j=1,9
 !       write(6,*) zz(j)
 !    end do
 ! else if (searchmode .eq. S_ELC) then
 !    write(16,*)betae,gne,beta,wimax,dreal(zz(emax))
 !    write(6,*)"Vector found is:"
 !    do j=1,9
 !       write(6,*) zz(j)
 !    end do
 ! end if

!  stop

End Subroutine w20solv
