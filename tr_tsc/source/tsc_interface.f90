
      MODULE Wall_Input_Mod
      USE kind_spec_mod
      IMPLICIT none

      INTEGER :: nwall, istat
      REAL(KIND=rspec) :: Rwall(501), Zwall(501)

      CONTAINS
 
      subroutine get_wall_rz
      IMPLICIT none

      NAMELIST /Wall_Data/  nwall, Rwall, Zwall
      LOGICAL :: ex
     

      Nwall = 0
      Rwall = 0.0
      Zwall = 0.0

      inquire(file="wall_data", exist=ex)
      if(.not.ex) return
      close(51)
      open(51,file="wall_data",status="old",form="formatted",           &
     &     iostat=istat)
      if(istat .ne. 0) then
      write(6,*) "wall_data open error, use default instead"
      else
      read(51,nml=wall_data,iostat=istat)
      if (istat .ne. 0) then
      write(6,*) "wall_data read error, use default instead"
      Nwall = 0
      Rwall = 0.0
      Zwall = 0.0
      endif
      endif
      close(51)
      return

      END subroutine get_wall_rz

      END MODULE Wall_Input_Mod

      subroutine wrgeqdsk
      USE trtsc, ONLY : suffix2
      USE EZspline_obj
      USE EZspline
      USE Wall_Input_Mod
      USE CLINAM
      USE SAPROP
      USE SCR10
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nxpest,nzpest,n,i,jpest,j,jj,k,idum,nzdim,nwalliter
      INTEGER nlimiter,kk
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rwalliter,zwalliter,rlimiter,zlimiter,alrpest
      REAL*8 alzpest,rzpest,rpest,pmpest,plpest,dps,ps,pval,ppval
      REAL*8 gval,gpval,gppval,t1,t3,t5,rdim,zdim,zmid,bcentr,xdum
      REAL*8 zip,psival,pprimesav2,pprimesav1,ffprimsav2,ffprimsav1
!============
      DIMENSION rwalliter(41),zwalliter(41),rlimiter(14),zlimiter(14)
      CHARACTER*256 :: geqfile
!...first wall R, m
      data nwalliter /41/
      data rwalliter /                                                  &
     &     0.389699993E+01, 0.389699993E+01, 0.389699993E+01,           &
     &     0.389699993E+01, 0.389699993E+01,                            &
     &     0.389699993E+01, 0.390590332E+01, 0.422439195E+01,           &
     &     0.488341711E+01, 0.579541201E+01,                            &
     &     0.659269730E+01, 0.753578650E+01, 0.804297693E+01,           &
     &     0.841374349E+01, 0.862634702E+01,                            &
     &     0.844925396E+01, 0.802537144E+01, 0.734628092E+01,           &
     &     0.648383942E+01, 0.619649187E+01,                            &
     &     0.641517387E+01, 0.628078406E+01, 0.612278301E+01,           &
     &     0.597535286E+01, 0.585236581E+01,                            &
     &     0.576516494E+01, 0.541100001E+01, 0.542365547E+01,           &
     &     0.517086686E+01, 0.507371412E+01,                            &
     &     0.498017040E+01, 0.488168675E+01, 0.473817533E+01,           &
     &     0.454194599E+01, 0.420778989E+01,                            &
     &     0.426303195E+01, 0.433206076E+01, 0.432631467E+01,           &
     &     0.427385267E+01, 0.413935633E+01,                            &
     &     0.389699993E+01/
!...first wall Z, m
      data zwalliter /                                                  &
     &    -0.247399998E+01,-0.152300000E+01,-0.507000029E+00,           &
     &     0.509999990E+00, 0.152699995E+01,                            &
     &     0.254399991E+01, 0.361090915E+01, 0.442029307E+01,           &
     &     0.485028371E+01, 0.465868987E+01,                            &
     &     0.405595509E+01, 0.322615987E+01, 0.249932699E+01,           &
     &     0.169408710E+01, 0.644794793E+00,                            &
     &    -0.502509739E+00,-0.145731983E+01,-0.242725511E+01,           &
     &     -0.307962131E+01,-0.319189131E+01,                           &
     &    -0.331806128E+01,-0.335194210E+01,-0.336413279E+01,           &
     &     -0.342225498E+01,-0.352023336E+01,                           &
     &    -0.364202560E+01,-0.390221210E+01,-0.466697304E+01,           &
     &     -0.405860664E+01,-0.388478514E+01,                           &
     &    -0.380548145E+01,-0.375761886E+01,-0.376581725E+01,           &
     &     -0.380923193E+01,-0.395054759E+01,                           &
     &    -0.329789073E+01,-0.310480773E+01,-0.297296136E+01,           &
     &     -0.285143298E+01,-0.273474946E+01,                           &
     &    -0.247399998E+01/
!...limiter R, m
      data rlimiter /8.271_R8,8.258_R8,8.269_R8,8.275_R8,                &
     &               8.282_R8,8.283_R8,8.281_R8,8.276_R8,                &
     &               8.267_R8,8.255_R8,8.240_R8,8.222_R8,                &
     &               8.201_R8,8.205_R8/
!...limiter Z, m
      data zlimiter /-0.436_R8,-0.274_R8,-0.112_R8,0.000_R8,             &
     &               0.213_R8,0.375_R8,0.537_R8,0.700_R8,0.862_R8,       &
     &               1.024_R8,1.185_R8,1.347_R8,1.508_R8,1.670_R8/
!============
!     INTEGER :: istat = 0
      INTEGER :: ngeq, ii, ij
      REAL*8, ALLOCATABLE, DIMENSION(:) :: te2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ane2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ti2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: avez2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zeff2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: fpol
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pres
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ffprim
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pprime
      REAL*8, ALLOCATABLE, DIMENSION(:) :: qpsi
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: psiaux
!============

      TYPE (EZspline1_r8) :: spln
      real*8 tmpa(ppsi),tmp1,tmp2
      integer :: bdycon(2)
      logical :: l_first=.true.
!============
      IF(.not.ALLOCATED(te2)) ALLOCATE( te2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ane2)) ALLOCATE( ane2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ti2)) ALLOCATE( ti2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(avez2)) ALLOCATE( avez2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(zeff2)) ALLOCATE( zeff2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(fpol)) ALLOCATE( fpol(pnx), STAT=istat)
      IF(.not.ALLOCATED(pres)) ALLOCATE( pres(pnx), STAT=istat)
      IF(.not.ALLOCATED(ffprim)) ALLOCATE( ffprim(pnx), STAT=istat)
      IF(.not.ALLOCATED(pprime)) ALLOCATE( pprime(pnx), STAT=istat)
      IF(.not.ALLOCATED(qpsi)) ALLOCATE( qpsi(pnx), STAT=istat)
      IF(.not.ALLOCATED(psiaux)) ALLOCATE( psiaux(pnx,2*pnz),            &
     &                                 STAT=istat)
!============
      if (istat .ne. 0) stop 'Allocation Error : wrpest  '
!============
!
 2001 format(20x,8a10)
 6100 format(5i10)
 6101 format(1p5e20.12)
 460  format(2f20.4)
!.....write GEQDSK ASCII file in the format used by EFIT (added 6/04/03)
 2000 format(6a8,3i4)
 2020 format(5e16.9)
 2022 format(2i5)

      if(l_first) then
      call get_wall_rz
      if(nwall .eq. 0) then
!     if not defined, use input limiter data
      do i=1,nlim
      rwall(i)=xlima(i)
      zwall(i)=zlima(i)
      enddo
      nwall= nlim
      if(isym .eq. 1) then
      do j=1, nlim-2
      ii=nlim+j
      ij=nlim-j
      rwall(ii)= rwall(ij)
      zwall(ii)=-zwall(ij)
      enddo
      nwall=2*(nlim-1)
      endif
      endif
      if(nwall .gt. 501)                                                &
     &          stop "No. of 1st wall points exceeded allocation in geq"
      l_first=.false.
      endif

      ngeq = 52
      filename=trim(suffix2)//"_ps.geq"
      open(ngeq,file=trim(filename),form="formatted",status="unknown",   &
     &                           iostat=istat)
      
      rdim = alx - ccon
      zdim = 2*alz
      zmid = 0._R8      
      bcentr = gzero/xplas
      xdum = 0._R8
      idum = 13
      nzdim = nz + isym*(nz-1)
        do i=2,nx+1
          if(isym.ne.0) then
             do j=3,nz+1
        psiaux(i,j) = psi(i,j)
               psiaux(i,j+nz-1) = psi(i,j)
               psiaux(i,nz+3-j) = psi(i,j)
              enddo
               psiaux(i,nz+1) = psi(i,2)
          else
             do j=2,nz+1
        psiaux(i,j) = psi(i,j)
               psiaux(i,j) = psi(i,j)
             enddo
          endif
        enddo
      write(ngeq,2000) (name(i),i=1,6),idum,nx,nzdim
      write(ngeq,2020) rdim,zdim,xplas,ccon,zmid
      write(ngeq,2020) xmag,zmag,psimin,psilim,bcentr
!    
!     global(75) is beta based on vacuum toroidal field
!     global(76) is central beta x 0.5
!     global(174) is troyon coefficient
!
      zip = tcurdtp*tpi*udsi
      write(ngeq,2020) zip,psimin,global(75),xmag,global(76)*2._R8
      write(ngeq,2020) zmag,global(174),psilim,xdum,xdum

!     bdycon = 0
!     call EZspline_init(spln,npsit,bdycon,ierr)
!     call EZspline_error(ierr)
!     spln%x1(1:npsit) = xsv2(1:npsit)
!     tmpa(1:npsit)=qprof2(1:npsit)
!     call EZspline_setup(spln,tmpa,ierr)
!     call EZspline_error(ierr)
   
      do 250 i=1,nx
      psival = psimin + (i-1)*(psilim-psimin)/(nx-1)
      call geval(psival,2,fpol(i),gpval,ffprim(i),imag,jmag)
      call peval(psival,2,pres(i),pprime(i),imag,jmag)
!     tmp1=psival
!     call EZspline_interp(spln,tmp1,tmp2,ierr)
!     qpsi(i)=tmp2
!     call EZspline_error(ierr)
      do 249 jj=2,npsit
      j = jj
      if(psival .ge. xsv2(j-1) .and. psival .le. xsv2(j)) go to 251
 249  continue
 251  continue
      qpsi(i) = qprof2(j-1) + (psival-xsv2(j-1))                        &
     &        /(xsv2(j)-xsv2(j-1))*(qprof2(j)-qprof2(j-1))
!     if(j .le. 2) then
!     tmp1=(psival-xsv2(3))*(psival-xsv2(4))
!     tmp3=(psival-xsv2(4))*(psival-xsv2(2))
!     tmp5=(psival-xsv2(2))*(psival-xsv2(3))
!     tmp2=(xsv2(2)-xsv2(3))*(xsv2(2)-xsv2(4))
!     tmp4=(xsv2(3)-xsv2(4))*(xsv2(3)-xsv2(2))
!     tmp6=(xsv2(4)-xsv2(2))*(xsv2(4)-xsv2(3))
!     qpsi(i)=tmp1/tmp2*qprof2(2)+                                      &
!    &        tmp3/tmp4*qprof2(3)+                                      &
!    &        tmp5/tmp6*qprof2(4)
!     endif
 250  continue
!     call EZspline_free(spln,ierr)
!     call EZspline_error(ierr)
!
!.....temporary patch
      pprimesav2 = pprime(2)
      pprimesav1 = pprime(1)
      ffprimsav2 = ffprim(2)
      ffprimsav1 = ffprim(1)
!
      pprime(2) = 2._R8*pprime(3) - pprime(4)
      pprime(1) = 2._R8*pprime(2) - pprime(3)
      ffprim(2) = 2._R8*ffprim(3) - ffprim(4)
      ffprim(1) = 2._R8*ffprim(2) - ffprim(3)
!     write(nterm,6666) pprimesav1,pprimesav2,ffprimsav1,ffprimsav2
!     write(nterm,6666) pprime(1) ,pprime(2) ,ffprim(1) ,ffprim(2)
!6666 format("***diag***  ",1p4e12.4)
!
!.....end of temporary patch
      write(ngeq,2020) (fpol(i),i=1,nx)
      write(ngeq,2020) (udsp*pres(i),i=1,nx)
      write(ngeq,2020) (ffprim(i),i=1,nx)
      write(ngeq,2020) (udsp*pprime(i),i=1,nx)
      write(ngeq,2020) ((psiaux(i,j),i=2,nxp),j=2,nzdim+1)

!     temporary fix
      if(qpsi(1) .le. 0.0) then
      qpsi(1)=2.0*qpsi(2)-qpsi(3)
      endif

      write(ngeq,2020) (qpsi(i),i=1,nx)
      nlimiter = kmax
      if(isym.eq.1) then
        nlimiter = 2*(kmax-2)
          do kk=kmax+2,nlimiter
          k = kk-kmax
          xplot(1,kmax+k) = xplot(1,kmax-k)
          zplot(1,kmax+k) =-zplot(1,kmax-k)
          enddo
        endif
!
      write(ngeq,2022) nlimiter,nwall
      write(ngeq,2020) (xplot(1,k),zplot(1,k),k=1,nlimiter)
      write(ngeq,2020) (rwall(i),zwall(i),i=1,nwall)
      close(ngeq)
!!cj write wall_data for plasma state
!      write(76,"('&wall_data')")
!      write(76,"('nwall=',I3)") nwall
!      write(76,"('rwall=')") 
!      write(76,7991) (rwall(i),i=1,nwall)
!      write(76,"('zwall=')")
!      write(76,7991) (zwall(i),i=1,nwall)
! 7991 format(1p3(2x,e15.7))
!      write(76,"('/')")

      filename=trim(suffix2)//"_ps.jso"
      open(ngeq,file=trim(filename),form="formatted",status="unknown",   &
     &                           iostat=istat)
      do 100 j=1,npsit
      ps = xsv2(j)
      call peval(ps,1+isurf,pval,ppval,imag,jmag)
      prpest2(j) = pval
      pppest2(j) = ppval
      te2(j)    = .5_R8*(te(j)   +te(j+1))
      ane2(j)   = .5_R8*(ane(j)  +ane(j+1))
      ti2(j)    = .5_R8*(ti(j)   +ti(j+1))
      avez2(j)  = .5_R8*(avez(j) +avez(j+1))
      zeff2(j) = .5_R8*(zeffa(j)+zeffa(j+1))
      if(j.eq.1 .or. j.eq.npsit) go to 100
      t1 = rdpsi*tpi*(qprof2(j))**2/(xmja2(j))**2
      t3 = gxmja(j+1)*xmja(j+1)/(.5_R8*(qprof2(j)+qprof2(j+1)))
      t5 = gxmja(j)*xmja(j)/(.5_R8*(qprof2(j)+qprof2(j-1)))
      ajpest2(j) = t1*(t3-t5)
  100 continue
      ajpest2(1) = 2._R8*ajpest2(2) - ajpest2(3)
      ajpest2(npsit) = 0._R8
      te2(1)    = 1.5_R8*te(2)    - 0.5_R8*te(3)
      ane2(1)   = 1.5_R8*ane(2)   - 0.5_R8*ane(3)
      ti2(1)    = 1.5_R8*ti(2)    - 0.5_R8*ti(3)
      avez2(1)  = 1.5_R8*avez(2)  - 0.5_R8*avez(3)
      zeff2(1) = 1.5_R8*zeffa(2) - 0.5_R8*zeffa(3)

!......write equilibrium data in ASCII form (changed 9/22/02)
      write(ngeq,2001) (name(i),i=1,8)
      write(ngeq,6100) kcycle,isym,ipest,npsit,kmax
      write(ngeq,6101) times,xmag,zmag,gzero,apl,beta,                  &
     & betapol,ali2,qsaw,psimin,psilim
      write(ngeq,6101) (prpest2(j),j=1,npsit)
      write(ngeq,6101) (pppest2(j),j=1,npsit)
      write(ngeq,6101) (ajpest2(j),j=1,npsit)
      write(ngeq,6101) (xsv2(j),   j=1,npsit)
      write(ngeq,6101) (xplot(1,k),k=1,kmax+1)
      write(ngeq,6101) (zplot(1,k),k=1,kmax+1)
      write(ngeq,6101) (te2(j),    j=1,npsit)
      write(ngeq,6101) (ane2(j),   j=1,npsit)
      write(ngeq,6101) (ti2(j),    j=1,npsit)
      write(ngeq,6101) (avez2(j),  j=1,npsit)
      write(ngeq,6101) (zeff2(j), j=1,npsit)
      close(ngeq)


      RETURN
      END

      subroutine get_ps_data(nb,nicrf,nlhrf,necrf)
      use plasma_state_mod
      integer :: nb, nicrf, nlhrf, necrf

      nb=ps%nbeam
      nicrf=ps%nicrf_src
      nlhrf=ps%nlhrf_src
      necrf=ps%necrf_src

      return
      end  



      subroutine wrxpls
      USE trtsc
      USE tsc_ps_mod
      USE plasma_state_eq_mod
      USE plasma_state_pa_mod
!     USE plasma_state_mod 
      USE clinam
      USE saprop
      USE radtab
      USE specie
      IMPLICIT NONE
      real*8 :: ajtzon, abszon, afwzon, alhzon, apszon, acdzon, aeczon

      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      character*150 zfile
      integer ns,nt1,irz,ier,str_length
      character*20 label
      real*8 fbdy,crat
      real*8 zadjust

      REAL*8, ALLOCATABLE, DIMENSION(:) :: te2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ane2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ti2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zeff2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: savebre2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: savecyc2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: saveimp2

      integer :: i, j, k, istat, id, id_fun
      integer :: id_rho, num_rho, isize
      real*8, allocatable, dimension(:) :: tmp1, tmp2, tmp3, tmp4
      real*8 :: tol
      integer :: ii, itimes
      logical :: l_beam=.false., l_icrf=.false., l_lh=.false., l_ecrf=.false.
      logical :: l_first=.true.
      integer :: nchrgs, lsav, n, nspec
      real*8 :: beam_power,denom,ffact,dvol,total,dum1,dum2,ppval,pval
      real*8 :: gps,rmid
      integer :: ps_nbeam, ps_nicrf_src, ps_nlhrf_src, ps_necrf_src
      real*8 :: term, fac, curtot, powtot

      eq%t0 = acoef(ipos)
      if(acoef(ipos) .ge. 1000000.0) eq%t0=times
      eq%t1 = eq%t0+acoef(ipos+1)


!     scalar
      pa%vsur=global(13)
      pa%taup=acoef(93)

      call get_ps_data(ps_nbeam,ps_nicrf_src,ps_nlhrf_src, ps_necrf_src)
      if(ps_nbeam.gt.0 .and. .NOT.ALLOCATED(pa%power_nbi)) ALLOCATE(pa%power_nbi(ps_nbeam),STAT=istat)
      if(ps_nicrf_src.gt.0 .and. .NOT.ALLOCATED(pa%power_ic))  ALLOCATE(pa%power_ic(ps_nicrf_src),STAT=istat)
      if(ps_nlhrf_src.gt.0 .and. .NOT.ALLOCATED(pa%power_lh))  ALLOCATE(pa%power_lh(ps_nlhrf_src),STAT=istat)
      if(ps_necrf_src.gt.0 .and. .NOT.ALLOCATED(pa%power_ec))  ALLOCATE(pa%power_ec(ps_necrf_src),STAT=istat)
      pa%nbeam = ps_nbeam
      pa%power_nbi=0.0
      l_beam=.false.
      do j=1,ntpts
!cj    write(*,*) "--- tsc_interface L430 beamp(l) =", j, beamp(j) 
       if(beamp(j) .gt. 0) l_beam=.true.
      enddo
      if(l_beam) then
        do j=1,ntpts-1
          lsav = j
          if(tpro(j).le.time .and. tpro(j+1).gt.time) exit
        enddo
        if(ps_nbeam .gt. 0) then
          do j=1,ps_nbeam
!cj         write(*,*) "--- tsc_interface L440 beamp(l) =", lsav, beamp(lsav) 
            pa%power_nbi(j) = beamp(lsav)/ps_nbeam
      write(*,1810) kcycle, time, lsav,ps_nbeam,j,pa%power_nbi(j)
 1810 format("kcycle, time, lsav,ps_nbeam,j,pa%power_nbi(j)",i7,1pe12.4,3i5,1pe12.4)
          enddo
        else
          write(6,*) "error, ps_nbeam = 0"
          if(kcycle.gt.0) ineg=61
        endif
      endif
!
      pa%nicrf=max(1,ps_nicrf_src)
      pa%power_ic=0.0
      l_icrf=.false.
      do j=1,ntpts
      if(picrh(j).gt.0) then
      l_icrf=.true.
      exit
      endif
      enddo
      if(l_icrf) then
      do j=1,ntpts-1
      lsav = j
      if(tpro(j).le.time .and. tpro(j+1).gt.time) exit
      enddo
      denom = tpro(lsav+1)-tpro(lsav)
      if(denom.eq.0) denom = 1.
      ffact = (time-tpro(lsav))/denom
      do j=1,pa%nicrf
      pa%power_ic(j) = picrh(lsav)/pa%nicrf
      enddo
      endif

      pa%nlh=max(1,ps_nlhrf_src)
      pa%power_lh=0.0
      l_lh=.false.
      do j=1,ntpts
      if(plhamp(j).gt.0) then
      l_lh=.true.
      exit
      endif
      enddo
      if(l_lh) then
      do j=1,ntpts-1
      lsav = j
      if(tpro(j).le.time .and. tpro(j+1).gt.time) exit
      enddo
      denom = tpro(lsav+1)-tpro(lsav)
      if(denom.eq.0) denom = 1.
      ffact = (time-tpro(lsav))/denom
      do j=1,pa%nlh
      pa%power_lh(j) = plhamp(lsav)/pa%nlh
      enddo
      endif

      pa%necrf=max(1,ps_necrf_src)
      pa%power_ec=0.0
      l_ecrf=.false.
      do j=1,ntpts
      if(pecrh(j).gt.0) then
      l_ecrf=.true.
      exit
      endif
      enddo
      if(l_ecrf) then
      do j=1,ntpts-1
      lsav = j
      if(tpro(j).le.time .and. tpro(j+1).gt.time) exit
      enddo
      denom = tpro(lsav+1)-tpro(lsav)
      if(denom.eq.0) denom = 1.
      ffact = (time-tpro(lsav))/denom
      do j=1,pa%necrf
      pa%power_ec(j) = pecrh(lsav)/pa%necrf
      enddo
      endif


!     profiles

      IF(.NOT.ALLOCATED(te2)) ALLOCATE(te2(ppsi))
      IF(.NOT.ALLOCATED(ane2)) ALLOCATE(ane2(ppsi))
      IF(.NOT.ALLOCATED(ti2)) ALLOCATE(ti2(ppsi))
      IF(.NOT.ALLOCATED(zeff2)) ALLOCATE(zeff2(ppsi))
      IF(.NOT.ALLOCATED(savebre2)) ALLOCATE(savebre2(ppsi))
      IF(.NOT.ALLOCATED(savecyc2)) ALLOCATE(savecyc2(ppsi))
      IF(.NOT.ALLOCATED(saveimp2)) ALLOCATE(saveimp2(ppsi))
      te2=0.0d0
      ti2=0.0d0
      ane2=0.0d0
      zeff2=0.0d0
      savebre2=0.0d0
      savecyc2=0.0d0
      saveimp2=0.0d0
      IF(.NOT.ALLOCATED(tmp2)) ALLOCATE(tmp2(ppsi))
      tmp2=0.0

      IF (l_first) THEN
      isize=ppsi
      ALLOCATE(pa%rho(isize),stat=istat)
      pa%rho=0.0
      ALLOCATE(pa%ne(isize),stat=istat)
      pa%ne=0.0
      ALLOCATE(pa%nhy(isize),stat=istat)
      pa%nhy=0.0
      ALLOCATE(pa%nhe(isize),stat=istat)
      pa%nhe=0.0
      ALLOCATE(pa%nimp(isize,imp%nimpchg),stat=istat)
      pa%nimp=0.0
      ALLOCATE(pa%te(isize),stat=istat)
      pa%te=0.0
      ALLOCATE(pa%ti(isize),stat=istat)
      pa%ti=0.0
      ALLOCATE(pa%zeff(isize),stat=istat)
      pa%zeff=0.0
      ALLOCATE(pa%vpars(isize),stat=istat)
      pa%vpars=0.0
      ALLOCATE(pa%prad(isize),stat=istat)
      pa%prad=0.0
      ALLOCATE(pa%prad_br(isize),stat=istat)
      pa%prad_br=0.0
      ALLOCATE(pa%prad_cy(isize),stat=istat)
      pa%prad_cy=0.0
      ALLOCATE(pa%prad_li(isize),stat=istat)
      pa%prad_li=0.0
      ALLOCATE(pa%adi(isize),stat=istat)
      pa%adi=0.0
      ALLOCATE(pa%chie(isize),stat=istat)
      pa%chie=0.0
      ALLOCATE(pa%chii(isize),stat=istat)
      pa%chii=0.0
      ALLOCATE(pa%chio(isize),stat=istat)
      pa%chio=0.0
      ALLOCATE(pa%d_s(isize),stat=istat)
      pa%d_s=0.0
      ALLOCATE(pa%eta_p(isize),stat=istat)
      pa%eta_p=0.0
      ALLOCATE(pa%vloop(isize),stat=istat)
      pa%vloop=0.0
      ALLOCATE(pa%omegat(isize),stat=istat)
      pa%omegat=0.0
      ALLOCATE(pa%pres(isize),stat=istat)
      pa%pres=0.0
      ALLOCATE(pa%qprof(isize),stat=istat)
      pa%qprof=0.0
      ALLOCATE(pa%gasfl_ion(eq%nspec),stat=istat)
      pa%gasfl_ion=0.0
      ALLOCATE(pa%recyc_ion(eq%nspec),stat=istat)
      pa%recyc_ion=0.0
      ALLOCATE(pa%pohm(isize),stat=istat)
      pa%pohm=0.0
      ALLOCATE(pa%curbs(isize),stat=istat)
      pa%curbs=0.0
      ALLOCATE(pa%qie(isize),stat=istat)
      pa%qie=0.0
      ALLOCATE(pa%curoh(isize),stat=istat)
      pa%curoh = 0.0
      ALLOCATE(pa%pelh(isize),stat=istat)
      pa%pelh = 0.0
      ALLOCATE(pa%curlh(isize),stat=istat)
      pa%curlh = 0.0
      ENDIF

      do j=2, npsit
      te2(j)    = .5*(te(j)   +te(j+1))
      ane2(j)   = .5*(ane(j)  +ane(j+1))
      ti2(j)    = .5*(ti(j)   +ti(j+1))
      zeff2(j)  = .5*(zeffa(j)+zeffa(j+1))
      dvol = vary(j)-vary(j-1)
      savebre2(j) = savebre2(j-1)+savebre(j)*dvol
      savecyc2(j) = savecyc2(j-1)+savecyc(j)*dvol
      saveimp2(j) = saveimp2(j-1)+(saveimp(j) +                         &
     &              sradion(j)*usdp/usdt)*dvol
      enddo
      te2(1)    = 2.0*te(2)-1.0*te(3)
      ane2(1)   = 2.0*ane(2)-1.0*ane(3)
      ti2(1)    = 2.0*ti(2)-1.0*ti(3)
      zeff2(1) = 2.0*zeffa(2) - 1.0*zeffa(3)

      te2 = te2*1.D-3     ! convert to keV
      ti2 = ti2*1.D-3     ! convert to keV
      savebre2 = savebre2*udsp/udst
      savecyc2 = savecyc2*udsp/udst
      saveimp2 = saveimp2*udsp/udst

      pa%nrho = npsit
      do i = 1, npsit
      pa%rho(i) = SQRT(FLOAT(i-1)/FLOAT(npsit-1))
      enddo
      pa%rho(1)=0.0d0
      pa%rho(npsit)=1.0d0

      pa%te(:npsit) = te2(:npsit)
      pa%ti(:npsit) = ti2(:npsit)
      pa%ne(:npsit) = ane2(:npsit)
      pa%qprof(:npsit) = qprof2(:npsit)
      pa%zeff(:npsit) = zeff2(:npsit)
      pa%vpars(:npsit) = 0.0
      pa%prad_br(:npsit) = savebre2(:npsit)
      pa%prad_cy(:npsit) = savecyc2(:npsit)
      pa%prad_li(:npsit) = saveimp2(:npsit)
      pa%prad(:npsit)=pa%prad_br(:npsit)+                               &
     &                pa%prad_cy(:npsit)+                               &
     &                pa%prad_li(:npsit)

      do j=2, npsit
      tmp2(j)    = .5*(anhy(j)   +anhy(j+1))
      enddo
      tmp2(1)    = 2.0*tmp2(2)-1.0*tmp2(3)
      pa%nhy(:npsit) = tmp2(:npsit)
      pa%frac_h = min(1.0, max(0.0, acoef(4955)))                           !H fraction
      pa%frac_d = (1.0-pa%frac_h)*acoef(113)                                !D fraction
      pa%frac_t = (1.0-pa%frac_h)*(1.0-acoef(113))                          !T fraction
!
!>>>>>debug
!     write(6,*) pa%nhy(:npsit)

      do j=2, npsit
      tmp2(j)    = .5*(anhe(j)   +anhe(j+1))
      enddo
      tmp2(1)    = 2.0*tmp2(2)-1.0*tmp2(3)
      pa%nhe(:npsit) = tmp2(:npsit)

      pa%gasfl_ion=0.0
      pa%recyc_ion=0.0
      total=0.0
      do j=2, npsit
      total= total+anhy(j)*(vary(j)-vary(j-1))
      enddo
      if(acoef(93) .gt. 0.0) then
      i=0
      if(eq%tsc_spec(1) .gt. 0) then
      i=i+1
      pa%recyc_ion(i)=pa%frac_h*total/acoef(93)
      endif
      if(eq%tsc_spec(2) .gt. 0) then
      i=i+1
      pa%recyc_ion(i)=pa%frac_d*total/acoef(93)
      endif
      if(eq%tsc_spec(3) .gt. 0) then
      i=i+1
      pa%recyc_ion(i)=pa%frac_t*total/acoef(93)
      endif
      endif
      total=0.0
      do j=2, npsit
      total= total+anhe(j)*(vary(j)-vary(j-1))
      enddo
      if(acoef(93) .gt. 0.0) then
      if(eq%tsc_spec(4) .gt. 0) then
      i=i+1
      pa%recyc_ion(i)=total/acoef(93)
      endif
      endif

      if(iimp .eq. 0) then
      tmp2=0.0
      total=0.0
      do k=2, npsit
      dum1=(ane(k) - zgas*anhy(k) - 2.*anhe(k))/zimp
      dum2=(ane(k+1) - zgas*anhy(k+1) - 2.*anhe(k+1))/zimp
      total=total+dum1*(vary(k)-vary(k-1))
      tmp2(k)=0.5*(dum1+dum2)
      enddo
      tmp2(1)=2.0*tmp2(2)-1.0*tmp2(3)
      pa%nimp(:npsit,1)=tmp2(:npsit)
      if(acoef(93) .gt. 0.0) then
      pa%recyc_ion(eq%nspec-imp%nimpchg+1)=total/acoef(93)
      endif

      else
      i=0
      do j=1, pimp
      nchrgs=nchrgsr(j)
      if(.not.imp_present(j)) cycle
!     if( j.ne.4 .and. j .ne. 2 .and. j.ne.7) cycle
!     if( j.ne.4 ) cycle
!     if( j.ne.4 .and. j .ne. 2 ) cycle
      do n=2, nchrgs
      tmp2=0.0
      total=0.0
!     do n=2, nchrgs
      do k=2, npsit
      dum1=nq(n,j,k)/vp(k)
      dum2=nq(n,j,k+1)/vp(k+1)
      tmp2(k)=tmp2(k)+0.5d0*(dum1+dum2)
      total = total + dum1*(vary(k)-vary(k-1))
      enddo                           !k
!     enddo                           !n
      tmp2(1)=2.0*tmp2(2)-1.0*tmp2(3)
      i=i+1
      pa%nimp(:npsit,i)=tmp2(:npsit)
      if(acoef(93) .gt. 0.0) then
      pa%recyc_ion(eq%nspec-imp%nimpchg+i)=total/acoef(93)
      endif
      enddo                           !n
      enddo                           !j
      endif

      do j=1, npsit
      call peval(xsv2(j),1+isurf,pval,ppval,imag,jmag)
      pa%pres(j) = pval*udsp
      enddo

      pa%chie(1:npsit) = chiesec(1:npsit)
      pa%chii(1:npsit) = chiisec(1:npsit)
      do j=1, npsit
      gps = gja2(j)/vp2(j)*(tpi*qprof2(j))
      rmid = 0.5*(adn(j+1)/vp(j+1)+adn(j)/vp(j))
      pa%d_s(j)=0.0
      if(gps .ne. 0.0) then
      pa%d_s(j)=cs1(j)*rmid/(gps*udst)
      endif
      enddo
      pa%d_s(1)=pa%d_s(2)
      
      do j=1, npsit
      pa%eta_p(j)=etpara(j)*udsr
      enddo
      do j=1, npsit
      pa%vloop(j)=as(j)*udsv
      enddo
      
      do j=2, npsit
      tmp2(j)    = .5*(adi(j)   +adi(j+1))
      enddo
      tmp2(1)    = 2.0*tmp2(2)-1.0*tmp2(3)
      pa%adi(:npsit) = tmp2(:npsit)

      pa%chie(:npsit) = chiesec(:npsit)
      pa%chii(:npsit) = chiisec(:npsit)

!     addition to 2.019
      pa%q95 = global(61)

      pa%pohm(1:npsit)=0.0d0
      do j=2, npsit
      term = (.5*(as(j-1)+as(j)) +vlooph(j) )/tpi                           &
     &     *(gxmja2(j)-gxmja2(j-1))
      if(term .lt. 0.0d0) term=0.0
      pa%pohm(j)=pa%pohm(j-1)+term
      enddo
      pa%pohm(1:npsit)=pa%pohm(1:npsit)*udsp/udst
!     write(6,*) " pohm : ", global(47), global(167)

      pa%qie(1:npsit)=0.0d0
      do j=2, npsit
      dvol=vary(j)-vary(j-1)
      term =                                                                &
     &  equila(j)*(-(1.+avez(j))*ade(j)+avez(j)*adp(j))/vpg(j)
      term = 1.5d0*term*dvol
      pa%qie(j)=pa%qie(j-1)+term
      enddo
      pa%qie(1:npsit)=pa%qie(1:npsit)*udsp/udst

      pa%curbs(1:npsit)=0.0d0
      do j=2, npsit
      term = 0.5d0*(ajavbs(j)+ajavbs(j-1))*udsi
      fac  = 0.5d0*(ajbtsq(j)/ajbsq(j) + ajbtsq(j-1)/ajbsq(j-1))
      pa%curbs(j)=pa%curbs(j-1)+term*fac*dpsi
!     second way, but is used in curdrive when reading from ps
!     call geval(xsv(j),2,gval,gpval,gppval,imag,jmag)
!     fac = xmja(j)*gval*(xsv2(j)-xsv2(j-1))
!     fac = fac/tpi
!     term = 0.5d0*(ajavbs(j)+ajavbs(j-1))*udsi
!     pa%curbs(j)=pa%curbs(j-1)+term*fac
      enddo
      pa%curoh(1:npsit)=0.0d0
      do j=2, npsit
      ajtzon = (gxmja2(j)-gxmja2(j-1))*rdpsi*udsi/tpi
      fac = .5*(ajbtsq(j)/ajbsq(j) + ajbtsq(j-1)/ajbsq(j-1))
      abszon = .5*(ajavbs(j)+ajavbs(j-1))*udsi
      acdzon = .5*(ajavcd(j)+ajavcd(j-1))*udsi
      afwzon = .5*(ajavfw(j)+ajavfw(j-1))*udsi
      aeczon = .5*(ajavec(j)+ajavec(j-1))*udsi
      alhzon = ajavlh(j)*udsi
      apszon = ajpary(j)
      pa%curoh(j)=pa%curoh(j-1)+dpsi*(ajtzon-fac*      &
     &            (abszon+acdzon+afwzon+alhzon+aeczon))     &
     &            -apszon
      enddo

!cj aug-15_2011 pilh==0, pelh (electron heating power by LH), curlh (LH current drive)
!fmp may 2013 - add IF statement to avoid TSC to overwrite LH profiles when use_tsc_lhh=.F.
      write(6,*) 'I am about to write the LH to the plasma state'
      if (from_xplasma .and. use_tsc_lhh) then
        write(6,*) 'I am writing powtsc and currtsc to the plasma state'
        pa%pelh(:npsit)=powtsc(:npsit) * 1.e6
        pa%curlh(:npsit)=currtsc(:npsit) * 1.e4
      endif
!cj mar14_2012 debug
      powtot=0.
      curtot=0.
      do j=2, npsit
         powtot=powtot+powtsc(j) * 1.e6 * (vary(j)-vary(j-1))
         curtot=curtot+currtsc(j) * 1.e4 * (vary(j)-vary(j-1)) / (tpi * xplas)
         !write(*,*) "mypacurlh=", j, currtsc(j)*1.e4, currtsc(j)*1.e4*(vary(j)-vary(j-1)) / (tpi*xplas)
      enddo
      write(*,*) "mypacurlh powtot curtot=", powtot, curtot

      l_first=.false.

      DEALLOCATE(tmp2)
      DEALLOCATE(te2,ti2,ane2,zeff2,savebre2,savecyc2,saveimp2)

      RETURN
      END

