       MODULE TRTSCJ
       IMPLICIT NONE

       REAL*8 :: tr_ts1, tr_ts2
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_ibeam
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_ifusn
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_iicrf, tr_iecrf
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_ilh
       REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tr_ife, tr_djde
       INTEGER :: first_read_j=1

       
       CONTAINS
       subroutine get_dataJ_from_xplasma(filenamx, nxin )
       USE trtsc
       USE plasma_state_mod
       USE sawtooth_mod
       USE clinam, only : xsv, xsv2, imag, jmag, xmja, tpi, dpsi
!      USE scr3
!      USE saprop
      USE clinam, only : kcycle, ineg, npsit



       implicit none
       INTEGER, intent (in) :: nxin
       REAL*8, ALLOCATABLE, DIMENSION(:) :: x, xh
       REAL*8, ALLOCATABLE, DIMENSION(:) :: pint, pdens
       REAL*8, ALLOCATABLE, DIMENSION(:) :: tmp, tmp2
       REAL*8, ALLOCATABLE, DIMENSION(:) :: b_dot_gradphi
       CHARACTER*32 :: list_name
       CHARACTER*32, ALLOCATABLE, DIMENSION(:) :: enames
       CHARACTER*80, ALLOCATABLE, DIMENSION(:) :: elabels
       REAL*8, ALLOCATABLE, DIMENSION(:) :: rvals
       INTEGER, ALLOCATABLE, DIMENSION(:) :: itype, ivals
       CHARACTER*(*) :: filenamx

       INTEGER :: ix, iwant, id, list_id, inum, istat
       INTEGER :: it, ii, len_ename
       REAL*8 :: con, sum, sum1,sum2,factor1, factor2
       REAL*8 :: gppval, gpval, gval
       logical :: smooth=.true., norm=.true.
       integer :: icount=0


       call ps_get_plasma_state(ierr,filenamx)
       if(ierr.ne.0) then
         write(6,*) " in get_dataJ_from_xplasma : failed to read ",     &
                      trim(filenamx)
         from_xplasma_j=.false.
         from_xplasma=.false.
         ineg=61
         return
       endif

       if(allocated(tr_ibeam)) deallocate(tr_ibeam)
       if(allocated(tr_ifusn)) deallocate(tr_ifusn)
       if(allocated(tr_iicrf)) deallocate(tr_iicrf)
       if(allocated(tr_iecrf)) deallocate(tr_iecrf)
       if(allocated(tr_ilh)) deallocate(tr_ilh)
       if(allocated(tr_ife)) deallocate(tr_ife)
       if(allocated(tr_djde)) deallocate(tr_djde)
       allocate(tr_ibeam(nxin+1,2),stat=ierr)
       allocate(tr_ifusn(nxin+1,2),stat=ierr)
       allocate(tr_iicrf(nxin+1,2),stat=ierr)
       allocate(tr_iecrf(nxin+1,2),stat=ierr)
       allocate(tr_ilh(nxin+1,2),stat=ierr)
       allocate(tr_ife(nxin+1,2),stat=ierr)
       allocate(tr_djde(nxin+1,2),stat=ierr)
       tr_ibeam = 0.0d0
       tr_ifusn = 0.0d0
       tr_iicrf = 0.0d0
       tr_iecrf = 0.0d0
       tr_ilh = 0.0d0
       tr_ife = 0.0d0
       tr_djde = 0.0d0
!
       if(allocated(x)) deallocate(x)
       if(allocated(xh)) deallocate(xh)
!      if(allocated(b_dot_gradphi)) deallocate(b_dot_gradphi)
       if(allocated(pint)) deallocate(pint)
       if(allocated(pdens)) deallocate(pdens)
       allocate(x(nxin),stat=ierr)
       allocate(xh(nxin-1),stat=ierr)
!      allocate(b_dot_gradphi(nxin+1),stat=ierr)
       allocate(pint(nxin-1),stat=ierr)
       allocate(pdens(nxin+1),stat=ierr)

       do ix = 1, nxin
         x(ix) = float(ix-1)/float(nxin-1)
         x(ix) = sqrt(x(ix))
       enddo
       x(1) = 0.D0
       x(nxin) = 1.D0

       tpi=2.0d0*acos(-1.0d0)
       iwant=0

!      sum=0.0
!      do i=1,ps%nrho_fi-1
!      sum=sum+ps%curbeam(i)
!      enddo
!      write(6,*) "beam current from TRANSP: ", sum

       if(.not. use_tsc_beam) then
         id=ps%id_curbeam
         if(saw_tr) id=saw1%id_curbeam
         call ps_rho_rezone1(x,id,pint,ierr,curdens=norm,zonediff=xh)
         if( ierr .ne. 0 ) then
           write(6,*) "---> curbeam read error"  ,id,ierr,kcycle
         else
           write(6,*) "---> curbeam read -ok- ",id,kcycle,npsit
         endif
!
!........compute current in each zone and sum
         sum=0.0
         do ix=1,nxin-1
           pdens(ix+1) =  pint(ix)*xh(ix)
           sum=sum+pdens(ix+1)
         enddo
         write(6,*) "beam current to TSC: ", sum
!  
!........convert from current in zone to current per unit toroidal flux

         do ix=2,nxin
           pdens(ix)=pdens(ix)/dpsi
         enddo
         pdens(1) = pdens(2)
         pdens(nxin+1) = 2.*pdens(nxin)-pdens(nxin-1)
         do ix = 1, nxin+1
           tr_ibeam(ix,1)=pdens(ix)
           tr_ibeam(ix,2)=pdens(ix)
         enddo
!
!        debug printout
         write(6,4010)
         do ix=1,nxin
         factor2  = 1./dpsi
         write(6,4011) ix,pdens(ix),pint(ix),xh(ix),factor2
         enddo
 4010    format(" ix          pdens       pint        xh       factor   &
     &                1/dpsi")
 4011    format(i5,1p6e12.4)
  
       endif

       if(.not. use_tsc_fp) then
         id=ps%id_curfusn
         if(saw_tr) id=saw1%id_curfusn
         call ps_rho_rezone1(x,id,pint,ierr,                               &
     &                     curdens=norm,zonediff=xh)
         if( ierr .ne. 0 ) then
           write(6,*) "---> curfusn read error"
           from_xplasma_j = .false.
           pint = 0.0d0
         endif
         do ix=1,nxin-1
           pdens(ix+1) =  pint(ix)*xh(ix)
         enddo
         do ix=2,nxin
           call geval(xsv(ix),2,gval,gpval,gppval,imag,jmag)
           pdens(ix)=pdens(ix)*tpi/xmja(ix)/gval/(xsv2(ix)-xsv2(ix-1))
         enddo
         pdens(1) = pdens(2)
         pdens(nxin+1) = 2.*pdens(nxin)-pdens(nxin-1)
         do ix = 1, nxin+1
           tr_ifusn(ix,1)=pdens(ix)
           tr_ifusn(ix,2)=pdens(ix)
         enddo
       endif

       if(.not. use_tsc_ich) then
         id=ps%id_curich
         if(saw_tr) id=saw1%id_curich
         call ps_rho_rezone1(x,id,pint,ierr,                               &
     &                     curdens=norm,zonediff=xh)
         if( ierr .ne. 0 ) then
           write(6,*) "---> curich read error"
           from_xplasma_j = .false.
           pint = 0.0d0
         endif
         do ix=1,nxin-1
           pdens(ix+1) =  pint(ix)*xh(ix)
         enddo
         do ix=2,nxin
           call geval(xsv(ix),2,gval,gpval,gppval,imag,jmag)
           pdens(ix)=pdens(ix)*tpi/xmja(ix)/gval/(xsv2(ix)-xsv2(ix-1))
         enddo
         pdens(1) = pdens(2)
         pdens(nxin+1) = 2.*pdens(nxin)-pdens(nxin-1)
         do ix = 1, nxin+1
           tr_iicrf(ix,1)=pdens(ix)
           tr_iicrf(ix,2)=pdens(ix)
         enddo
         sum=0.0
         do ix=1,nxin-1
           pdens(ix+1) =  pint(ix)*xh(ix)
           sum=sum+pdens(ix+1)
         enddo
         write(6,*) "IC current to TSC: ", sum
       endif

       if(.not. use_tsc_ech) then
         id=ps%id_curech
         if(saw_tr) id=saw1%id_curech
         call ps_rho_rezone1(x,id,pint,ierr,                               &
     &                     curdens=norm,zonediff=xh)
         if( ierr .ne. 0 ) then
           write(6,*) "---> curech read error"
           from_xplasma_j = .false.
           pint = 0.0d0
         endif
         do ix=1,nxin-1
           pdens(ix+1) =  pint(ix)*xh(ix)
         enddo
         do ix=2,nxin
           call geval(xsv(ix),2,gval,gpval,gppval,imag,jmag)
           pdens(ix)=pdens(ix)*tpi/xmja(ix)/gval/(xsv2(ix)-xsv2(ix-1))
         enddo
         pdens(1) = pdens(2)
         pdens(nxin+1) = 2.*pdens(nxin)-pdens(nxin-1)
         do ix = 1, nxin+1
           tr_iecrf(ix,1)=pdens(ix)
           tr_iecrf(ix,2)=pdens(ix)
         enddo
!........compute current in each zone and sum
         sum=0.0
         do ix=1,nxin-1
           pdens(ix+1) =  pint(ix)*xh(ix)
           sum=sum+pdens(ix+1)
         enddo
         write(6,*) "EC current to TSC: ", sum
       endif

       if(.not. use_tsc_lhh) then
         if(allocated(ps%curlh)) then
           write(6,*) 'reading LHCD from plasma state'
           id=ps%id_curlh
           if(saw_tr) id=saw1%id_curlh
           call ps_rho_rezone1(x,id,pint,ierr,                                &
                             curdens=norm,zonediff=xh)
           if( ierr .ne. 0 ) then
             write(6,*) "---> curlh read error"
             from_xplasma_j = .false.
             pint = 0.0d0
           endif
           do ix=1,nxin-1
             pdens(ix+1) =  pint(ix)*xh(ix)
           enddo
           do ix=2,nxin
             call geval(xsv(ix),2,gval,gpval,gppval,imag,jmag)
             pdens(ix)=pdens(ix)*tpi/xmja(ix)/gval/(xsv2(ix)-xsv2(ix-1))
           enddo
           pdens(1) = pdens(2)
           pdens(nxin+1) = 2.*pdens(nxin)-pdens(nxin-1)
           do ix = 1, nxin+1
             tr_ilh(ix,1)=pdens(ix)
             tr_ilh(ix,2)=pdens(ix)
           enddo
         else
           write(6,*) "---> curlh not exists in ps use_tsc_lhh=", use_tsc_lhh
         endif
         sum1=0.0
         do ix=1, ps%nrho_lhrf-1
           sum1=sum1+ps%curlh(ix)
         enddo
         sum2=0.0
         do ix=1,nxin-1
           pdens(ix+1) =  pint(ix)*xh(ix)
           sum2=sum2+pdens(ix+1)
         enddo
         write(6,*) "LH current to TSC: ", sum1,sum2
      endif

      if(ierr .ne. 0 .or. (.not. from_xplasma_j)) then
        tr_ibeam = 0.0d0
        tr_ifusn = 0.0d0
        tr_iicrf = 0.0d0
        tr_iecrf = 0.0d0
        tr_ilh = 0.0d0
        from_xplasma_j = .false.
!        from_xplasma = .false.
      endif
      return

       end subroutine get_dataJ_from_xplasma


       end MODULE TRTSCJ

      subroutine curdrive
!
!.....define quantities needed for current drive
!.....(neutral beam and bootstrap current)
!
      USE trtsc
      USE trtscj
      USE CLINAM
      USE SAPROP
      USE SCR3
!     USE nclass_tsc
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ifircd,ios2,ios11,l,lsav,j,iskip10,iclh,iraytd,ns
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pval,ppval,eval,epval,rval,rpval
      REAL*8 rm2avx,epsx,tempxe,tempxpe,pvali,pvalimin,ppvali
      REAL*8 tempxi,tempxpi,ajavbsn,factor,zbeam,vbeam,alam0,alam1
      REAL*8 anz,amz,coef1,exp1,coef2,exp2,alphz,exp3,vte,vbbar,xf
      REAL*8 fneo,fbeam,tause,pfac,denom,ffact,acfw,dcfw,a1cfw
      REAL*8 a2cfw,ar,f1,f2,f3,gval,gpval,gppval,fp2,btmax,fcyc
      REAL*8 fcyc2,fcyci,fcyci2,fpi2,y2x,y2rat,accn,anecm,tpsh
      REAL*8 radpos,w2lim,zdiff,tfs,ane13,aloglam,aclh,dclh,a1clh
      REAL*8 a2clh,fac,bpowsmo,bpowsmn,diff,anorm,powsum,cursum
      REAL*8 ajavlhm,aterm,ajnew,cursum2,rathlh,backa,fudgefac
      REAL*8 sum,dum
!============
      data ifircd/0/
      character*120 :: filenamx
      integer :: npsit_old
      save :: npsit_old
!
!     dimension wrat(ppsi),anucoll(ppsi),tauconf(ppsi),disedge(ppsi)
!     dimension floss(ppsi),tauf(ppsi),form(ppsi),ajavlho(ppsi)
!     dimension djdets2(ppsi),djdelhn(ppsi)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: wrat
      REAL*8, ALLOCATABLE, DIMENSION(:) :: anucoll
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tauconf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: disedge
      REAL*8, ALLOCATABLE, DIMENSION(:) :: floss
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tauf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: form
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ajavlho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: djdets2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: djdelhn
!============      
      logical :: first_read_beam=.true., first_read_fw=.true.
      logical :: first_read_lhh=.true.
      logical :: first_read_ec=.true.
      logical :: ex_beam=.false., ex_lhh=.false., ex_fw=.false.
      logical :: ex_ec=.false.
      real*8 :: tint, rint, xjb1, xjb2, xjf1, xjf2, xjl1, xjl2,xje1,xje2
      real*8, dimension(:),allocatable :: time_beam, rad_beam
      real*8, dimension(:),allocatable :: time_lhh, rad_lhh
      real*8, dimension(:),allocatable :: time_fw, rad_fw
      real*8, dimension(:),allocatable :: time_ec, rad_ec
      real*8, dimension(:,:),allocatable :: xjb, xjf, xjl, xje
      real*8, dimension(:),allocatable :: xib, xif, xil, curden, xie
      integer :: j1, j2, kk, ii
      integer :: ntime_beam, nrad_beam
      integer :: ntime_lhh, nrad_lhh
      integer :: ntime_fw, nrad_fw
      integer :: ntime_ec, nrad_ec
      integer :: ibootst_temp
      real*8 :: tfluxd, curr, sumc 
!============      
      IF(.not.ALLOCATED(wrat)) ALLOCATE( wrat(ppsi), STAT=istat)
      IF(.not.ALLOCATED(anucoll)) ALLOCATE( anucoll(ppsi), STAT=istat)
      IF(.not.ALLOCATED(tauconf)) ALLOCATE( tauconf(ppsi), STAT=istat)
      IF(.not.ALLOCATED(disedge)) ALLOCATE( disedge(ppsi), STAT=istat)
      IF(.not.ALLOCATED(floss)) ALLOCATE( floss(ppsi), STAT=istat)
      IF(.not.ALLOCATED(tauf)) ALLOCATE( tauf(ppsi), STAT=istat)
      IF(.not.ALLOCATED(form)) ALLOCATE( form(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ajavlho)) ALLOCATE( ajavlho(ppsi), STAT=istat)
      IF(.not.ALLOCATED(djdets2)) ALLOCATE( djdets2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(djdelhn)) ALLOCATE( djdelhn(ppsi), STAT=istat)
!fmp  if( use_user_beam .or. use_user_fw .or. use_user_lhh .or. use_tsc_ech) then
      if( use_user_beam .or. use_user_fw .or. use_user_lhh .or. use_user_ech) then
      IF(.not.ALLOCATED(curden)) ALLOCATE( curden(ppsi), STAT=istat )
      endif
!============      
      if (istat .ne. 0) stop 'Allocation Error : curdrive  ' 
!============      
      if(.not. from_xplasma) from_xplasma_j=.false.
      if(from_xplasma) then
        filenamx=trim(xplasma_filename)
        if( .not. use_tsc_fp .or. .not. use_tsc_beam .or.                     &
     &      .not. use_tsc_ich .or. .not. use_tsc_lhh .or.                     &
     &      .not. use_tsc_ech ) then

          if(read_xplasma_j .or. (npsit .ne. npsit_old)) then
            call get_dataJ_from_xplasma(filenamx,npsit)
            npsit_old=npsit
            read_xplasma_j = .false.
          endif
        endif
      endif

!
      ifircd=0
      if(kcycle.le.0) ifstrt=0
      if(ilhcd.eq.0 .or. ifk.ne.2 .or. ifircd .ne. 0) go to 49
!     iouttsc(1:6) = 'tscout'
!     iouttsc(7:7) = isuffix(1:1)
      if (numargs .lt. 1) then
         filename = 'tscout' // isuffix(1:1)
      else
         filename = 'tscout'//'.'//trim(suffix)
      end if
      iouttsc = trim(filename)
      open(nlhcd2,file=trim(iouttsc),status='unknown',iostat=ios2)
!     ilhcdou(1:6) = 'lhcdou'
!     ilhcdou(7:7) = isuffix(1:1)
      if (numargs .lt. 1) then
         filename = 'lhcdou' // isuffix(1:1)
      else
         filename = 'lhcdou' // '.' //trim(suffix)
      end if
      ilhcdou = trim(filename)
      open(nlhcd,file=trim(ilhcdou),status='unknown',iostat=ios11)
      ifircd = 1
   49 continue
!
      do 50 l=2,npsit
      ajavbs(l) = 0._R8
   50 continue
      if(ibootst.eq.0) go to 101
      ibootst_temp = ibootst
      if(ibootst .eq. 8) ibootst = 3
!
!.(1).first calculate ajavbs(l) ... bootstrap current density
!                  Reference:   S.P.Hirshman,PF,Vol 31,p. 3150
!                  Reference:   G.R.Harris,EUR-CEA-FC-1436
!
      do 100 l=2,npsit
      call peval(xsv2(l),2,pval,ppval,imag,jmag)
      call eeval(xsv2(l),2,eval,epval,imag,jmag)
      call reval(xsv2(l),idens,1,rval,rpval,imag,jmag)
      rm2avx = xmja2(l)/(tpi*vp2(l)*qprof2(l))
      epsx = sqrt(vary(l)/(2.0_R8*pi*pi*xmag))/xmag
      tempxe = (eval/rval)*udsh
      tempxpe= (epval*rval - rpval*eval)*udsh/rval**2
      pvali = pval - eval
      pvalimin = fracn0*r0*smallt*(acoef(882)-1._R8)
      if(pvali.le.pvalimin) pvali = pvalimin
      ppvali = ppval - epval
      tempxi = (pvali/rval)*udsh
      tempxpi= (ppvali*rval - rpval*pvali)*udsh/rval**2
      if(tempxi .lt.0._R8 ) then
      write(nterm,6001) l,npsit,tempxi,pval,eval,pvali,rval
 6001 format(2i5,1p5e12.4)
      endif
!
!
      call bootstrap(ajavbsn,ibootst, ierr, nout,                        &  
     & zeff,rm2avx, epsx, xmag, qprof2(l), ftrap(l),                     &  
     & eval, epval, tempxe, tempxpe,                                     &  
     & pvali,ppvali,tempxi, tempxpi )
      if(ierr .ne. 0) ineg=48
!
!.....compute <B dot grad phi> / <B dot B> from info available
      factor = acoef(103)*1.E8_R8*tpi / (vp2(l)*bsqar(l))
!
!     alpha = <Jboot dot B> / < B dot B>
      ajavbs(l) = factor * ajavbsn
  100 continue
      ajavbs(1) = ajavbs(2)
!     ajavbs(1) = 0.

!     ibootst = ibootst_temp
!     if(ibootst .eq. 8) then
!     call tsc_nclass
!     ajavbs(2:npsit) = nclass_bs(2:npsit)
!     ajavbs(1) = ajavbs(2)
!     endif

  101 continue
      if(acoef(296).eq.5._R8) return

      do 102 l=2,npsit
      ajavcd(l) = 0._R8
  102 continue

      if(from_xplasma_j .and. .not. use_tsc_beam) then
        fac=usdi
        dum = (times-tr_ts1)/(tr_ts2-tr_ts1)
        dum = 0.0
        do j=1,npsit
          ajavcd(j)=tr_ibeam(j,1)+dum*(tr_ibeam(j,2)-tr_ibeam(j,1))
          ajavcd(j)=fac*ajavcd(j)
        enddo
        goto 299
      endif

      if(use_user_beam .and. first_read_beam) then
        inquire(file="user_beam_current", exist=ex_beam)
        if (ex_beam) then
          open(55,file="user_beam_current",form="formatted",status="old")
          read(55,*) ntime_beam, nrad_beam 
          allocate(time_beam(ntime_beam), rad_beam(nrad_beam))
          allocate(xib(ntime_beam),xjb(nrad_beam,ntime_beam))
          read(55,*) time_beam(1:ntime_beam) 
          read(55,*) rad_beam(1:nrad_beam)
          read(55,*) xib(1:ntime_beam) 
          read(55,*) xjb(1:nrad_beam,1:ntime_beam)
          close (55) 
        endif
        first_read_beam=.false.
      endif

      if(use_user_beam .and. ex_beam) then
        do  kk = 1,ntime_beam-1
          if(times .gt. time_beam(kk) .and. times .le. time_beam(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_beam(j1))/(time_beam(j2)-time_beam(j1))
            exit
          endif
        enddo
        if(times .le. time_beam(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif
        if(times .gt. time_beam(ntime_beam)) then
          j1 = ntime_beam-1
          j2 = ntime_beam
          tint = 1.
        endif
 
        sumc = 0.
!fmp  curr = xib(j1)+(xib(j2)-xib(j1))*tint*usdp/usdt
        curr = xib(j1)+(xib(j2)-xib(j1))*tint
        if (xib(j1).eq.0.0) curr = 0.0_R8
        do j = 2, npsit
          tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
          do  ii=1,nrad_beam-1
            if(tfluxd.gt.rad_beam(ii) .and. tfluxd.le.rad_beam(ii+1)) then
              rint = (tfluxd-rad_beam(ii))/(rad_beam(ii+1)-rad_beam(ii))
              xjb1 = xjb(ii,j1) + (xjb(ii+1,j1)-xjb(ii,j1))*rint
              xjb2 = xjb(ii,j2) + (xjb(ii+1,j2)-xjb(ii,j2))*rint
!fmp          curden(j) = (xjb1 + (xjb2-xjb1)*tint)*usdp/usdt
              curden(j) = xjb1 + (xjb2-xjb1)*tint
              call geval(xsv(j),2,gval,gpval,gppval,imag,jmag)
              sumc = sumc+curden(j)*(xsv2(j)-xsv2(j-1))*xmja(j)*gval/tpi
              exit
            endif
          enddo
          if(tfluxd .le. rad_beam(1)) then
!fmp        curden(j) = (xjb(1,j1) + (xjb(1,j2)-xjb(1,j1))*tint)*usdp/usdt
            curden(j) = xjb(1,j1) + (xjb(1,j2)-xjb(1,j1))*tint
          endif
        enddo
        do j=2,npsit
          if(sumc .eq. 0.) then
            ajavcd(j) = 0.
          else
            ajavcd(j) = curr*curden(j)*usdi/sumc
          endif
        enddo
        ajavcd(1) = ajavcd(2)
        goto 299         !fmp 02/26/2012
      endif

!     if(.not. from_xplasma_j .or. use_tsc_beam .or.                    &
!    &    (use_user_beam .and. .not. ex_beam))  then
      if(fracpar.eq.0._R8) go to 201
!
!.(2).calculate ajavcd(l) ... beam driven current
!
!
!.....reference:  S.P.Hirshman,PF,Vol.23,p.1238(1980)
!
      zbeam = 1.0_R8
      vbeam = sqrt((2._R8*ebeamkev*1.6022E-16_R8)/(ambeam*1.6726E-27_R8)  &  
     & )
      alam0 = (1._R8+ 0.374_R8*zeff)*(1.085_R8+ 2.783_R8*zeff +          &  
     & 0.424_R8*zeff**2)                                                 &  
     & /(zeff*(1._R8+ 0.292_R8*zeff)*(1.0_R8+ 1.161_R8*zeff + 0.16_R8*   &  
     & zeff**2))
      alam1 = (0.877_R8+ 0.695_R8*zeff + 0.095_R8*zeff**2)               &  
     & /((1._R8+0.292_R8*zeff)*(1.0_R8+ 1.161_R8*zeff + 0.16_R8*zeff**2)  &  
     &  )
      anz = 1.65_R8+ 0.22_R8/(zeff-0.55_R8)
      amz = 3.37_R8- 0.85_R8*(zeff+0.48_R8)
      if(amz.lt.2.0_R8) amz = 2.0_R8
      coef1 = 0.75_R8*sqrt(pi)*zeff*(alam0+3.0_R8*alam1)
      exp1 = 1._R8/(3._R8+2._R8*anz)
      coef2 = 5._R8*anz/6._R8
      exp2 = anz/(3._R8+2._R8*anz)
      alphz = coef1**exp1*coef2**exp2
      exp3 = (3._R8+2._R8*anz)/amz
!
      do 200 l=2,npsit
!
!.....electron thermal velocity (m/sec)
      vte = 4.19E5_R8*sqrt(te(l))
      vbbar = vbeam/vte
!
!.....neoclassical correction to reverse electron current (06/07/2003)
      xf = ftrap(l)/(1._R8-ftrap(l))
      fneo = xf*((0.754_R8+2.21_R8*zeff+zeff**2)+xf*(0.348_R8+1.243_R8*  &  
     & zeff                                                              &  
     & +zeff**2))/(1.414_R8*zeff+zeff**2+xf*(0.754_R8+2.657_R8*zeff+     &  
     & 2._R8*zeff**2)+                                                   &  
     & xf**2*(0.348_R8+1.243_R8*zeff+zeff**2))
 
!.....fraction of beam ion current not cancelled by plasma
      fbeam = 1._R8-((zbeam/zeff)*(1._R8+1.2_R8*(vbbar**2/anz))**anz     &  
     &      / (1._R8+(vbbar/alphz)**amz)**exp3)*(1._R8-fneo)
!
!.....slowing down time(dimensionless units)
      tause = (1.17E18_R8*ambeam*(te(l)*.001_R8)**1.5_R8/ane(l))*usdt
!
!.....proportionality factor
      pfac = tause*(1.6E-19_R8)/sqrt(2._R8*1.672E-27_R8*ambeam*ebeamkev*  &  
     & 1.6E-16_R8)
!
!.....dimensionless current density
      ajavcd(l) = pfac*(savee(l)+savei(l))*fbeam*fracpar*xplas/gzero
!
  200 continue
      ajavcd(1) = ajavcd(2)
  201 continue
!     endif
!
  299 continue      !fmp 02/26/2012
!.... calculate ajavfw(l) ... fast wave current
!
      if(ifwcd .le. 0) go to 699
      if(from_xplasma_j .and. .not. use_tsc_ich) then
        fac = usdi
        dum = (times-tr_ts1)/(tr_ts2-tr_ts1)
        dum = 0.0
        do j=1,npsit
          ajavfw(j)=tr_iicrf(j,1)+dum*(tr_iicrf(j,2)-tr_iicrf(j,1))
!         ajavfw(j)=tr_iecrf(j,1)+dum*(tr_iecrf(j,2)-tr_iecrf(j,1)) +       &
!    &          ajavfw(j)
          ajavfw(j)=fac*ajavfw(j)
        enddo
        goto 699
      endif

      if(use_user_fw .and. first_read_fw) then
         inquire(file="user_fw_current", exist=ex_fw)
         if (ex_fw) then
           open(55,file="user_fw_current",form="formatted",status="old")
           read(55,*) ntime_fw, nrad_fw 
           allocate(time_fw(ntime_fw), rad_fw(nrad_fw))
           allocate(xif(ntime_fw),xjf(nrad_fw,ntime_fw))
           read(55,*) time_fw(1:ntime_fw) 
           read(55,*) rad_fw(1:nrad_fw)
           read(55,*) xif(1:ntime_fw) 
           read(55,*) xjf(1:nrad_fw,1:ntime_fw)
           close (55) 
         endif
         first_read_fw=.false.
      endif

      if(use_user_fw .and. ex_fw) then
        do  kk = 1,ntime_fw-1
          if(times .gt. time_fw(kk) .and. times .le. time_fw(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_fw(j1))/(time_fw(j2)-time_fw(j1))
            exit
          endif
        enddo
        if(times .le. time_fw(1)) then
           j1 = 1
           j2 = 2
           tint = 0.
        endif
        if(times .gt. time_fw(ntime_fw)) then
           j1 = ntime_fw-1
           j2 = ntime_fw
           tint = 1.
        endif
 
        sumc = 0.
!fmp    curr = xif(j1)+(xif(j2)-xif(j1))*tint*usdp/usdt
        curr = xif(j1)+(xif(j2)-xif(j1))*tint
        if (xif(j1).eq.0.0) curr = 0.0_R8
        do j = 2, npsit
           tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
           do  ii=1,nrad_fw-1
              if(tfluxd .gt. rad_fw(ii) .and. tfluxd .le. rad_fw(ii+1)) then
                 rint = (tfluxd-rad_fw(ii))/(rad_fw(ii+1)-rad_fw(ii))
                 xjf1 = xjf(ii,j1) + (xjf(ii+1,j1)-xjf(ii,j1))*rint
                 xjf2 = xjf(ii,j2) + (xjf(ii+1,j2)-xjf(ii,j2))*rint
!fmp             curden(j) = (xjf1 + (xjf2-xjf1)*tint)*usdp/usdt
                 curden(j) = xjf1 + (xjf2-xjf1)*tint
                 call geval(xsv(j),2,gval,gpval,gppval,imag,jmag)
                 sumc = sumc + curden(j)*(xsv2(j)-xsv2(j-1))*xmja(j)*gval/tpi
                 exit
               endif
           enddo
           if(tfluxd .le. rad_fw(1)) then
!fmp         curden(j) = (xjf(1,j1) + (xjf(1,j2)-xjf(1,j1))*tint)*usdp/usdt
             curden(j) = xjf(1,j1) + (xjf(1,j2)-xjf(1,j1))*tint
           endif
         enddo
         do j=2,npsit
            if(sumc .eq. 0.) then
               ajavfw(j) = 0.
            else
               ajavfw(j) = curr*curden(j)*usdi/sumc
            endif
         enddo
         ajavfw(1) = ajavfw(2)
         goto 699          !fmp 02/26/2012
      endif

!fmp 2013/05/22    if(.not. from_xplasma_j .or. use_tsc_ich .or.                     &
!    &    (use_user_fw .and. .not. ex_fw)) then
!     if( (use_tsc_ich .and. .not. use_user_fw) .or.                    &
!    &    (use_user_fw .and. .not. ex_fw)) then
      do 490 l=1,ntpts-1
        lsav = l
        if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 497
  490 continue
  497 continue
!.....linear time point interpolation
      denom = tpro(lsav+1)-tpro(lsav)
      if(denom.eq.0) denom = 1._R8
      ffact = (time-tpro(lsav))/denom
      acfw = acfwd(lsav) + ffact * (acfwd(lsav+1) - acfwd(lsav))
      dcfw = dcfwd(lsav) + ffact * (dcfwd(lsav+1) - dcfwd(lsav))
      a1cfw = a1cfwd(lsav) + ffact * (a1cfwd(lsav+1) - a1cfwd(lsav))
      a2cfw = a2cfwd(lsav) + ffact * (a2cfwd(lsav+1) - a2cfwd(lsav))
      sum = 0._R8
      do 680 l=2,npsit
        ar = (xsv(l)-psimin)/(psilim-psimin)
        if(abs(ar).ge.1.0_R8) go to 698
        if(ar.lt.0._R8) ar = 0._R8
        f1 = dcfw**2/((ar-acfw)**2 + dcfw**2)
        f2 = 1._R8
        if(a1cfw .gt. 0._R8) f2 = ar**a1cfw
        f3 = 1._R8
        if(a2cfw .gt. 0._R8) f3 = (1._R8-ar)**a2cfw
        form(l) = f1*f2*f3
        go to 679
  698   form(l) = 0._R8
  679   continue
        call geval(xsv(l),2,gval,gpval,gppval,imag,jmag)
        sum = sum + form(l)*(xsv2(l)-xsv2(l-1))*xmja(l)*gval / tpi
  680 continue
      do 681 j=2,npsit
        ajavfw(j) = fwcd(lsav)*usdi*form(j)/sum
  681 continue
      ajavfw(1) = ajavfw(2)
!     endif

  699 continue 
!cj...  699 continue
!cj...      do 289 l=2,npsit
!cj...      IF ( chiesec(l) .LE. 1.E-6_R8) chiesec(l) = 1.E-6_R8
!cj...  289 continue
!cj...!
!cj...      if(ilhcd .le. 0) go to 701
!cj...      if(ifk.gt.0) go to 650
!cj...
!cj...      if(from_xplasma_j .and. .not. use_tsc_lhh) then
!cj...      fac = usdi
!cj...      dum = (times-tr_ts1)/(tr_ts2-tr_ts1)
!cj...      dum = 0.0
!cj...      do j=1,npsit+1
!cj...      ajavlh(j)=tr_ilh(j,1)+dum*(tr_ilh(j,2)-tr_ilh(j,1))
!cj...      ajavlh(j)=fac*ajavlh(j)
!cj...      enddo
!cj...      go to 725
!cj...      endif

!
!.... calculate ajavec(l) ... eccd current
!
!      if(ieccd .le. 0) go to 699
      if(ieccd .le. 0) go to 799
!fmp      if(from_xplasma_j .and. .not. use_tsc_ich) then
      if(from_xplasma_j .and. .not. use_tsc_ech) then
        fac = usdi
        dum = (times-tr_ts1)/(tr_ts2-tr_ts1)
        dum = 0.0
        do j=1,npsit    
!     ajavec(j)=tr_iicrf(j,1)+dum*(tr_iicrf(j,2)-tr_iicrf(j,1))
!     ajavec(j)=tr_iecrf(j,1)+dum*(tr_iecrf(j,2)-tr_iecrf(j,1)) +       &
!    &          ajavec(j)
          ajavec(j)=tr_iecrf(j,1)+dum*(tr_iecrf(j,2)-tr_iecrf(j,1))
          ajavec(j)=fac*ajavec(j)
        enddo
        goto 799  !fmp 05/23/2013
      endif

!fmp  if( use_tsc_ech .and. first_read_ec) then
      if( use_user_ech .and. first_read_ec) then
        inquire(file="user_ec_current", exist=ex_ec)
        if (ex_ec) then
          open(55,file="user_ec_current",form="formatted",status="old")
          read(55,*) ntime_ec, nrad_ec 
          allocate(time_ec(ntime_ec), rad_ec(nrad_ec))
          allocate(xie(ntime_ec),xje(nrad_ec,ntime_ec))
          read(55,*) time_ec(1:ntime_ec) 
          read(55,*) rad_ec(1:nrad_ec)
          read(55,*) xie(1:ntime_ec) 
          read(55,*) xje(1:nrad_ec,1:ntime_ec)
          close (55) 
        endif
        first_read_ec=.false.
      endif

!fmp  if( use_tsc_ech .and. ex_ec) then
      if( use_user_ech .and. ex_ec) then
        do  kk = 1,ntime_ec-1
          if(times .gt. time_ec(kk) .and. times .le. time_ec(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_ec(j1))/(time_ec(j2)-time_ec(j1))
            exit
          endif
        enddo
        if(times .le. time_ec(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif
        if(times .gt. time_ec(ntime_ec)) then
          j1 = ntime_ec-1
          j2 = ntime_ec
          tint = 1.
        endif
 
        sumc = 0.
!fmp  curr = xie(j1)+(xie(j2)-xie(j1))*tint*usdp/usdt
        curr = xie(j1)+(xie(j2)-xie(j1))*tint
        if (xie(j1).eq.0.0) curr = 0._R8
        do j = 2, npsit
          tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
          do  ii=1,nrad_ec-1
            if(tfluxd .gt. rad_ec(ii) .and. tfluxd .le. rad_ec(ii+1)) then
              rint = (tfluxd-rad_ec(ii))/(rad_ec(ii+1)-rad_ec(ii))
              xje1 = xje(ii,j1) + (xje(ii+1,j1)-xje(ii,j1))*rint
              xje2 = xje(ii,j2) + (xje(ii+1,j2)-xje(ii,j2))*rint
!fmp  curden(j) = (xje1 + (xje2-xje1)*tint)*usdp/usdt
              curden(j) = xje1 + (xje2-xje1)*tint
              call geval(xsv(j),2,gval,gpval,gppval,imag,jmag)
              sumc = sumc+curden(j)*(xsv2(j)-xsv2(j-1))*xmja(j)*gval/tpi
              exit
            endif
          enddo
          if(tfluxd .le. rad_ec(1)) then
!fmp  curden(j) = (xje(1,j1) + (xje(1,j2)-xje(1,j1))*tint)*usdp/usdt
            curden(j) = xje(1,j1) + (xje(1,j2)-xje(1,j1))*tint
          endif
        enddo
        do j=2,npsit
          if(sumc .eq. 0.) then
            ajavec(j) = 0.
          else
            ajavec(j) = curr*curden(j)*usdi/sumc
          endif
        enddo
        ajavec(1) = ajavec(2)
        GOTO 799            !fmp 02/26/2012
      endif


!fmp  if(.not. from_xplasma_j .or. use_tsc_ich .or.                     &
!fmp &    (use_tsc_ech .and. .not. ex_ec)) then
!     if( (use_tsc_ech .and. .not. use_user_ech) .or.               &
!    &    (use_user_ech .and. .not. ex_ec)) then
      do 590 l=1,ntpts-1
        lsav = l
        if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 597
  590 continue
  597 continue
!.....linear time point interpolation
      denom = tpro(lsav+1)-tpro(lsav)
      if(denom.eq.0) denom = 1._R8
      ffact = (time-tpro(lsav))/denom
      acfw = aecd(lsav) + ffact * (aecd(lsav+1) - aecd(lsav))
      dcfw = decd(lsav) + ffact * (decd(lsav+1) - decd(lsav))
      a1cfw = a1ecd(lsav) + ffact * (a1ecd(lsav+1) - a1ecd(lsav))
      a2cfw = a2ecd(lsav) + ffact * (a2ecd(lsav+1) - a2ecd(lsav))
      sum = 0._R8
      do 780 l=2,npsit
        ar = (xsv(l)-psimin)/(psilim-psimin)
        if(abs(ar).ge.1.0_R8) go to 798
        if(ar.lt.0._R8) ar = 0._R8
        f1 = dcfw**2/((ar-acfw)**2 + dcfw**2)
        f2 = 1._R8
        if(a1cfw .gt. 0._R8) f2 = ar**a1cfw
        f3 = 1._R8
        if(a2cfw .gt. 0._R8) f3 = (1._R8-ar)**a2cfw
        form(l) = f1*f2*f3
        go to 779
  798   form(l) = 0._R8
  779   continue
        call geval(xsv(l),2,gval,gpval,gppval,imag,jmag)
        sum = sum + form(l)*(xsv2(l)-xsv2(l-1))*xmja(l)*gval / tpi
  780 continue
      do 781 j=2,npsit
        ajavec(j) = eccd(lsav)*usdi*form(j)/sum
  781 continue
      ajavec(1) = ajavec(2)
!     endif

! 699 continue      
  799 continue
      do 289 l=2,npsit
      IF ( chiesec(l) .LE. 1.E-6_R8) chiesec(l) = 1.E-6_R8
  289 continue
!
      if(ilhcd .lt. 1) go to 701
      if(ifk.gt.0) go to 650

      if(from_xplasma_j .and. .not. use_tsc_lhh) then
        fac = usdi
        dum = (times-tr_ts1)/(tr_ts2-tr_ts1)
        dum = 0.0
        do j=1,npsit+1
          ajavlh(j)=tr_ilh(j,1)+dum*(tr_ilh(j,2)-tr_ilh(j,1))
          ajavlh(j)=fac*ajavlh(j)
        enddo
        go to 725
      endif

      if(use_user_lhh .and. first_read_lhh) then
        inquire(file="user_lhh_current", exist=ex_lhh)
        if (ex_lhh) then
          open(55,file="user_lhh_current",form="formatted",status="old")
          read(55,*) ntime_lhh, nrad_lhh 
          allocate(time_lhh(ntime_lhh), rad_lhh(nrad_lhh))
          allocate(xil(ntime_lhh),xjl(nrad_lhh,ntime_lhh))
          read(55,*) time_lhh(1:ntime_lhh) 
          read(55,*) rad_lhh(1:nrad_lhh)
          read(55,*) xil(1:ntime_lhh) 
          read(55,*) xjl(1:nrad_lhh,1:ntime_lhh)
          close (55) 
        endif
        first_read_lhh=.false.
      endif

      if(use_user_lhh .and. ex_lhh) then
        do  kk = 1,ntime_lhh-1
          if(times .gt. time_lhh(kk) .and. times .le. time_lhh(kk+1)) then
            j1 = kk
            j2 = kk + 1
            tint = (times-time_lhh(j1))/(time_lhh(j2)-time_lhh(j1))
            exit
          endif
        enddo
        if(times .le. time_lhh(1)) then
          j1 = 1
          j2 = 2
          tint = 0.
        endif
        if(times .gt. time_lhh(ntime_lhh)) then
          j1 = ntime_lhh-1
          j2 = ntime_lhh
          tint = 1.
        endif
 
        sumc = 0.
!fmp  curr = xil(j1)+(xil(j2)-xil(j1))*tint*usdp/usdt
        curr = xil(j1)+(xil(j2)-xil(j1))*tint
        if (xil(j1).eq.0.0) curr = 0._R8
        do j = 2, npsit
          tfluxd = sqrt((float(j-1)*dpsi)/(float(npsit-1)*dpsi))
          do  ii=1,nrad_lhh-1
            if(tfluxd .gt. rad_lhh(ii) .and. tfluxd .le. rad_lhh(ii+1)) then
              rint = (tfluxd-rad_lhh(ii))/(rad_lhh(ii+1)-rad_lhh(ii))
              xjl1 = xjl(ii,j1) + (xjl(ii+1,j1)-xjl(ii,j1))*rint
              xjl2 = xjl(ii,j2) + (xjl(ii+1,j2)-xjl(ii,j2))*rint
!fmp  curden(j) = (xjl1 + (xjl2-xjl1)*tint)*usdp/usdt
              curden(j) = xjl1 + (xjl2-xjl1)*tint
              call geval(xsv(j),2,gval,gpval,gppval,imag,jmag)
              sumc = sumc+curden(j)*(xsv2(j)-xsv2(j-1))*xmja(j)*gval/tpi
              exit
            endif
          enddo
          if(tfluxd .le. rad_lhh(1)) then
!fmp        curden(j) = (xjl(1,j1) + (xjl(1,j2)-xjl(1,j1))*tint)*usdp/usdt
            curden(j) = xjl(1,j1) + (xjl(1,j2)-xjl(1,j1))*tint
          endif
        enddo
        do j=2,npsit
          if(sumc .eq. 0.) then
            ajavlh(j) = 0.
          else
            ajavlh(j) = curr*curden(j)*usdi/sumc
          endif
        enddo
        ajavlh(1) = ajavlh(2)
        go to 725
      endif

!fmp 2013/05/22     if(.not. from_xplasma_j .or. use_tsc_lhh .or.                     &
!    &   (use_user_lhh .and. .not. ex_lhh) ) then
!fmp 2013/05/22 - replace the above with the line below.
! This is to avoid LHCD calculations if reading from user files or
! using the LH component in the IPS 
      if( (use_tsc_lhh .and. .not. use_user_lhh ) .or.                     &
     &   (use_user_lhh .and. .not. ex_lhh) ) then
!
!.(3).calculate current driven by Lower Hybrid waves
        fp2 = 806.2_R8* r0 / freqlh**2
        btmax = ABS(gs(imag,jmag))/xmag
        fcyc = 28._R8 * btmax / freqlh
        fcyc2 = fcyc * fcyc
        fcyci = fcyc * zion / (aion * 1836._R8)
        fcyci2 = fcyci * fcyci
        fpi2 = fp2 * zion**2 / (aion * 1836._R8)
        y2x = 1.0_R8/ (fcyci * fcyc)
        y2rat = y2x / (1.0_R8- y2x)
        IF ( y2x .GE. 1.0_R8.OR. fpi2 .LT. y2rat ) THEN
           accn = SQRT(fp2)/fcyc + SQRT(1._R8+ fp2/fcyc2 - fpi2)
        ELSE
          accn = SQRT( 1._R8/(1._R8- y2x) )
        ENDIF
        do 300 l=2,npsit
!
!.....calculate confinement param needed in calc of loss factor below
!             ref:  Mynick and Strachan,Phys.Fluids,Vol 24 p695 (1981)
!
        anecm = ane(l)*1.E-6_R8
        anucoll(l) = 3.0E-4_R8*anecm / te(l)**1.5_R8
        tpsh = 1._R8 / (pi * SQRT(2.0_R8* xplas))
        radpos = tpsh * sqrt(vary(l))
        disedge(l) = tpsh * ( SQRT(vary(npsit)) +                          &  
     &     0.5_R8* (SQRT(vary(npsit)) - SQRT(vary(npsit-1))) ) - radpos
        tauconf(l) = 0.5_R8*disedge(l)**2/max(chiesec(l),0.1)
!
!.....calculate w2lim based on lh accessibility condition.
!
!.....Normalize w2lim by v-thermal
        w2lim = 7.15E2_R8/(accn*sqrt(te(l)) )
!
!.....Since raytracing results show that w2lim become no larger than 10
!.....when w2lim is greater than 10 we replace it by 10.  Larger values
!.....of w2lim are likely to correspond to a runaway region.
!
        IF (w2lim .GT. 10.0_R8) w2lim = 10.0_R8
!
!.....Check that w2lim is no less than w1lim + 1.0e-3
        zdiff = ABS(w1lim)+ 1.0E-3_R8
        IF (w2lim .LT. zdiff) w2lim = zdiff
        IF (w1lim .LT. 0._R8) w2lim = -w2lim
!
!.....calculate loss factor F
!            ref: S.C.Luckhardt, MIT-PFC/5A-86-4
!.....Originally the expression below that depends on anpar (the
!.....parallel wave number) which was read in (on card 38) was used
!.....to calculate wrat(l).  This value was then used to compute tfs.
!.....Note that the Luckhardt formular used below assumes Zeff=1.
!.....These equations are below.  The value now used is w2lim (the high
!.....velocity end of the lower hybrid spectrum).  This value of
!.....w2lim is calculated above utilizing the accessibility condition.
!     wrat(l) = 7.15e2/(anpar*sqrt(te(l)) )
!     tfs = tauf(l)/wrat(l)**3
!
        tauf(l) = tauconf(l)*anucoll(l)
        tfs = ABS( tauf(l)/w2lim**3 )
        floss(l) = .75_R8*tfs*(1._R8+3._R8*tfs*(1._R8-(1._R8+1._R8/tfs)*   &  
     &   exp(-2._R8/(3._R8*tfs))))
        if(floss(l).gt.1._R8) floss(l) = 1._R8
!
!.....now calculate current driven using Fish formula
!
        ane13 = anecm/1.E13_R8
        aloglam = 9.03_R8- log ( sqrt(ane13)/te(l))
!
!.....note that 15.35e16 = 16 pi eps0**2 / e**2   The value of w1lim,
!.....given by damping, is estimated and read as input on card 38.
!.....The value of w2lim is calculated below based on accessibility
!.....of the lower hybrid waves.
        ajavlh(l) = acoef(106)*floss(l)*15.35E-3_R8*savelh(l)*te(l)*       &  
     &         w2lim**2 * (1._R8-(w1lim/w2lim)**2)/(2._R8                &  
     &         *log(w2lim/w1lim))                                        &  
     &          /((5._R8+zeff)*aloglam*ane13*udst)*xplas/gzero
  300 continue
      ajavlh(1) = ajavlh(2)
      if(cprof.eq.1.0_R8) go to 401
      go to 725
!
!
 401  continue
      do 400 l=1,ntpts-1
        lsav = l
        if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 407
  400 continue
  407 continue
!.....linear time point interpolation
      denom = tpro(lsav+1)-tpro(lsav)
      if(denom.eq.0) denom = 1._R8
      ffact = (time-tpro(lsav))/denom
      aclh = aclhd(lsav) + ffact * (aclhd(lsav+1) - aclhd(lsav))
      dclh = dclhd(lsav) + ffact * (dclhd(lsav+1) - dclhd(lsav))
      a1clh = a1clhd(lsav) + ffact * (a1clhd(lsav+1) - a1clhd(lsav))
      a2clh = a2clhd(lsav) + ffact * (a2clhd(lsav+1) - a2clhd(lsav))
      fac = 0._R8
      sum = 0._R8
      do 610 l=2,npsit
        ar = (xsv(l)-psimin)/(psilim-psimin)
        if(abs(ar).ge.1.0_R8) go to 608
        if(ar.lt.0._R8) ar = 0._R8
        f1 = dclh**2/((ar-aclh)**2 + dclh**2)
        f2 = 1._R8
        if(a1clh .gt. 0) f2 = ar**a1clh
        f3 = 1._R8
        if(a2clh .gt. 0) f3 = (1._R8-ar)**a2clh
        form(l) = f1*f2*f3
        go to 609
  608   form(l) = 0._R8
  609   continue
!.....added FEB 18 , 2005
        call geval(xsv(l),2,gval,gpval,gppval,imag,jmag)
        sum = sum + form(l)*(xsv2(l)-xsv2(l-1))*xmja(l)*gval / tpi
        fac = fac + ajavlh(l) * (xsv2(l) - xsv2(l-1))*xmja(l)*gval/tpi
  610 continue
      do 611 j=2,npsit
        ajavlh(j) = fac*form(j)/sum
  611 continue
      ajavlh(1) = ajavlh(2)
      endif
      go to 725
!
!.....use values read from disk file tscouta
!                          currtsc is in amps/cm**2
  650 continue
      if(ifk.eq.1) go to 651
!
!...NOTE:   nlhcd is lu which TSC writes equilibrium info
!           nlhcd2 is "  "    "   reads current and power
!           irayt=1 to trace rays,  2 otherwise
!           iplot=1 to make plots,  0 otherwisw
!
      if(iplt.gt.nskipl) iplot=1
      if(mod(kcycle,nslhrt).eq.0)irayt=1
      if(ifstrt.eq.0) then
                      irayt=1
                      iplot=1
                      ifstrt=1
                      endif
      do 600 l=1,ntpts-1
      lsav = l
      if(tpro(l).le.time .and. tpro(l+1).gt.time) go to 607
  600 continue
  607 continue
      bpowsmo = bpowsmw
      bpowsmn = plhamp(lsav)*1.E-6_R8
      iskip10 = 0
      diff = bpowsmn - bpowsmo
      anorm = max(bpowsmo,bpowsmn)
      if(diff.gt. 0.1_R8*anorm) then
        diff = 0.1_R8*anorm
        iskip10 = 1
                             endif
      if(iskip10.eq.1) go to 619
!
!.....call LSC every nslhpc cycles only
      if(mod(kcycle,nslhpc).gt.0) go to 700
  619 continue
      bpowsmw = bpowsmo + diff
      do 652 j=1,npsi
      powtsc(j)=0._R8
      currtsc(j)=0._R8
      djdetsc(j) = 0._R8
      djdelh(j) = 0._R8
  652 continue
      if(bpowsmw .le. 0._R8) then
      do 654 j=1,npsi
  654 ajavlh(j) = 0._R8
      go to 700
               endif
 
 
!
!     set counter keeping track of LSC/TSC iterations
!
      ICLH = 0
!     re-cycle point if current is not well-converged
!
 1654 continue
 
      call lhwrite(nlhcd,0)
      call lsc(bpowsmw,nlhcd,nlhcd2,nlscgnu,nout,nlhcdin,nterm,irayt     &  
     &        ,iplot,ierr)
      iraytd = irayt
      if(ierr.gt.0) ineg=45
      if(ierr.lt.0) write(nout,1991) ierr
 1991 format("LSC returned ierr=",i5)
      if(ierr.lt.0) write(nterm,1991) ierr
      if(ineg.ne.0) return
      rewind nlhcd2
      read(nlhcd2,1001) (powtsc(l),l=2,npsit)
      read(nlhcd2,1001) (currtsc(l),l=2,npsit)
      read(nlhcd2,1001) (djdetsc(l),l=2,npsit)
      read(nlhcd2,1001,end=651) (djdets2(l),l=2,npsit)
 1001 format(1p5e16.6)
!
      do 655 ns=1,nspc
      read(nlhcd2,1001) (powtsci(ns,l),l=2,npsit)
  655 continue
!
  651 continue
!
      powsum = 0._R8
      cursum = 0._R8
      do 653 l=2,npsi
      powsum = powsum + powtsc(l)*vp(l)*dpsi
      cursum = cursum +                                                  &  
     & acoef(106)*currtsc(l)*1.E4_R8*vp(l)*dpsi/(tpi*xplas)
  653 continue
      write(91,9662) kcycle
      write(91,9661) (powtsc(l),l=2,npsit)
      write(91,9961)
      write(91,9663) kcycle
      write(91,9661) (currtsc(l),l=2,npsit)
      write(91,9961)
      write(91,9664) kcycle
      write(91,9661) (djdetsc(l),l=2,npsit)
      write(91,9961)
      write(91,9665) kcycle
      write(91,9661) (voltlp(l),l=2,npsit)
      write(91,9961)
      write(91,9666) kcycle
      write(91,9661) (.5*udsv*(as(l)+as(l-1)),l=2,npsit)
 9661 format(1p5e12.4)
 9662 format(" powtsc at cycle", i10)
 9663 format(" currtsc at cycle", i10)
 9664 format(" djdetsc at cycle", i10)
 9665 format(" voltlp at cycle", i10)
 9666 format(" udsv*as  at cycle", i10)
 9961 format(1x)
      iplot=0
      irayt=2
!
!
      do 660 j=2,npsit
      ajavlho(j) = ajavlh(j)
      djdelh(j) = 0._R8
      ajavlh(j) = acoef(106)*currtsc(j)*1.E4_R8*usdi*xplas/gzero
!     if( abs(voltlp(j)) .ge. .01)
!    1djdelh(j) = djdetsc(j)*ajavlh(j)*tpi / (voltlp(j)*usdv)
      djdelh(j) = djdets2(j)*1.E4_R8*usdi/(usdv)
  660 continue
      ajavlh(1) = ajavlh(2)
      djdelh(1) = djdelh(2)
      djdetsc(1) = djdetsc(2)
      if(ifk.eq.1) go to 700
!
!     correct current by Newton's iteration
!
      AJAVLHM = 0._R8
      do 710 j=2,npsit
      aterm = .5_R8*(etpara(j)+etpara(j-1))*xplas*djdelh(j)
      if(aterm.lt.0) aterm=0._R8
      fac = aterm / (1._R8+aterm)
      ajnew = ajavlh(j) + fac*(ajavlho(j) - ajavlh(j))
      ajavlh(j) = ajnew
      AJAVLHM = max(AJAVLHM,abs(ajavlh(j)))
      currtsc(j) = ajavlh(j) / (1.E4_R8*usdi*xplas/gzero)
  710 continue
      ajavlh(1) = ajavlh(2)
      currtsc(1) = currtsc(2)
!     write(nout,9663)
!     write(nout,9661) (currtsc(l),l=2,npsit)
!
      cursum2 = 0._R8
      do 753 l=2,npsi
      cursum2 = cursum2 + currtsc(l)*1.E4_R8*vp(l)*dpsi/(tpi*xplas)
  753 continue
      write(nout,9660) kcycle,iraytd,bpowsmw,powsum,cursum,cursum2
      write(nterm,9660) kcycle,iraytd,bpowsmw,powsum,cursum,cursum2
 9660 format(" LSC:N=",i7," iraytd=",i2," ,powin,out,cur",               &  
     &     1p4e 9.1)
!
!     see if the correction to the current is reasonably big,
!     and if so, rewrite the equilibrium file and recall LSC;
!     but not more than 50 times
      RATHLH = 0._R8
      do 1770 j=2,npsit
      DENOM = ABS(ajavlh(j))
      IF(DENOM .LE. .01_R8*AJAVLHM) GO TO 1770
      RATHLH = max(RATHLH,ABS(ajavlh(j)-ajavlho(j))/DENOM)
 1770 continue
      write(nterm,4700) kcycle,iclh,rathlh
 4700 format(i7,i3,1pe12.4)
      IF(RATHLH .LE. 0.1_R8) GO TO 750
      ICLH = ICLH+1
!
!     under-relax to avoid oscillations in current
!
      write(nout,1790) kcycle,iclh
      do 1780 j=2,npsit
      DENOM = ABS(ajavlh(j))
      IF(DENOM .LE. .01_R8*AJAVLHM) GO TO 1780
      diff = (ajavlh(j)-ajavlho(j))/DENOM
      if(abs(diff).le.0.1_R8) go to 1780
      write(nout,1781) j,ajavlh(j),ajavlho(j),voltlp(j),djdelh(j)        &  
     &     ,diff
 1780 continue
 1790 format(" kcycle=",i4,"  iclh=",i2,/,                               &  
     &"  j      ajavlh     ajavlho      voltlp      djdelh diff")
 1781 format(i3,1p6e12.4)
      backa = 0.5_R8
!     if(iclh.gt.10) backa = .25
!     if(iclh.gt.15) backa = .10
!     if(iclh.gt.25) backa = .05
      DO 1700, J=2, NPSIT
      AJAVLH(J) = AJAVLHO(J) + backa*( AJAVLH(J) - AJAVLHO(J) )
 1700 CONTINUE
!
!.....set switch to indicate e-field iteration
      irayt = 0
!
      IF(ICLH .LT. 50) GO TO 1654
      INEG = 45
  750 continue
!.....adjust power to be equal to desired input power
      if(powsum.ne.0) then
      fudgefac = bpowsmw/powsum
      write(91,9971) times,fudgefac,cursum,cursum2
 9971 format(1p4e12.4)
      do l=2,npsit
      powtsc(l) = powtsc(l)*bpowsmw/powsum
      enddo
      endif
  725 continue
      ajavlh(npsit+1) = ajavlh(npsit)
      djdelh(npsit+1) = djdelh(npsit)
      djdetsc(npsit+1)= djdetsc(npsit)
!
!
!.....define surface centered arrays
      do 720 j=1,npsit
        ajavlh2(j) = .5_R8*(ajavlh(j)+ajavlh(j+1))
        djdetsc2(j) = .5_R8*(djdetsc(j)+djdetsc(j+1))
  720 continue
  700 continue
!
      return
!
  701 continue
      do 702 j=1,npsit+1
      powtsc(j) = 0._R8
      currtsc(j) = 0._R8
      djdetsc(j) = 0._R8
      ajavlh(j) = 0._R8
      djdelh(j) = 0._R8
      ajavlh2(j) = 0._R8
      djdetsc2(j) = 0._R8
  702 continue
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
