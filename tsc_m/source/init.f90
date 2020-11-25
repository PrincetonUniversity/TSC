      subroutine init
!
!.....initialize everything for time loop
!
      USE CLINAM
      USE SAPROP
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!...set values for ifrst array (9 elements currently used, 7/26/85 )
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,ios2,l,nskip2o,nrecordo,nthvar,kcount,n,j,ig,ii,ngr
      INTEGER iabs, index, kskip
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 aplfb,fac,anum,adenom,rojmin,psimid,prmin,ppmin,qmin
      REAL*8 gsval,gpval,gppval,rval,rpval
!============
      do 5 i = 1,50
      ifrst(i) = 1
    5 continue
      ifrst(2) = 0
      ifrst(5) = 0
      ifrst(6) = 0
      ifrst(7) = 0
!
!.....set default values for input
      call defolt
!
!.....read in all input
      call inpt
      if(ineg.ne.0) return
      if(ifk.ne.1) go to 50
!
!.....read data file tscouta for ifk .ne. 0    (a is suffix)
!     iouttsc(1:6) = 'tscout'
!     iouttsc(7:7) = isuffix(1:1)
      if( numargs .lt. 1 ) then
         iouttsc = 'tscout' // isuffix(1:1)
      else
         iouttsc = 'tscout' // '.' // trim(suffix)
      end if
      open(nlhcd2,file=trim(iouttsc),status='old',iostat=ios2)
      read(nlhcd2,1001) (powtsc(l),l=2,npsit)
      read(nlhcd2,1001) (currtsc(l),l=2,npsit)
      close(nlhcd2)
 1001 format(1p5e16.6)
   50 continue
!
!
!.....for a start run, jump around
      if(irst1 .ne. 1) go to 10
      if(kcycle.lt.0) go to 141
      call prset
      nskip2o =nskip2
      nrecordo = nrecord
      nthvar = pnsave
      nskip2 = (ncycle)/nthvar + 1
!
!.....rearrange time storage for restart run
      lenscr = pgls+1+3*pncoil+10*pngroup+6
      kcount = 0
      nrecord = 0
      if(acoef(5) .gt. 0) go to 102
      kskip = 1
      if(kcycle.gt.0) kskip = ncycle*nrecordo/(kcycle*pnsave)
      if(kskip.lt.1) kskip = 1
      write(*,1002) nrecordo,nskip2,nskip2o,ncycle,nthvar,kskip
 1002 format(" restart run: nrecordo,nskip2,nskip2o,ncycle,nthvar,kskip =",6i8)
      do 100 n=1,nrecordo,kskip
        nrecord = nrecord + 1
        apld0sv(nrecord) = apld0sv(n)
        do l=1,lenscr
          index = (nrecord-1)*lenscr+l
          if(index .ge. 2*pnsave*pglobp) then
            ineg=63
            write(*,*) "pltsav index exceeds dimension", index
            return
          endif
          pltsav(index)=pltsav((n-1)*lenscr+l)
        enddo
        do nn=1,nobs-1,2
          fluxu(nn,nrecord) = fluxu(nn,n)
          fluxl(nn,nrecord) = fluxl(nn,n)
        enddo
        if(acoef(3001) .eq. 1.0_R8) then
          do nn=1,30
            do j=1,pnx
              ufdata(nrecord,j,nn) = ufdata(n,j,nn)
            enddo
          enddo
        endif
  100 continue

      go to 103
  102 continue
      nskip2 = (ncycle-kcycle)/nthvar + 1
      do 101 i=1,pgls
      globmax(i) = -1.E30_R8
      globmin(i) =  1.E30_R8
  101 continue
      gvmin = 0._R8
      gvmax = 0._R8
      gvmin0 = 0._R8
      gvmax0 = 0._R8
      gpmin = 0._R8
      gpmax = 0._R8
      gemin = 0._R8
      gemax = 0._R8
      gpmint = 0._R8
      gpmaxt = 0._R8
      gemint = 0._R8
      gemaxt = 0._R8
  103 continue
      iplt2 = nskip2
      nplrec = 0
      if(icirc.eq.1) call cirequi
!
!.....modify things at restart time
      call rsmod
      return
!
   10 continue
!
!
!.....define constants and initialize counters
      call define
!
!.....zero out arrays
!
      call zero
!
!
!.....generate grid
!
      call setup
!
!.....initialize circuit equations
      if(icirc.eq.1) call cirequi
      if(idata.eq.7) call tcv1
      if(ncycle .lt. 0) return
      if(ineg.ne.0) return
!
      aplfb = 0._R8
      do 805 ig=1,pngroup
  805 gcurfb(ig) = 0._R8
!
      if(idata.eq.0) go to 140
      if(idata.eq.1) call spdpbx
      if(idata.eq.2) call spdtftr
      if(idata.eq.3) call spdd3d
      if(idata.eq.4) call spdpbxm
      if(idata.eq.5) call spdpbxn
      if(idata.eq.6) call fedtsc
      if(idata.eq.8) call spdasdex
      if(idata.eq.9 .or. idata.eq.11) call spdnstx
      if(idata.eq.10) call spdmast
      p0 = ppres(istart)*usdp
      tcuro = pcur(istart)*usdi
      fac = sqrt(nx*(nz-1)/((2-isym)*4._R8*(alx-ccon)*alz))
      anum = (vloopv(istart)-acoef(850))*ndiv
      adenom = etav*fac*acoef(11)*udsv*tpi
      if(adenom .ne. 0._R8) tcuro = tcuro + anum/adenom
!
      r0 = rnorm(istart)
      gzero = gzerov(istart)
!
!......define coil current arrays for idata.ne.0
      do 138 ii=1,nwire
      ngr = iabs(igroupw(ii))
      n = ncoil-nwire+ii
      cwire0(ii)=aturnsw(ii)*gcur(istart,ngr)*usdi
      ccoil(n) = cwire0(ii)
      ccoils(n) = cwire0(ii)*udsi
  138 continue
!
      if(ncoil.eq.nwire) go to 140
      do 139 ii=1,ncoil-nwire
      ngr = iabs(igroupc(ii))
      ccoil(ii) = aturnsc(ii)*gcur(istart,ngr)*usdi
  139 ccoils(ii) = ccoil(ii)*udsi
  140 continue
!
!.....check if reading initial equilibrium from disk
!
      if(irst1.ne.2) go to 141
      call eqread
      xmag = xplas
      zmag = zplas
      imag = (xplas-ccon)/deex + 2
      jmag = (zplas-zzero)/deez + 2.0000001_R8
      if(ineg.ne.0) return
      go to 142
  141 continue
!
!...............................................................
!.....solve for initial equilibrium
!...............................................................
      call initeq
!
      psib = 0._R8
      psiob = 0._R8
      if(ineg .ne. 0) return
!.....define q and g arrays
!
      rojmin = fracn0*r0
      do 18 i=2,nxp+1
      if(ajey(i).lt.0) ajey(i) = abs(ajey(i))
      if(ajey(i).eq.0) ajey(i) = 1._R8
      do 18 j=2,nzp+1
      rdens(i,j) = 1._R8
!
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      prmin = tevv*rojmin*udsd/(.5_R8*udsh)
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 17
      if(psimid.gt.psilim) go to 17
      call peval(psimid,1,prmin,ppmin,i,j)
   17 continue
!
      qmin = prmin**(3._R8/5._R8)*ajey(i)
      q(i,j) = qmin
      call geval(psimid,1,gsval,gpval,gppval,i,j)
      if(xsqoj(i).ne.0) g(i,j) = gsval/xsqoj(i)
!
   19 call reval(psimid,idens,0,rval,rpval,i,j)
      roj(i,j) = rval
      r(i,j) = rval*ajey(i)
      ro(i,j) = r(i,j)
   18 continue
  142 continue
!
      do 20 i=2,nxp+1
      if(ajey(i).lt.0) ajey(i) = abs(ajey(i))
      if(ajey(i).eq.0) ajey(i) = 1._R8
      do 20 j=2,nzp+1
      bmagy(i,j) = 1
      qo(i,j) = q(i,j)
      pr(i,j) = (q(i,j)/ajey(i))**(5._R8/3._R8)
      pro(i,j) = pr(i,j)
      go(i,j) = g(i,j)
   20 psio(i,j) = psi(i,j)
!
!......define gs array
      call gsdef
!
!
      if(lrswtch.ne.0) go to 21
      call limpsi
      call magaxis
      call delpdef
      if(psilim.le.psimin) ineg=34
      if(ineg.ne.0) return
      if(isvd.ge.1) call xptcalc
!
      call ajdef
      call icalc
      npts = nx/2
!
      if(npsi.gt.1 .and. isurf.gt.0) then
      npts = npsi-1
      call flxvol(1)
        endif
!
      if(ineg.ne.0) return
      if(isurf.eq.0) go to 21
      call flxvol(2)
      if(ineg.ne.0) return
   21 continue
      call prdef
      call auxdef
      call boundc
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
