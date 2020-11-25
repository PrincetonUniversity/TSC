!                     ./tsc_a/source/spdnstxo.f90          31 May 2010
!
      subroutine spdnstxo(ishot,stime,numdatapoints,                     &  
     &     si,ncurrent,                                                  &  
     &     iflux_active,rfl,zfl,sfl,numfluxloops,                        &  
     &     mirnov_active,bmirnov_r,bmirnov_z,bmirnov_tangle,             &  
     &     bmirnov_pangle, bmirnov_data,nmirnov,                         &  
     &     icurrent_active,scurrent,numcurrent ,                         &  
     &     halosensor_current, nhalosensors,                             &  
     &     expnedl,nnedl , expdifl,ndifl ,nbipower,nbi1p,maxte,tffactor,imachine)
!
!			          Rewrite to use new NSTX shot data files.
!			          ROS.  20 Nov 2009 
!
!			          Added output of TSC params from 18 Jan 2010 version.
!			          ROS.  28 May 2010 
      USE CLINAM
      USE POLCUR
      USE SAPROP
      USE SPDMOD
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER numdatapoints,nmirnov,nhalosensors,ncurrent,numfluxloops
      INTEGER numcurrent,nnedl,ndifl,nbi1p,ishot
      INTEGER icount,l,i,ios22,n,indx,iz,ip,lz,lp,k,ii, imachine
!
	integer nptpts
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 turnsa,tint,dtms,vesppe2,tffactor
      REAL*8 rline,difle,rdis,pval
      REAL*8 rlinae,abeamp,fac,pone,r1,z1,graddum,dpxarg,dpzarg
      REAL*8 gsdum,psarg,pxzarg,pxxarg,pzzarg,p1,anumfl
!============
!
      REAL*8 stime(numdatapoints),                                       &  
     & rfl(numfluxloops),zfl(numfluxloops),                              &  
     & bmirnov_r(nmirnov),bmirnov_z(nmirnov),                            &  
     & bmirnov_tangle(nmirnov),bmirnov_pangle(nmirnov),                  &  
     & expnedl(nnedl),expdifl(ndifl),nbipower(nbi1p),                    &  
     & maxte(numdatapoints)
 
      REAL*8 si(numdatapoints,ncurrent),                                 &  
     & bmirnov_data(numdatapoints,nmirnov),                              &  
     & sfl(numdatapoints,numfluxloops),                                  &  
     & scurrent(numdatapoints,numcurrent),				 &
     & halosensor_current (numdatapoints,nhalosensors)
 
      integer mirnov_active(nmirnov),iflux_active(numfluxloops),         &  
     &          icurrent_active(numcurrent)
      dimension turnsa(pngroup)
!
	real*8  qzer, qedge, preskpa, adii, ajt, fabs
	real*8  rmserrtime, avgrmserr, sumrmserr
!
!============
!From the coil currents data file:  NSTX naming conventions
! tinterp   Ip1 IOH  IPF1AU  IPF2U IPF3U IPF4 IPF5 IPF3L  IPF2L IPF1b  IPF1AL  ITF
!            1   2     3       4     5    6     7    8      9     10     11     12
!
	integer :: icnstxplas = 1, icnstxoh   = 2, icnstxpf1au = 3, icnstxpf2u = 4
	integer :: icnstxpf3u = 5, icnstxpf4  = 6, icnstxpf5   = 7, icnstxpf3l = 8
	integer :: icnstxpf2l = 9, icnstxpf1b =10, icnstxpf1al =11
      integer :: icnstxtf
!
	REAL*8 totcureo, ztimes, zapl, zpolcurchi
!
	integer :: debug = 0,  nparout = 12, numtermout, ntest = 11
!
        ncurrentg = ncurrent-2
        numfluxloopsg = numfluxloops
!
!......tf coil current
      icnstxtf = ncurrent
!
!.....check if initialization is necessary, otherwise proceed to 100
      if(ifrst(9).ne.1) go to 100
!------------------------------------------------------------------------------------
      ifrst(9) = 0
      icount = 48
      avgrmserr = 0.d0
      numtermout = 0
      write(nterm,8888) ishot
      write(nout, 8888) ishot
 8888 format ("spdnstxo: NSTX data files read for shot", i7)
      write(nout,8889) numdatapoints, numfluxloops ,numcurrent
 8889 format ("spdnstxo: numdatapoints, numfluxloops ,numcurrent = ",3i4)
!
!...........................................................
!
!...save trajectories from input file
!
!...........................................................
      allocate(pflux(numfluxloops))
      allocate(pfluxo(numfluxloops))
      allocate(eflux(numfluxloops))
      allocate(csuma(pngroup))
      allocate(ppcur(pngroup))
      allocate(fbcur(pngroup))
!...........................................................
      if(numfluxloops.gt.0)	then
				do i=1,numfluxloops
				pflux(i) = 0._R8
				enddo
				endif
!
!.....open the file to write the simulation and experimental data
      if( numargs .lt. 1 ) 	then
				filename = 'tscnstxdata_out'
      else
         filename = 'tscnstxdata_out' // '.' // trim(suffix)
				endif
      open (no167a,file=trim(filename),status='unknown',iostat=ios22)
!
!.....open the file to write the simulation parameters
      if( numargs .lt. 1 ) 	then
				filename = 'tscnstx_params_out'
      else
         filename = 'tscnstx_params_out' // '.' // trim(suffix)
				endif
      open (nparout,file=trim(filename),status='unknown',iostat=ios22)
      write (nparout,'("       Time Iplasma_kA      Betat       li/2       Xmag",&
     &"       Zmag         r0          a      kappa      delta       xsep",&
     &"       zsep         q0      qedge    betapol       liGA      te(1)      ne(1)    preskpa",&
     &"       ffac       tevv      thalo      whalo     fbchia")')
!
      write (ntest,'(/"spdnstxo: t         Ip    Sig_RMS  <Sig_RMS>")') 
!
       if(irst1.ne.1) go to 100
!					restart run, read to end of no167a
   98 continue
      read(no167a,7010, end=99, err=7012) ztimes,zapl,atot,csum,vesppe,vescure,totcure,    &
     &                   chicur,vescure2,zpolcurchi
      read(no167a,7010) (csuma(k),k=1,ncurrent-2)
      read(no167a,7010) (ppcur(k),k=1,ncurrent-2)
      read(no167a,7010) (fbcur(k),k=1,ncurrent-2)
      if(numfluxloops.gt.0) read(no167a,7010) (pflux(i),i=1,numfluxloops)
      if(numfluxloops.gt.0) read(no167a,7010) (eflux(i),i=1,numfluxloops)
      read(no167a,7011)
      go to 98
!
!.....end of initialization section
!------------------------------------------------------------------------------------
!
   99 continue
      backspace no167a
  100 continue
      do  n=2,numdatapoints
        indx = n
	if(stime(n) .ge. times) exit
	enddo
        iz     = indx - 1
        ip     = indx
        tint = stime(ip)-stime(iz)
        dtms = (times - stime(iz))/tint
!
        if(indx .ge. numdatapoints) dtms = 0._R8
!
!-----------------------------------------------------------------------------------------------
        if(imachine.eq.0) then    ! for NSTX
!
!.......OH
        gcur(istart,1) = (1._R8-dtms)*si(iz,icnstxoh)  + dtms*si(ip,icnstxoh)
!.......PF5
        gcur(istart,5) = (1._R8-dtms)*si(iz,icnstxpf5) + dtms*si(ip,icnstxpf5)
!
        if(isym.eq.1) 	then
!
!			.....isym=1
!.........PF3U and PF3L
          gcur(istart,3) = .5_R8*((1._R8-dtms)*si(iz,icnstxpf3u) + dtms*si(ip,icnstxpf3u))     &  
     &                   + .5_R8*((1._R8-dtms)*si(iz,icnstxpf3l) + dtms*si(ip,icnstxpf3l))
!.........PF1AU and PF1AL
          gcur(istart,8) = .5_R8*((1._R8-dtms)*si(iz,icnstxpf1au) + dtms*si(ip,icnstxpf1au))     &  
     &                   + .5_R8*((1._R8-dtms)*si(iz,icnstxpf1al) + dtms*si(ip,icnstxpf1al))
!.........PF2U AND PF2L
          gcur(istart,2) = .5_R8*((1._R8-dtms)*si(iz,icnstxpf2u) + dtms*si(ip,icnstxpf2u))     &  
     &                   + .5_R8*((1._R8-dtms)*si(iz,icnstxpf2l) + dtms*si(ip,icnstxpf2l))
!.........PF4
          gcur(istart,4) = (1._R8-dtms)*si(iz,icnstxpf4)  + dtms*si(ip,icnstxpf4)
!.........PF1B
          gcur(istart,6) = (1._R8-dtms)*si(iz,icnstxpf1b) + dtms*si(ip,icnstxpf1b)
!
        else
!	  write (nterm, '(" here is isym = 0:")')
!
!			.....isym=0
!.........PF3U
          gcur(istart,3) = ((1._R8-dtms)*si(iz,icnstxpf3u) + dtms*si(ip,icnstxpf3u))
!.........PF3L
          gcur(istart,10)= ((1._R8-dtms)*si(iz,icnstxpf3l) + dtms*si(ip,icnstxpf3l))
!.........PF1AU
          gcur(istart,8) = ((1._R8-dtms)*si(iz,icnstxpf1au)+ dtms*si(ip,icnstxpf1au))
!.........PF1AL
          gcur(istart,9) = ((1._R8-dtms)*si(iz,icnstxpf1al)+ dtms*si(ip,icnstxpf1al))
!.........PF2U
          gcur(istart,2) = ((1._R8-dtms)*si(iz,icnstxpf2u) + dtms*si(ip,icnstxpf2u))
!.........PF2L
          gcur(istart,6) = ((1._R8-dtms)*si(iz,icnstxpf2l) + dtms*si(ip,icnstxpf2l))
!.........PF4
!         gcur(istart,4) = (1.-dtms)*si(iz,icnstxpf4) + dtms*si(ip,icnstxpf4)
!.........PF1B
          gcur(istart,7) = (1._R8-dtms)*si(iz,icnstxpf1b) + dtms*si(ip,icnstxpf1b)
	endif
!
!.......Plasma Current
        vesppe  = (1._R8-dtms)*si(iz,icnstxplas) + dtms*si(ip,icnstxplas)
!
!--------------------------------------------------------------------------------------
      else   ! for DIII-D
!cj     si(iz,1) plasma current (ip)
!cj     si(iz,2:21) coil current (ip)
!cj     si(iz,22) tf current
!cj     do i=1,20
        do i=1,ncurrentg
          gcur(istart,i) = (1._R8-dtms)*si(iz,i+1) + dtms*si(ip,i+1)
        enddo
!
!.......Plasma Current
        vesppe  = (1._R8-dtms)*si(iz,1) + dtms*si(ip,1)
        
      endif
!--------------------------------------------------------------------------------------

      vesppe2  = 0.
      vescure = 0.d0
      totcure = 0.d0
      vescure2= totcure - vesppe2
!
!.....Chi injector (supply) current
      chicur = 0.d0
!     								
!
      rline  = (1._R8-dtms)*expnedl(iz) + dtms*expnedl(ip)
      difle =-((1._R8-dtms)*expdifl(iz) + dtms*expdifl(ip))
!										5 Nov 2009
      gzerov(istart) = abs(tffactor*((1._R8-dtms)*si(iz,icnstxtf) + dtms*si(ip,icnstxtf)))
!
!
      if(numfluxloops.gt.0) then
        do i=1,numfluxloops
          eflux(i) = (1._R8-dtms)*sfl(iz,i) + dtms*sfl(ip,i)
!.........special for DIII
!
!.........I'm not sure why this needs to be done....scj
!         if(imachine.eq.1 .and. i .gt. 1) eflux(i) = eflux(i) - eflux(1)
        enddo
      endif
!
!.....debug
      write(85,1085) iz,ip,times,eflux(1),sfl(iz,1),sfl(ip,1),eflux(7),sfl(iz,7),sfl(ip,7)
 1085 format(2i5, 3x,e10.2,3x,1p3e10.2,3x,1p3e10.2)
!
!.....put maximum temperature in global(217)
      global(217) = (1._R8-dtms)*maxte(iz) + dtms*maxte(ip)
!
      rdis = 0._R8
      do  i=iminn,imaxx
        if(iexv(i,nh).eq.1 .or. iexs(i,nh).eq.1) cycle
        pval = psi(i,nh)
        if(pval.gt.psilim) cycle
        rdis = rdis + deex
      enddo
      if(rdis.le.deex) rdis = deex
!				.....line averaged density
      rlinae = abs(rline)/rdis
      global(215) = rlinae*1.E-20_R8
      global(216) = difle
!
      do i=1,nbi1p
        if(nbipower(i).ne.0) go to 34
      enddo
      go to 37
   34 continue
      abeamp  = (1._R8-dtms)*nbipower(iz) + dtms*nbipower(ip)
!
!
      do l=1,ntpts
!cj     write(*,*) "--- spdnstxo L276 beamp(l) =", l, beamp(l) 
        beamp(l) = abeamp
      enddo
   37 continue
!
      if (kcycle < 1) write(nterm, '("spdnstxo: beamp =", 1p5e12.4)')  &
     &                 (beamp(l), l=1,5)
!
!
      do  l=1,ntpts
        fact(l) = 0._R8
      enddo
      fact(istart) = 1._R8
!
      pcur(istart) = vesppe

!cj      write (*,'(A,i4)') "spdnstxo: istart pcur gcur gzerov = ", istart
!cj      write (*,'(2e12.4)') pcur(istart), gzerov(istart)
!cj      write (*,'(1p6e12.4)') (gcur(istart,i), i=1,ncurrentg)
!cj      k=534
!cj      write (*,'(A,i4,2e12.4)') "spdnstxo: ntime = ", k,stime(k),times
!cj      do k = 2, numdatapoints
!cj         if( abs(stime(k) - times) .lt. 1.e-10) then
!cj      write (*,'(A,i4,2e12.4)') "spdnstxo: ntime = ", k,stime(k), times
!cj      write (*,'(1p6e12.4)') (si(k,i), i=2,ncurrent-1)
!cj         endif
!cj      enddo

      if(kcycle.lt.0) return
!
!.....flux from simulation
      pone = 0._R8
      if(numfluxloops.gt.0) then
        do 300 i=1,numfluxloops
 
          r1 = rfl(i)
          z1 = zfl(i)
          if(z1.lt.0 .and. isym.eq.1) z1 = -zfl(i)
          call grap(2,z1,r1,graddum,dpxarg,dpzarg,gsdum,psarg,pxzarg,pxxarg,pzzarg,1)
          p1  = -tpi*psarg
!
! 10/13/2010 scj  ... subtract off first loop flux for DIII-D
!         if(imachine.eq.1 .and. i.gt.1) p1 = p1 - pflux(1)
!
          pfluxo(i)  = pflux(i)
          if(kcycle.le.1) go to 299
          pflux(i)  = 0.5_R8*(p1  + pfluxo(i))
          go to 300
  299     pflux(i) = p1
  300   continue
!
      endif
!     if(imachine.eq.0) then
!.......set flux offset for NSTX
!       fbcon(1) = (.5_R8*(eflux(33)+eflux(34))-eflux(53))/tpi
!     endif
!
!.....calculate total current inside vessel and plasma
      csum = 0._R8
      do k=1,pngroup
        csuma(k)=0._R8
        turnsa(k) = 0._R8
      enddo
!
      do 400 ii=1,nwire
        fac = 1._R8
        if(isym.eq.1 .and. zcoil(ncoil-nwire+ii) .gt. 0) fac=2._R8
!
        if(igroupw(ii).le.ncurrent-2) go to 399
          csum = csum + ccoil(ncoil-nwire+ii)*udsi*fac
        go to 400
 399    continue
        do k=1,ncurrent-2
          indx = ncoil-nwire+ii
          if(igroupw(ii).eq.k) csuma(k)=csuma(k)+ccoil(indx)*udsi*fac
          if(igroupw(ii).eq.k) turnsa(k)=turnsa(k)+aturnsc(indx)*fac
        enddo
!
 400  continue
      atot = apl + csum
      do k=1,ncurrent-2
        ppcur(k) = gcur(istart,k)*turnsa(k)
        fbcur(k) = gcurfb(k)*udsi*turnsa(k)
      enddo
!
      anumfl = numfluxloops
!
!             apl  = plasma current      - simulation 
!             atot = vv + plasma current - simulation 
!             csum = vv current          - simulation 
!             vesppe  = plasma current      - experiment 
!             vescure = vv current          - experiment 
!             totcure = vv + plasma current - experiment 
!
!.....only write to o167a every 50 time steps
      icount = icount + 1
      if(icount .lt. 50) return
      icount = 0
      write(no167a,7010) times,apl,atot,csum,vesppe,vescure,totcure,    &  
     &                   chicur,vescure2,polcurchi
      write(no167a,7010) (csuma(k),k=1,ncurrent-2)
      write(no167a,7010) (ppcur(k),k=1,ncurrent-2)
      write(no167a,7010) (fbcur(k),k=1,ncurrent-2)
      if(numfluxloops.gt.0) write(no167a,7010) (pflux(i),i=1,numfluxloops)
      if(numfluxloops.gt.0) write(no167a,7010) (eflux(i),i=1,numfluxloops)
      write(no167a,7011)
!..........				Write rms error
		numtermout = numtermout + 1
		rmserrtime = 0.d0
		do i=1,numfluxloops
		rmserrtime = rmserrtime + (pflux(i)-eflux(i))**2
		enddo
		rmserrtime = dsqrt (rmserrtime/dfloat(numfluxloops))
		sumrmserr = sumrmserr + rmserrtime
		avgrmserr = sumrmserr/dfloat(numtermout)
      write (ntest,'(1p8e11.3e1)')   times, 1.d-3*apl, rmserrtime, avgrmserr
		if (mod(numtermout,10) .eq. 0)	then
      write (nterm,'("spdnstxo: t, Ip, Sig_RMS, <Sig_RMS> =",1p8e11.3e1)')   &
     &            times, 1.d-3*apl, rmserrtime, avgrmserr
     						endif
!
!..........				Write out some TSC params to nparout
	qzer = qprof(1)
	qedge = qprof(npts)
	if (isurf.ne.0)	then
			qzer = qprof2p
			qedge = qprof2(npts)
			preskpa = adp(1)/vpg(1)*udsp*1.E-3_R8
			endif
!
      write (nparout,'(1p30e11.4e1)')  times, 1.d-3*apl ,betat0,ali2, xmag,zmag, rmajor, rminor, &
     &    shape5, shape3, shape6,shape7, qzer,qedge, betapol,aliga, te(1), 1.d-20*ane(1), preskpa, &
     &    ffac, tevv, thalov(istart), whalov(istart), fbchia(istart)
!
      return
 7012 write(nterm,'(" *** error in spdnstxo ***, ishot = ",2i6)') ishot
      ineg=37
      return
!
 7010 format(1p10e12.4)
 7011 format(1x)
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
