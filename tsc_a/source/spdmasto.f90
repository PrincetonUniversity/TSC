      subroutine spdmasto(ishot,stime,numdatapoints,                     &  
     &     si,ncurrent,                                                  &  
     &     iflux_active,rfl,zfl,sfl,numfluxloops,                        &  
     &     mirnov_active,bmirnov_r,bmirnov_z,bmirnov_tangle,             &  
     &     bmirnov_pangle, bmirnov_data,nmirnov,                         &  
     &     icurrent_active,scurrent,numcurrent ,                         &  
     &     expnedl,nnedl )
!
!
!
      USE CLINAM
      USE POLCUR
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER numdatapoints,ncurrent,numfluxloops,nmirnov
      INTEGER numcurrent,nnedl,ishot,nptpts,nflparam,icount,l,i
      INTEGER ios22,n,indx,iz,ip,lz,lp,k,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pprsav,rnormsv,fbchisv,pcursv,vlpsav,xmagzsv,zmagzsv
      REAL*8 tevv0sv,ffac0sv,zeffvsv,alphrsv,betarsv,frparsv
      REAL*8 thalosv,whalosv,heactsv,gzerosv,pflux,pfluxo,eflux
      REAL*8 csuma,turnsa,ppcur,fbcur,tint,dtms,vesppe,rline,rdis
      REAL*8 pval,rlinae,fac,pone,r1,z1,graddum,dpxarg,dpzarg,gsdum
      REAL*8 psarg,pxzarg,pxxarg,pzzarg,p1,csum,atot,anumfl
      REAL*8 vescure2,chicur,totcure,vescure
!============
      parameter(nptpts=60,nflparam=100)
!
!
!
      dimension pprsav(nptpts),rnormsv(nptpts), fbchisv(nptpts),         &  
     & pcursv(nptpts),vlpsav(nptpts),                                    &  
     & xmagzsv(nptpts),zmagzsv(nptpts),tevv0sv(nptpts),                  &  
     & ffac0sv(nptpts),zeffvsv(nptpts),                                  &  
     & alphrsv(nptpts),betarsv(nptpts),frparsv(nptpts),                  &  
     & thalosv(nptpts),whalosv(nptpts),heactsv(nptpts),gzerosv(nptpts)
!
 
!
      dimension pflux(nflparam),pfluxo(nflparam),eflux(nflparam)
!
      REAL*8 stime(numdatapoints),                                       &  
     & rfl(numfluxloops),zfl(numfluxloops),                              &  
     & bmirnov_r(nmirnov),bmirnov_z(nmirnov),                            &  
     & bmirnov_tangle(nmirnov),bmirnov_pangle(nmirnov),                  &  
     & expnedl(nnedl)
 
 
      REAL*8 si(numdatapoints,ncurrent),                                 &  
     & bmirnov_data(numdatapoints,nmirnov),                              &  
     & sfl(numdatapoints,numfluxloops),                                  &  
     & scurrent(numdatapoints,numcurrent)
 
      integer mirnov_active(nmirnov),iflux_active(numfluxloops),         &  
     &          icurrent_active(numcurrent)
      dimension csuma(pngroup),turnsa(pngroup),ppcur(pngroup),           &  
     &    fbcur(pngroup)
!
!
!
!.....check if initialization is necessary, otherwise proceed to 100
      if(ifrst(9).ne.1) go to 100
      ifrst(9) = 0
      icount = 48
      write(nterm,8888) ishot
      write(nout, 8888) ishot
 8888 format("netcdf data file read for shot",i7)
      write(nout,8889) numdatapoints, numfluxloops ,numcurrent
 8889 format(" numdatapoints, numfluxloops = ",3i4)
 
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
      gzerosv(l) = gzerov(l)
    5 continue
      do i=1,numfluxloops
      pflux(i) = 0._R8
      enddo
!
!.....open the file to write the simulation and experimental data
      if( numargs .lt. 1 ) then
        filename = 'tscmastdata_out'
      else
        filename = 'tscmastdata_out' // '.' // trim(suffix)
      end if
      open(no167a,file=trim(filename),status='unknown',iostat=ios22)
!
!
!.....end of initialization section
!
!
!
  100 continue
!
      do 888 n=2,numdatapoints
      indx = n
      if(stime(n) .ge. times) go to 110
!
  888 continue
 110  continue
!     write(nterm,9002) indx
!9002 format(" indx = ",i3)
      iz     = indx - 1
      ip     = indx
      tint = stime(ip)-stime(iz)
      dtms = (times - stime(iz))/tint
!
      if(indx .ge. numdatapoints) dtms = 0._R8
!.....OH
      gcur(istart,1) = (1._R8-dtms)*si(iz,15) + dtms*si(ip,15)
!
      if(isym.eq.0) go to 201
!
!.....isym=1
!.....P3A and P3B
      gcur(istart,3) = .5_R8*((1._R8-dtms)*si(iz,1) + dtms*si(ip,1))     &  
     &               + .5_R8*((1._R8-dtms)*si(iz,2) + dtms*si(ip,2))
!.....P4A and P4B
      gcur(istart,4) = .5_R8*((1._R8-dtms)*si(iz,3) + dtms*si(ip,3))     &  
     &               + .5_R8*((1._R8-dtms)*si(iz,4) + dtms*si(ip,4))
!.....P5A AND P5B
      gcur(istart,5) = .5_R8*((1._R8-dtms)*si(iz,5) + dtms*si(ip,5))     &  
     &               + .5_R8*((1._R8-dtms)*si(iz,6) + dtms*si(ip,6))
!.....P2IA and P2IB
      gcur(istart,2) = .5_R8*((1._R8-dtms)*si(iz,7) + dtms*si(ip,7))     &  
     &               + .5_R8*((1._R8-dtms)*si(iz,9) + dtms*si(ip,9))
!.....P2OA and P2OB
      gcur(istart,7) = .5_R8*((1._R8-dtms)*si(iz,8) + dtms*si(ip,8))     &  
     &               + .5_R8*((1._R8-dtms)*si(iz,10) + dtms*si(ip,10))
!.....P6UA and P6UB
      gcur(istart,6) = .5_R8*((1._R8-dtms)*si(iz,11) + dtms*si(ip,11))   &  
     &               + .5_R8*((1._R8-dtms)*si(iz,13) + dtms*si(ip,13))
!.....P6LA and P6LB
      gcur(istart,8) = .5_R8*((1._R8-dtms)*si(iz,12) + dtms*si(ip,12))   &  
     &               + .5_R8*((1._R8-dtms)*si(iz,14) + dtms*si(ip,14))
!
      go to 202
  201 continue
!     write(nterm,9004)
!9004 format(" after 201 ")
!
!.....isym=0
!.....P3A and P3B
      gcur(istart,3) = ((1._R8-dtms)*si(iz,1) + dtms*si(ip,1))
      gcur(istart,9) = ((1._R8-dtms)*si(iz,2) + dtms*si(ip,2))
!.....P4A and P4B
      gcur(istart,4) = ((1._R8-dtms)*si(iz,3) + dtms*si(ip,3))
      gcur(istart,10)= ((1._R8-dtms)*si(iz,4) + dtms*si(ip,4))
!.....P5A AND P5B
      gcur(istart,5) = ((1._R8-dtms)*si(iz,5) + dtms*si(ip,5))
      gcur(istart,11)= ((1._R8-dtms)*si(iz,6) + dtms*si(ip,6))
!.....P2IA and P2IB
      gcur(istart,2) = ((1._R8-dtms)*si(iz,7) + dtms*si(ip,7))
      gcur(istart,12)= ((1._R8-dtms)*si(iz,9) + dtms*si(ip,9))
!.....P2OA and P2OB
      gcur(istart,7) = ((1._R8-dtms)*si(iz,8) + dtms*si(ip,8))
      gcur(istart,13)= ((1._R8-dtms)*si(iz,10) + dtms*si(ip,10))
!.....P6UA and P6UB
      gcur(istart,6) = ((1._R8-dtms)*si(iz,11) + dtms*si(ip,11))
      gcur(istart,14)= ((1._R8-dtms)*si(iz,13) + dtms*si(ip,13))
!.....P6LA and P6LB
      gcur(istart,8) = ((1._R8-dtms)*si(iz,12) + dtms*si(ip,12))
      gcur(istart,15)= ((1._R8-dtms)*si(iz,14) + dtms*si(ip,14))
 
!
  202 continue
      vesppe = ((1._R8-dtms)*si(iz,16) + dtms*si(ip,16))
!
      rline  = (1._R8-dtms)*expnedl(iz) + dtms*expnedl(ip)
!.....note:  the following line was for NSTX
!     gzerov(istart)=abs(0.0000072*((1.-dtms)*si(iz,11)+dtms*si(ip,11)))
!
      do i=1,numfluxloops
      eflux(i) = (1._R8-dtms)*sfl(iz,i) + dtms*sfl(ip,i)
      enddo
!
 
      rdis = 0._R8
      do 38 i=iminn,imaxx
      if(iexv(i,nh).eq.1 .or. iexs(i,nh).eq.1) go to 38
      pval = psi(i,nh)
      if(pval.gt.psilim) go to 38
      rdis = rdis + deex
   38 continue
      if(rdis.le.deex) rdis = deex
!
!.....line averaged density
      rlinae = abs(rline)/rdis
      global(215) = rlinae*1.E-20_R8
!
!     write(nterm,9999) rdis,rline,rlinae
!9999 format(" rdis,rline,rlinae ", 1p3e12.4)
!
!
      do 35 l=1,ntpts
      lz = l
      lp = l + 1
      if((tpros(lz).le.times).and.(tpros(lp).gt.times)) go to 36
 35   continue
 36   continue
!
      fac = (times-tpros(lz))/(tpros(lp)-tpros(lz))
      ppres(istart) = pprsav(lz) + fac*(pprsav(lp)-pprsav(lz))
      vloopv(istart) = vlpsav(lz) + fac*(vlpsav(lp)-vlpsav(lz))
      rnorm(istart) = rnormsv(lz) + fac*(rnormsv(lp)-rnormsv(lz))
      fbchia(istart) = fbchisv(lz) + fac*(fbchisv(lp)-fbchisv(lz))
      xmagz(istart) = xmagzsv(lz) + fac*(xmagzsv(lp) - xmagzsv(lz))
      zmagz(istart) = zmagzsv(lz) + fac*(zmagzsv(lp) - zmagzsv(lz))
      tevv0(istart) = tevv0sv(lz) + fac*(tevv0sv(lp) - tevv0sv(lz))
      ffac0(istart) = ffac0sv(lz) + fac*(ffac0sv(lp) - ffac0sv(lz))
      zeffv(istart) = zeffvsv(lz) + fac*(zeffvsv(lp) - zeffvsv(lz))
      frcparv(istart) = frparsv(lz) + fac*(frparsv(lp) - frparsv(lz))
      betarv(istart) = betarsv(lz) + fac*(betarsv(lp) - betarsv(lz))
      alpharv(istart) = alphrsv(lz) + fac*(alphrsv(lp) - alphrsv(lz))
      thalov(istart) = thalosv(lz) + fac*(thalosv(lp) - thalosv(lz))
      whalov(istart) = whalosv(lz) + fac*(whalosv(lp) - whalosv(lz))
      heactv(istart) = heactsv(lz) + fac*(heactsv(lp) - heactsv(lz))
      gzerov(istart) = gzerosv(lz) + fac*(gzerosv(lp) - gzerosv(lz))
!
      do 40 l=1,ntpts
   40 fact(l) = 0._R8
      fact(istart) = 1._R8
!
      pcur(istart) = vesppe
      if(kcycle.lt.0) return
!
!.....flux from simulation
      pone = 0._R8
      do 300 i=1,numfluxloops
 
      r1 = rfl(i)
      z1 = zfl(i)
      if(z1.lt.0 .and. isym.eq.1) z1 = -zfl(i)
      call grap(2,z1,r1,graddum,dpxarg,dpzarg,gsdum,psarg,               &  
     &          pxzarg,pxxarg,pzzarg,1)
      p1  = -tpi*psarg
!
      pfluxo(i)  = pflux(i)
      if(kcycle.le.1) go to 299
      pflux(i)  = 0.5_R8*(p1  + pfluxo(i))
      go to 300
  299 pflux(i) = p1
  300 continue
!      write(nterm,1300) rfl(33),zfl(33),pflux(33),eflux(33),
!     1                  rfl(34),zfl(34),pflux(34),eflux(34),
!     2                  rfl(53),zfl(53),pflux(53),eflux(53)
! 1300 format(" flux loop 31 r and z",1p4e12.4,/,
!     1       " flux loop 36 r and z",1p4e12.4,/,
!     2       " flux loop 53 r and z",1p4e12.4)
      fbcon(1) = (.5_R8*(eflux(33)+eflux(34))-eflux(53))/tpi
!
      icount = icount + 1
      if(icount .lt. 50) return
      icount = 0
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
      if(igroupw(ii).le.10) go to 399
      csum = csum + ccoil(ncoil-nwire+ii)*udsi*fac
      go to 400
 399  continue
      do k=1,10
      indx = ncoil-nwire+ii
      if(igroupw(ii).eq.k) csuma(k)=csuma(k)+ccoil(indx)*udsi*fac
      if(igroupw(ii).eq.k) turnsa(k)=turnsa(k)+aturnsc(indx)*fac
      enddo
!
 400  continue
      atot = apl + csum
      do k=1,10
      ppcur(k) = gcur(istart,k)*turnsa(k)
      fbcur(k) = gcurfb(k)*udsi*turnsa(k)
      enddo
!
      anumfl = numfluxloops
      vescure2 = 0._R8
      chicur = 0._R8
      totcure = 0._R8
      vescure = 0._R8
      write(no167a,7010) times,apl,atot,csum,vesppe,vescure,totcure,     &  
     &                   chicur,vescure2,polcurchi
      write(no167a,7010) (csuma(k),k=1,10)
      write(no167a,7010) (ppcur(k),k=1,10)
      write(no167a,7010) (fbcur(k),k=1,10)
      write(no167a,7010) (pflux(i),i=1,numfluxloops)
      write(no167a,7010) (eflux(i),i=1,numfluxloops)
!
      write(no167a,7011)
!
      return
 7010 format(1p10e12.4)
 7011 format(1x)
 9500 continue
 9501 continue
      write(nterm,1009) ishot
 1009 format(" *** error in spdmast ***, ishot = ",2i6)
      ineg=37
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
