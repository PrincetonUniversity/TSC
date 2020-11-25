!#include "f77_dcomplx.h"
      subroutine globe
!
!.....calculate global integrals for output
!
      USE CLINAM
      USE POLCUR
      USE RUNAWAY
      USE SAPROP
      USE SCR1
      USE SPECIE
      USE SVDCOM
      USE SAWCOM
      USE SCR11


      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ifrstg,ipelav,i,j,icount,iplas,jplas,ipass,ig,ifail,n
      INTEGER iabs,jz,jm,npt95,nmax,indx,npsihm,l,ii,igr,izero,iw
      INTEGER jw,iii,i1,i2,i3
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 gapdw,dgapdw,xshape,fshape
!     REAL*8 grsum0,gvsum,gvsum0,gcur0ka,gcurfka,sumvsg,x1old
!     REAL*8 grsum0,gvsum,gvsum0,gcur0ka,gcurfka,sumvsg
      REAL*8 dtglobe,rojmin,pvac,sum5,sum6,sum7,sum8,psii,psij,ai
      REAL*8 aj,xe,ajac,vx,vz,gis,gjs,oi,oj,bx,bz,vt,bt,ajt,etat
      REAL*8 gsi,gsj,ajx,ajz,psimid,pval,ppval,eint,rdis,rbtave
      REAL*8 sum2,sum3,sum4,fac,alength,btor,rlin,rval,rpval,diamag
      REAL*8 tflux,sump,sumt,ratipol,oldipol,tauvvpol,oldtim,deldia
      REAL*8 olddia,delipol,flvvpol,f1,tauems2,tauemn,ellip,delta
      REAL*8 x1,x2,z1,z2,asp,psi1,psi2,psi3,psi4,psisn1,ans,anull
      REAL*8 adipol,aquad,ahex,aoct,adec,anull0,adipol0,aquad0
      REAL*8 ahex0,aoct0,adec0,zden,zr,za,zip,zk,zd,zp,zb,zq
      REAL*8 ztohmic,ztaux,tauegl,taueka,zl,tauerl,qmur,rmur,q95
      REAL*8 qcylin,fkd,qstar,denom,ancritd,ancrit,dt1s,dt2s,dt3s
      REAL*8 sum1,tcpu,rsaw,pohmp,pohmh,apld,xcons,zcons,gradsq
      REAL*8 dpsidx,dpsidz,gsval,psval,dum1,dum2,dum3,fac1
      REAL*8 dtpergl,sumw,facw,term,sumdn
      REAL*8 sum, AREAL
      real*8 qprofmin
      integer lqmin
!============
!     common /sawcom/  dwcore, dwcorec, dw, dwc, dwc2, ratio,scrit,
!    1dwbussac,dwelong,dwko,dwfast
      dimension gapdw(6), dgapdw(6)
      dimension xshape(4),fshape(4)
!    &         ,grsum(pngroup),grsum0(pngroup),
!    1 gvsum(pngroup),gvsum0(pngroup),gcur0ka(pngroup),gcurfka(pngroup),
!    2 sumvsg(pngroup)
!     data  x1old / 0._R8/
      REAL*8 :: x1old = 0._R8
      data ifrstg / 1 /
      data ipelav / 1 /
      save ifrstg, ipelav
!============      
      INTEGER :: istat = 0 
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: grsum
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: grsum0
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: gvsum
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: gvsum0
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: gcur0ka
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: gcurfka
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sumvsg
!============      
      INTERFACE

         SUBROUTINE GROUPCUR(grsum,grsum0,gvsum,gvsum0,gcur0ka,gcurfka)
         USE CLINAM
         USE SCRATCH
         USE SCR11, ONLY : gcurr
         IMPLICIT NONE
         INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
         INTEGER iig,ig,l,jg,iw1,iabs,ic1,ic
         REAL*8 sum0,sumv,sumt,sumve,fac,sumvi
         REAL*8 sum
         REAL*8, DIMENSION(*) :: grsum,grsum0,gvsum,gvsum0
         REAL*8, DIMENSION(*) :: gcur0ka, gcurfka

         END SUBROUTINE GROUPCUR

      END INTERFACE

!============

!     IF(.not.ALLOCATED(grsum)) ALLOCATE( grsum(pngroup), STAT=istat)
!     IF(.not.ALLOCATED(grsum0)) ALLOCATE( grsum0(pngroup), STAT=istat)
!     IF(.not.ALLOCATED(gvsum)) ALLOCATE( gvsum(pngroup), STAT=istat)
!     IF(.not.ALLOCATED(gvsum0)) ALLOCATE( gvsum0(pngroup), STAT=istat)
!     IF(.not.ALLOCATED(gcur0ka)) ALLOCATE( gcur0ka(pngroup), 
!    &                                 STAT=istat)
!     IF(.not.ALLOCATED(gcurfka)) ALLOCATE( gcurfka(pngroup), 
!    &                                 STAT=istat)
      IF(.not.ALLOCATED(sumvsg)) ALLOCATE( sumvsg(pngroup), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : globe  ' 
!============      
!
!.....calculate  volume integrals
!
!
!.....compute energys
!
!
      sum1o = sum1e
      sum2o = sum2e
      sum3o = sum3e
      sum4o = sum4e
      sum7o = sum7e
!     ener0 = enermj
!     enermj = uint*1.e-6
      dtglobe = times - tolds
!
      sum1e = 0._R8
      sum2e = 0._R8
      sum3e = 0._R8
      sum4e = 0._R8
      sum7e = 0._R8
      rojmin = fracn0*r0
      pvac = (tevv*udsd*rojmin/(.5_R8*udsh))
      sum5 = 0._R8
      sum6 = 0._R8
      sum7 = 0._R8
      sum8 = 0._R8
!
!.....calculate volume integrals of kinetic energy,
!.....internal energy, and ohmic dissipation
!
!
      do 200 i=3,nxp
      do 150 j=3,nzp
      psii = .5_R8*(psi(i,j)+psi(i,j-1)-psi(i-1,j)-psi(i-1,j-1))
      psij = .5_R8*(psi(i,j)+psi(i-1,j)-psi(i,j-1)-psi(i-1,j-1))
      ai = .5_R8*(abig(i,j)+abig(i,j-1)-abig(i-1,j)-abig(i-1,j-1))
      aj = .5_R8*(abig(i,j)+abig(i-1,j)-abig(i,j-1)-abig(i-1,j-1))
      xe = xarh(i)
      ajac = xe*dxdz
      vx = (aj*deex )/ajac
      vz = ( ai*deez)/ajac
!
      gis = (deez**2)*(xe/ajac)**2
      gjs = (deex**2)*(xe/ajac)**2
      oi = .5_R8*(omeg(i+1,j)-omeg(i-1,j))
      oj = .5_R8*(omeg(i,j+1)-omeg(i,j-1))
      vz = vz  +deez*oj*gjs
      vx = vx + deex*oi*gis
  100 continue
      bx = (psij*deex )/ajac
      bz = ( psii*deez)/ajac
      vt = w(i,j)/(ajac*xe)
      bt = g(i,j)*xe/ajac
      ajt = .25_R8*(ajphi(i,j)+ajphi(i,j-1)+ajphi(i-1,j)+ajphi(i-1,j-1))    
      etat = .25_R8*(etay(i,j)+etay(i,j-1)+etay(i-1,j)+etay(i-1,j-1))
      gsi = .5_R8*(g(i+1,j)*xsqoj(i+1)-g(i-1,j)*xsqoj(i-1))
      gsj = .5_R8*(g(i,j+1)*xsqoj(i)-g(i,j-1)*xsqoj(i))
      ajx = (gsj*deex )/ajac
      ajz = ( gsi*deez)/ajac
!
      sum1e = sum1e + (1+isym)*.5_R8*(vx**2+vz**2)*ajac*tpi
      sum2e = sum2e + (1+isym)*.5_R8*vt**2*ajac*tpi
      sum3e = sum3e + (1+isym)*.5_R8*(bx**2+bz**2)*ajac*tpi
      sum5=sum5+(1+isym)*etay(i,j)*(ajx**2+(1+isym)*ajz**2)*ajac*tpi
      sum6 = sum6 + (1+isym)*etat*ajt**2*ajac*tpi
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      pval = pvac
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 4
      sum4e = sum4e + (1+isym)*.5_R8*(bt**2-(gzero/xe)**2)*ajac*tpi
      if(psimid.gt.psilim) go to 4
      call peval(psimid,1+isurf,pval,ppval,i,j)
    4 continue
      sum7e = sum7e + (1+isym)*1.5_R8*ajac*tpi*pval
  150 continue
  200 continue
      ekino = ekin
      ekin = (sum1e + sum2e)
      eint = (sum3e + sum4e)
      ekinp = sum1e
      ekint = sum2e
      eintp = sum3e
      eintt = sum4e
!
      if(dtglobe.le.0._R8) go to 501
      ener3 = (sum3e-sum3o)/dtglobe
      ener4 = (sum4e-sum4o)/dtglobe
      ener7 = (sum7e-sum7o)/dtglobe
!     if(ener0 .eq. 0) go to 501
!     wdotmw = (enermj-ener0)/dtglobe
  501 continue
      ener5 = sum5
      ener6 = sum6
      ener9 = 0
!
!
      if(amach.gt.acoef(21) .or. ekin.gt.acoef(22)) then
      write(nout,4777) kcycle,times,dts,dtmins,amach,acoef(21),          &  
     &              ekin,acoef(22)
 4777 format(" kcycle,times,dts,dtmins,amach,acoef(21),ekin,acoef(22) ",  &  
     &                                                                   &  
     &      /,i3,1p7e12.4)
      ineg=4
      endif
      rdis = deex
      rbtave = gzero
      if(rbtave .le. 0._R8 ) rbtave = 1._R8
!
!
      sum = 0._R8
      sum2 = 0._R8
      sum3 = 0._R8
      sum4 = 0._R8
      sum5 = 0._R8
      sum6 = 0._R8
      sum7 = 0._R8
      icount = 0
      if(lrswtch.eq.0) then
      do 99 j=jminn,jmaxx
      fac=1.0_R8
      if(isym.eq.1 .and. j.ne.2) fac = 2.0_R8
      do 99 i=iminn,imaxx
      alength = tpi*xary(i)
      if(psi(i,j).ge.psilim) go to 99
      sum = sum + dxdz/(etay(i,j)*alength)*fac
      icount = icount + 1
   99 continue
!
      sum2 = 0._R8
      fac = 1.0_R8
      if(isym.eq.1 ) fac=2.0_R8
      do 90 i=iminn,imaxx
      do 90 j=jminn+1,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 89
      xe = xarh(i)
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      if(psimid.ge.psilim)go to 89
      sum5 = sum5 + fac*g(i,j)
      sum2 = sum2 + fac*(g(i,j)-gzero/xsqoj(i))
      btor = g(i,j)*xsqoj(i)/xe
!****      sum3 = sum3 + fac*.5*(btor**2)*ajey(i)
        sum3 = sum3 + fac * ajey(i)
      sum4 = sum4 + fac*pr(i,j)*ajey(i)
      sum6 = sum6 + fac*btor*xe
      sum7 = sum7 + fac
   89 continue
   90 continue
      rbtave = gzero
      if(sum7.ne.0) rbtave = sum6/sum7
      if(rbtave.eq.0) rbtave = 1._R8
!
!.....line average density
      rlin = 0._R8
      rdis = 0._R8
      do 38 i=iminn,imaxx
      if(iexv(i,nh).eq.1 .or. iexs(i,nh).eq.1) go to 38
      pval = psi(i,nh)
      if(pval.gt.psilim) go to 38
      call reval(pval,idens,isurf,rval,rpval,i,nh)
      rlin = rlin + rval*deex*udsd
      rdis = rdis + deex
   38 continue
      if(rdis.le.deex) rdis = deex
!
!
      enmjs = 1.5_R8*sum4*tpi*udsp
      diamag = sum2
      tflux = sum5
      if(icount.lt. .02_R8*nx*nz) sum = 0._R8
      if(sum3.ne.0.and.gzero.ge.0)  beta = 2._R8*sum4*(xplas/gzero)**2/  &  
     & sum3
!
      endif
      sump = 0._R8
      sumt = 0._R8
!
      global(1) = amach
      global(2) = ekin
      global(3) = ekinp
      global(4) = ekint
      global(5) = iplim
      global(6) = zmag
      global(7) = xmag
      global(8) = apl
      global(9) = tcurdtp*tpi*udsi
      global(10) = (psilim-psimin)*tpi
      global(11) = (psimin)*tpi
      global(12) = diamag
        if (isym.eq.1 .and. acoef(61).eq.0)   then
                global(6) = gzero
                                              endif
!
!.....added 4/24/05...scj
        if(acoef(79) .eq. 1.) global(7) = deltmks
        if (isym.eq.1 .and. acoef(61).eq.0 .and. lrswtch.gt.0) then
                global(5) = 5.E6_R8* (gzero - gs(iplas,jplas))
                if (kcycle.le.0) global(5) = 0
                ipass = ipass + 1
                if (mod(ipass,5).eq.0 .and. lrswtch.gt.0)  then
                if (global(5).gt.1.E-8_R8)   ratipol = oldipol /         &  
     & global(5)
                if (ratipol  .gt.1.001_R8)   tauvvpol= (times - oldtim)  &  
     &                                              /log(ratipol)
                deldia  = diamag - olddia
                delipol = global(5) - oldipol
                if (abs(delipol).gt.1.E-8_R8)   flvvpol = -deldia/       &  
     & delipol
!
                oldipol = global(5)
                olddia  = diamag
                oldtim  = times
                write (nout,4000)   times, gzero, gs(iplas,jplas),       &  
     &             global(5), tauvvpol, diamag, deldia, delipol,         &  
     & flvvpol
 4000         format(' t, g0, gpl, I, TAU, dia, deld, deli, L =',9e10.2)    
                                                            endif
                                                                endif
      if(irfp.eq.1 .or. acoef(504).ne.0) global(12) = helic
!     global(13) = reboun*udsv*tpi
      global(13) = 0._R8
      if(tolds2.eq.times .or. tolds2.le.0) go to 339
      global(13) = tpi*(psilim-pold2)/(times-tolds2)
      if(pold.eq.0) go to 339
      vsec = vsec + tpi*(psilim-pold)
      vsec2 = vsec2 + reboun*udsv*tpi*(dtglobe)
  339 pold2 = pold
      tolds2 = tolds
      pold = psilim
      tolds = times
      if(global(13).gt.acoef(14)) global(13) = acoef(14)
      if(global(13).lt.-acoef(14)) global(13) = -acoef(14)
!
      global(14) = qprof(1)
      if(isurf.ne.0) global(14) = 2._R8*qprof2(2) - qprof2(3)
      global(15) = abs(qprof(npts))
      if(irfp.eq.1) global(15) = qprof(npts)
      if(isurf.ne.1) go to 340
      if(psis1 .eq. psis2) go to 340
      f1 = (psis3-psis1)/(psis1-psis2)
      if(f1.gt.1._R8) f1 = 1._R8
      if(f1.lt.0._R8) f1 = 0._R8
      global(15) = abs((1._R8+f1)*qprof2(npsit) - f1*qprof2(npsit-1))
  340 continue
!
      global(16) = dts
      global(17) = betat
      if(isurf.ne.1) global(17) = beta
      if(global(18).ge.10._R8) global(18) = 10._R8
      global(20) = betapol+ali2
      global(21) = betapol
      global(22) = ali2
      global(23) = ali2
      do 6666 ig=20,22
      if(global(ig).lt.0) global(ig)=0._R8
 6666 if(global(ig).gt.10._R8) global(ig)=10._R8
!.....define global(24) at end to be kaye-goldstone confinement time
!
      if(kcycle.le.0 .and. tauems.gt.1.E6_R8) tauems = 0.5_R8*1000._R8*  &  
     & acoef(46)
      global(25) = tauems
      tauems2 = tauems
      if(hfluxav.ne.0._R8 .and. kcycle.gt.nskipsf) tauems2 = enerst2/    &  
     & hfluxav                                                           &  
     &                                                 *1000._R8
      tauemn = 2._R8*tauems
      if(kcycle.lt.10*nskipsf.and.tauems2.gt.tauemn) tauems2=tauemn
      global(26) = tauems2
      if(global(25).lt.0) global(25)=0._R8
      if(global(26).lt.0) global(26)=0._R8
      if(global(25).gt.1000._R8*acoef(46)) global(25) = 1000._R8*        &  
     & acoef(46)
      if(global(26).gt.1000._R8*acoef(46)) global(26) = 1000._R8*        &  
     & acoef(46)
      global(27) = te(2)
      global(28) = ti(2)
      global(29) = teform
      global(30) = tiform
      global(31) = chiauxs
      global(32) = chiohms
      global(33) = hfluxav*1.E-6_R8
!
!
!.....define initial guesses for shape parameters
      xshape(1) = xmag
      xshape(2) = sqrt(vary(npsit)/(tpi*pi*xmag))
      xshape(3) = shape3
      xshape(4) = shape4
      if(isurf.le.0.or.lrswtch.ne.0 .or. kmax .le. 0) go to 137
!
      call shapenp( 4,xshape,fshape,1.E-8_R8,ifail)
!     if(ifail .ne. 0)       ineg=15
      if(ineg.ne.0) return
      rmajor = xshape(1)
      rminor = xshape(2)
      shape3 = xshape(3)
      shape4 = xshape(4)
      if(rminor .ne. 0) shape5 = ellmom/rminor
      shape6 = xsep(1)
      shape7 = zsep(1)
      global(34) = rmajor
      global(35) = rminor
      global(36) = shape3
      global(37) = shape5
      global(38) = shape6
      global(39) = shape7
      go to 139
  137 continue
      call shape(ellip,delta,x1,x2,z1,z2)
      asp = (alx+ccon)/(alx-ccon)
      if(x2.ne.x1) asp = abs((x1+x2)/(x2-x1))
      global(34) = .5_R8*(x1+x2)
      global(35) = abs(.5_R8*(x2-x1))
      global(36) = delta
      global(37) = ellip
      global(38) = xsep(1)
      global(39) = zsep(1)
      if(irfp.ne.1) go to 139
      call fluxmod(psi1,psi2,psi3,psi4,psisn1,2.40_R8,0.40_R8)
      global(36) = psi1
      global(37) = psi2
      global(38) = psi3
      global(39) = psi4
      rmajor = global(34)
      rminor = global(35)
      shape3 = global(36)
      shape5 = global(37)
      shape6 = global(38)
      shape7 = global(39)
!
!
  139 continue
      global(40) = vsec
      if(kcycle.le.0) pzero = psimin
      global(41) = tpi*(psimin-pzero)
!.....upper half plane
!
      sum = 0
      do ig=1,pngroup
      sumvsg(ig) = 0
      enddo
!
      do 78 n=1,nwire
      ans = 0._R8
      if(xplas.eq.xary(iwire(n)).and.zplas.eq.zary(jwire(n)))            &  
     &    go to 78
      call gf(ineg,nmult,xplas,zplas,xary(iwire(n)),zary(jwire(n)),ans)
      ig = iabs(igroupw(n))
      sumvsg(ig) = sumvsg(ig) + ans*ccoil(ncoil-nwire+n)
   78 sum = sum + ans*ccoil(ncoil-nwire+n)/tpi
      if(ncoil.eq.nwire) go to 80
      do 79 n=1,ncoil-nwire
      call gf(ineg,nmult,xplas,zplas,xcoil(n),zcoil(n),ans)
      ig = iabs(igroupc(n))
      sumvsg(ig) = sumvsg(ig) + ans*ccoil(n)
   79 sum = sum + ans*ccoil(n)/tpi
   80 continue
      if(isym.eq.0) go to 180
!
!.....lower half plane
!
      do 178 n=1,nwire
      call gf(ineg,nmult,xplas,zplas,xary(iwire(n)),-zary(jwire(n)),ans)    
      if(zary(jwire(n)).eq.0) ans = 0
      ig = iabs(igroupw(n))
      sumvsg(ig) = sumvsg(ig) + ans*ccoil(ncoil-nwire+n)
  178 sum = sum + ans*ccoil(ncoil-nwire+n)/tpi
      if(ncoil.eq.nwire) go to 180
      do 179 n=1,ncoil-nwire
      call gf(ineg,nmult,xplas,zplas,xcoil(n),-zcoil(n),ans)
      if(zcoil(n).eq.0) ans = 0._R8
      ig = iabs(igroupc(n))
      sumvsg(ig) = sumvsg(ig) + ans*ccoil(n)
  179 sum = sum + ans*ccoil(n)/tpi
  180 continue
!
!...mod 05/13/01...group v-sec (sumvsg) stored in pgls+3*pncoil+5*pngroup+i
      global(42) = sum*tpi
      global(43) = reboun*udsv*tpi
      global(44) = vsec2
      global(45) = 1.E-6_R8*(pohmic + paux + palpha)
      global(46) = 1.E-6_R8*palphap
      global(47) = 1.E-6_R8*pohmic
      global(48) = 1.E-6_R8*paux
!
!.....multipolar moments of external field
      call multip2(anull,adipol,aquad,ahex,aoct,adec,                    &  
     &          anull0,adipol0,aquad0,ahex0,aoct0,adec0)
      global(49) = .001_R8*udsi*anull
      global(50) = .001_R8*udsi*anull0
      global(51) = .001_R8*udsi*adipol
      global(52) = .001_R8*udsi*adipol0
      global(53) = .001_R8*udsi*aquad
      global(54) = .001_R8*udsi*aquad0
      global(55) = .001_R8*udsi*ahex
      global(56) = .001_R8*udsi*ahex0
      global(57) = .001_R8*udsi*aoct
      global(58) = .001_R8*udsi*aoct0
      global(59) = .001_R8*udsi*adec
      global(60) = .001_R8*udsi*adec0
!
!....calculate kaye-goldston scaling...define global(24) here
!
!....first define density(10**20), major radius(m), minor radius(m),
!    plasma current(Ma),ellipticity,power(Mw),toroidal field,safety factor
      if(irfp.eq.1) go to 191
      if(isurf.le.0.or.lrswtch.ne.0.or.kmax.le.0.or.rdis.le.0) go to     &  
     & 237
      zden = (rlin/rdis)*1.E-20_R8
      zr = rmajor
      za = rminor
      zip = tcurdtp*tpi*udsi*1.E-6_R8
      if( zr  .le.0 .or. zip .le.0) go to 237
      zk = shape5
!     zk = vol/(2.*pi*rmajor)/(pi*rminor**2)
      zd = shape3
      zp = (pohmic+paux+prad+palpha)*1.E-6_R8
      zb = gzero/zr
      zq = 2.5_R8*za**2*zb*(1._R8+zk**2)/(zip*zr)
!
!....note that in these formula, density is in 10**20, not in
!         10**19 as Kaye sometime uses.  Also, effective mass is 2.5
!
!
      ztohmic = 0.07_R8*zden*zr**2*za*zq
      if(zp.gt.0 .and. za.gt.0 .and. zb.gt.0)                            &  
     &ztaux = 0.071_R8*zip**1.24_R8*zr**1.65_R8*zk**0.28_R8*zden**       &  
     & 0.26_R8                                                           &  
     &      /(zp**0.58_R8*za**0.49_R8*zb**0.09_R8)
      if(ztohmic.gt. 0 .and. ztaux .gt. 0)                               &  
     &tauekg = sqrt(1._R8/((1._R8/ztohmic)**2 + (1._R8/ztaux)**2))
!
      global(24) = tauekg*1000._R8
      if(global(24).gt.1000._R8*acoef(46)) global(24) = 1000._R8*        &  
     & acoef(46)
      if(global(24).lt.0.0_R8) global(24) = 0.0_R8
!
!.....Goldston  L-Mode
      if(za.gt.0 .and. zp.gt.0)                                          &  
     &ztaux = .048_R8*zk**0.5_R8*zip*zr**1.75_R8/(za**0.37_R8*zp**       &  
     & 0.50_R8)
      if(ztohmic.gt. 0 .and. ztaux .gt. 0)                               &  
     &tauegl = sqrt(1._R8/((1._R8/ztohmic)**2 + (1._R8/ztaux)**2))
!
      global(131) = tauegl*1000._R8
      if(global(131).gt.1000._R8*acoef(46)) global(131) = 1000._R8*      &  
     & acoef(46)
      if(global(131).lt.0.0_R8) global(131) = 0.0_R8
!
!.....ITER98(y,2) scaling(1998)  (corrected 3/14/02..scj)
      taueka = ztohmic
      if(zp.gt.0)                                                        &  
     & taueka = .172_R8*zk**.78_R8*zip**0.93_R8*zden**.41_R8*zb**.15_R8  &  
     &       *za**.58_R8*zr**1.39_R8/zp**.69_R8
!      in NSTX
!    & taueka = .1442*zk**.78*zip**0.93*zden**.41*zb**.15*(2.0)**0.19    &
!    &       *za**.58*zr**1.39/zp**.69
!
      global(132) = taueka*1000._R8
      if(global(132).gt.1000._R8*acoef(46)) global(132) = 1000._R8*      &  
     & acoef(46)
      if(global(132).lt.0.0_R8) global(132) = 0.0_R8
!
!.....Rebut-Lallia
      zl = (zr*za**2*zk)**0.333_R8
      if(zp.gt.0 .and. zeff.gt.0)                                        &  
     & ztaux = 0.324_R8*zip**0.5_R8*zden**0.75_R8*zb**0.5_R8*zeff**      &  
     & 0.25_R8                                                           &  
     &       *zl**(2.75_R8)/zp + .027_R8*zip*zl**1.5_R8/zeff**0.5_R8
      if(ztohmic.gt. 0 .and. ztaux .gt. 0)                               &  
     &tauerl = sqrt(1._R8/((1._R8/ztohmic)**2 + (1._R8/ztaux)**2))
!
      global(133) = tauerl*1000._R8
      if(global(133).gt.1000._R8*acoef(46)) global(133) = 1000._R8*      &  
     & acoef(46)
      if(global(133).lt.0.0_R8) global(133) = 0.0_R8
!
      go to 193
!
 191  continue
!
!.....for RFP only
      global(23) = tflux
      global(24) = gzero/rbtave
      global(25) = asp*tcurdtp/rbtave
      go to 237
 193  continue
!
!.....compute murikami density limit
      qmur = 5._R8*za**2*zk*zb/(zr*zip)
      if(qmur.le.1) qmur = 1._R8
      rmur = 2.0E20_R8*gzero/(xplas**2*qmur)
      global(18) = rlin/(rdis*rmur)
      global(19) = (rlin/rdis)*zr/zb
      global(149) = (rlin/rdis)*pi*za**2*1.E-20_R8/zip
!
!
!.....calculate q95
      qprofmin = 1.0e30
!     do 701 j=2,npsit
      do 701 j=1,npsit
      jz = j
      if(xsv2(j).gt.psi95) go to 702
      if(qprof2(j) .lt. qprofmin) then
      qprofmin = qprof2(j)
      lqmin = j
      endif
  701 continue
  702 continue
      jm = jz-1
      q95 = qprof2(jm)+(psi95-xsv2(jm))/(xsv2(jz)-xsv2(jm))              &  
     &                *(qprof2(jz)-qprof2(jm))
      qcylin = zq
      fkd = 1.24_R8-0.54_R8*zk+0.3_R8*(zk**2+zd**2)+.13_R8*zd
      qstar = zq*(1._R8+(za/zr)**2*(1._R8+.5_R8*(betapol+ali2)**2))*fkd
      qstar = zq
      denom = (zeff-1._R8)*zr*za**2*zk*(zq-2._R8/zk)
      if(denom.gt.0)                                                     &  
     & ancritd=2.5E19_R8*sqrt(1.0E-6_R8*abs(pohmic+paux+palpha)          &  
     &                                   *zq/denom)
      ancrit = 0._R8
      if(ancritd.ne.0 .and. kcycle.gt.0) ancrit = (rlin/rdis)/ancritd
      global(61) = abs(q95)
      global(62) = abs(qstar)
      global(63) = abs(qcylin)
      global(64) = del95
      global(65) = el95
      global(66) = xsepcal
      global(67) = zsepcal
      global(68) = psisep
      global(69) = psepcal
      global(70) = ancrit
      if(global(70) .gt. 10._R8) global(70) = 10._R8
      if(qcylin.ne.0) global(74) = 1._R8/qcylin
!     global(74) = qprofmin
      global(75) = betat0
!
      call peval(psimin,1+isurf,pval,ppval,imag,jmag)
      denom = gzero**2/xplas**2
      global(76) = pval/denom
      global(81) = el90
      global(82) = del90
      global(83) = vol
      global(84) = rtot
  237 continue
      if(isurf.eq.0) then
        npt95 = .95_R8*npts
        if(npt95 .lt.1) npt95 = 1
        global(61) = abs(qprof(npt95))
                     endif
!
      global(77) = rlin/rdis*1.E-20_R8
        if(idata.eq.2) global(77)=rlinv*1.E-20_R8
!
      global(78) = rvolav*1.E-20_R8
      if(lrswtch.eq.0)                                                   &  
     &call reval(psimin,idens,isurf,rval,rpval,imag,jmag)
      global(79) = rval*udsd*1.E-20_R8
      global(80) = uint*1.E-6_R8
      if(irfp.eq.1) global(80) = AREAL(nloop)
      dt1s = dt1*udst
      dt2s = dt2*udst
      dt3s = dt3*udst
      global(71) = max(dtmins,dt1s)
      global(72) = max(dtmins,dt2s)
      global(73) = max(dtmins,dt3s)
!
      global(85) = eps1a
      global(86) = eps1c
      global(87) = eps10
      global(88) = eps2a
      global(89) = eps2c
      global(90) = eps20
      global(91) = eps3a
      global(92) = eps3c
      global(93) = eps30
      global(94) = eps4a
      global(95) = eps4c
      global(96) = eps40
!
      if(kcycle.lt.50) go to 779
      global(97) = ener7*udsp*1.E-6_R8
      global(98) = ener3*udsp*1.E-6_R8
      global(99) = ener4*udsp*1.E-6_R8
      global(100) = tave1*1.E-6_R8
      global(101) = tave2*1.E-6_R8
      global(102) = tave3*1.E-6_R8
      global(103) = tave4*1.E-6_R8
      global(104) = tave5*1.E-6_R8
      global(105) = 0._R8
  779 continue
      do 778 n=97,104
  778 global(105) = global(105)+global(n)
      global(106) = tevv
      global(107) = ffac
      global(108) = npsit
      global(109) = resid
!
!.......compute global energys from powers
      do 781 n=1,9
      ebsav(n) = ebsav(n) + dtglobe*global(96+n)
  781 global(109+n) = ebsav(n)
!
!.....calculate energy sources(sinks) in Watts
      if(isurf.eq.0 .or. lrswtch.ne.0) go to 783
      sum1 = 0._R8
      sum2 = 0._R8
      sum3 = 0._R8
      sum4 = 0._R8
      sum5 = 0._R8
      sum6 = 0._R8
      sum7 = 0._R8
      do 782 j=2,npsit
      sum1 = sum1 + vp(j)*savebre(j)
      sum2 = sum2 + vp(j)*savecyc(j)
      sum3 = sum3 + vp(j)*(saveimp(j) + sradion(j)*usdp/usdt)
      sum4 = sum4 + vp(j)*(savefw(j)+savifw(j))
      sum5 = sum5 + vp(j)*(savee(j)+savei(j))
      sum6 = sum6 + vp(j)*(savelh(j)+savilh(j))
      sum7 = sum7 + vp(j)*(savebm(j)+savibm(j))
  782 continue
      global(119) = -sum1*dpsi*udsp/udst*1.E-6_R8
      global(120) = -sum2*dpsi*udsp/udst*1.E-6_R8
      global(121) = -sum3*dpsi*udsp/udst*1.E-6_R8
      global(122) = sum4*dpsi*udsp/udst*1.E-6_R8
!..
!..(possibly) temporary change....output Ploss in global(122)
!..instead of ech
!..
!     global(122)=global(33)-global(120)-global(119)
      global(241) = sum7*dpsi*udsp/udst*1.E-6_R8
      global(123) = sum5*dpsi*udsp/udst*1.E-6_R8
      global(124) = sum6*dpsi*udsp/udst*1.E-6_R8
  783 continue
!
      call second(tcpu)
      x1 = tcpu/60._R8
      global(125) = tlastg + (x1-x1old)
      if(dtglobe.le.0 .or. kcycle .le. 1) go to 784
      global(126) = (global(125)-tlastg)/dtglobe
  784 continue
      tlastg = global(125)
      x1old = x1
!
      global(129) = alphades
      global(130) = alphabes
      global(134) = alphafrs
      global(135) = zeff
      global(136) = rvolavi*1.E-20_R8
      global(137) = rvolavh*1.E-20_R8
      global(138) = wdotmw
      if(global(132) .gt. 0._R8) global(139)=global(25)/global(132)
      if(global(139) .gt. 4._R8) global(139)=4._R8
      if(global(132) .gt. 0._R8) global(140)=global(26)/global(132)
      if(global(140) .gt. 4._R8) global(140)=4._R8
!
      global(141) = bsitot
      global(142) = cditot
      global(198) = fwitot
      global(242) = ecitot
      global(143) = ain99
      global(144) = ain95
      global(145) = ain90
      global(146) = alhitot
      global(147) = xsep(2)
      global(148) = zsep(2)
!
!.....note: global(131)-(133) defined above
!           global(149) defined above
      iplas = (xplas-ccon)/deex + 2
      jplas = (zplas-zzero)/deez + 2
      global(150) = gs(iplas,jplas)
      call sawrad(rsaw)
      if(rsaw .le. 0.)then
      global(151) = rminora(lqmin)
      else 
      global(151) = rsaw
      endif
      global(152) = aliga
!
!  ....special for ZTH
!     if(irfp .ne. 1) go to 1153
!     global(153) = voltlp1
!     global(154) = voltlp2
!     global(155) = voltlp3
!     global(156) = voltlp4
!     global(157) = voltlp5
!     global(158) = vtran
!     global(159) = radius1
!     global(160) = radius2
!     global(161) = radius3
!     global(162) = radius4
!     global(163) = radius5
      global(164) = thalos
      global(165) = whalos
      pohmp = 0._R8
      pohmh = 0._R8
!
!.....calculate the ohmic power dissipated in halo and plasma regions
      if(whalos.ge.0) call ohmcalc(pohmp,pohmh)
      global(166) = pohmh
      global(167) = pohmp
!1153 continue
      if(isurf.le.0 .or. lrswtch.ne.0 .or. kmax.le.0) go to 789
!
!DEFINE:  168 - internal volt-sec  -- (Poynting)
!         169 - resistive "
!         170 - Ejima Coeff.
!
      apld = tcurdtp*tpi
      global(168) = apld*xmag*ali2
      if(kcycle.le.0) vzero = global(168)
      global(169) = vsec - (global(168)-vzero)
      if(global(169).lt.0) global(169) = 0._R8
      global(170) = global(169) / (apld*rmajor)
  789 continue
      global(171) = qprof2(2)
      global(172) = qprof2(3)
!.....put li and Ct in global(173) and (174)
      global(173) = 2._R8*global(23)
      if(global(9) .ne. 0) global(174) = global(75) * global(35)         &  
     &                     *gzero*1.E8_R8/ (global(34)*global(9))
      global(175) = pohmic + paux + palpha
      global(176) = palpha
      global(177) = acoef(113)
!
!.....global(178) - global(197) used to store shape control info
      if(acoef(296) .eq. 6._R8) then
      ncnt = 6
      call gaps(gapdw(1),gapdw(2),gapdw(3),gapdw(4),gapdw(5),            &  
     & gapdw(6),dgapdw(1),dgapdw(2),dgapdw(3),dgapdw(4),dgapdw(5),       &  
     & dgapdw(6))
      endif
      if(acoef(296) .eq. 7._R8) then
      ncnt = 6
      call gaps2(gapdw(1),gapdw(2),gapdw(3),gapdw(4),gapdw(5),           &  
     & gapdw(6),dgapdw(1),dgapdw(2),dgapdw(3),dgapdw(4),dgapdw(5),       &  
     & dgapdw(6))
      endif
      if(ncnt .le. 0) go to 1609
      nmax = min(10,ncnt)
      do 1608 i=1,nmax
      indx = 178 + 2*(i-1)
      if(acoef(296).eq. 6._R8.or. acoef(296).eq. 7._R8) go to 1602
      xcons = 0
      zcons = 0
      do 1601 j=1,ntpts
      xcons = xcons + fact(j)*xcon0(j,i)
      zcons = zcons + fact(j)*zcon0(j,i)
 1601 continue
      call grap(2,zcons,xcons,gradsq,dpsidx,dpsidz,gsval,psval,          &  
     &          dum1,dum2,dum3,1)
      global(indx)   = 0
      global(indx+1) = 0
      if(psisep .ge. 1.E20_R8.or. psisep.eq.0) go to 1608
      global(indx)   =  psval - psisep
      global(indx+1) = (psval - psisep) / max(sqrt(gradsq),1.E-8_R8)
      go to 1608
 1602 continue
      global(indx) = gapdw(i)
      global(indx+1)=dgapdw(i)
 1608 continue
 1609 continue
!....impurity pellet radius
      global(200) = dndtpel
      global(201) = xpel
      global(202) = tepel
      global(203) = totradp
      global(204) = radpel
      global(205) = sumpel
      global(206) = global(78)*global(83)*1.E20_R8
      global(207) = iqtrubmax
      global(208) = nskipsf
      global(209) = recur
      global(210) = apl
      global(211) = sreav
      global(212) = sumre
      global(213) = global(141)+global(142)+global(146)+global(198)      &  
     &            + ajpress + global(242)
      call extfield(global(214))
!.....global(215) is defined in subroutine spnstxo for idata=9 or idata=11
!.....global(216) is defined in subroutine spnstxo for idata=9 or idata=11
!
!.....global(217) defined in subroutine spnstxo for idata=9 or idata=11
      npsihm = int(pwidthc * AREAL(npsit))
      global(218) = te(npsihm)
      global(219) = ti(npsihm)
      global(220) = polcurchi
      global(221) = polcurchi2
      global(222) =-dwcore
      global(223) = dwcorec
      global(224) =-dw
      global(225) = dwc
      global(226) = dwc2
      global(227) = ratio
      global(228) = scrit
      global(229) = dwbussac
      global(230) = dwelong
      global(231) = dwko
      global(232) = dwfast
      global(233) = ajpress
      global(234) = (psisep-psimin)*tpi
      global(235) = tcurdtpcl*tpi*udsi
      global(236) = 0.
      if(whalos.gt.0.) global(236) = (phalo-psimin)*tpi
!
!......calculate area in plasma and within halo region
      call acalc(psilim,global(237))
      global(238) = 0.
      if(whalos.gt.0) call acalc(phalo,global(238))
      global(239) = global(238) - global(237)
!     write(nterm,8866) kcycle, ajpress
!8866 format(" kcycle=",i4,"ajpress=",1pe12.4)
!
!
!.....global(243) - global(250) not presently used
!
!..............................................................
!
!
!
      do 601 n=1,ppsi
      global(pglobs+n) = AREAL(ibaloon(n))
  601 continue
!
!.....divertor strike points
      if(iplate.ne.0 .and. nplate.ne.0 .and. lrswtch.eq.0)call divplat
      indx = pglobs + ppsi
      do 602 n=1,pnplat
      indx = indx + 1
      global(indx) = strike(n,1)
      indx = indx + 1
      global(indx) = strike(n,2)
      indx = indx + 1
      global(indx) = dsep(n,1)
      indx = indx + 1
      global(indx) = dsep(n,2)
      do 603 l=1,nseg(n)
      indx = indx + 1
  603 global(indx) = hplate(n,l)
  602 continue
      if(ifrstg.eq.0 .or. irst1.eq.1) go to 238
      do 338 i=1,pgls
      globmax(i) = -1.E30_R8
      globmin(i) = 1.E30_R8
  338 continue
  238 continue
      ifrstg=0
!
      do 138 i=1,pgls
      if(global(i).eq.0._R8) go to 138
      globmax(i) = max(global(i),globmax(i))
      globmin(i) = min(global(i),globmin(i))
  138 continue
!
      do 30 i=pgls+1,pgls+ncoil
      ii = i-pgls
      global(i) = ccoil(ii)*udsi*.001_R8
      if (kcycle.lt.2) go to 30
      if(global(i).lt.curmin) curmin = global(i)
      if(global(i).gt.curmax) curmax = global(i)
      igr = iabs(igroupc(ii))
      gigpmax(igr) = max(global(i),gigpmax(igr))
      gigpmin(igr) = min(global(i),gigpmin(igr))
   30 continue
!
      izero = pgls+ncoil
      do 40 i=izero+1,izero+ncoil
      n = i - izero
      ii = n - (ncoil-nwire)
      if(ii.ge.1) go to 41
!
!.....voltage for external coil   --> modify here for icirc=1 <--
      if(icirc.eq.0) global(i) = rscoils(n)*ccoil(n)*udsi*.001_R8
      if(icirc.eq.1) global(i) = vegplot(n)*udsv*.001_R8
      go to 45
   41 continue
!.....voltage for internal coil
      iw = iwire(ii)
      jw = jwire(ii)
      global(i) = resave(ii)*tpi*udsv*.001_R8
   45 continue
!
!.....power
      global(i+ncoil) = global(i)*ccoil(n)*udsi*.001_R8
      igr = iabs(igroupc(n))
      gvgpmax(igr) = max(global(i),gvgpmax(igr))
      gvgpmin(igr) = min(global(i),gvgpmin(igr))
      gpgpmax(igr) = max(global(i+ncoil),gpgpmax(igr))
      gpgpmin(igr) = min(global(i+ncoil),gpgpmin(igr))
   40 continue
!
      call groupcur(grsum,grsum0,gvsum,gvsum0,gcur0ka,gcurfka)
      izero = pgls + 3*pncoil
      n = 0
      if(kcycle.le.0) gcmin = 0._R8
      if(kcycle.le.0) gcmax = 0._R8
      do 42 i=izero+1,izero+2*pngroup,2
      n = n + 1
!
!.....note: Group 1 currents divided by 10
      fac1 = 1._R8
      if(n.eq.1) fac1 = 0.1_R8
      global(i) = grsum(n)*fac1
      global(i+1) = grsum0(n)*fac1
      if(grsum(n) .eq. 0 .and. gcurfka(n) .ne. 0)   then
      global(i) = gcurfka(n)
                                                    endif
!
!.....only include primary groups in  sum
      do 9912 nn=1,ncoil
      l = iabs(igroupc(nn))
      if(l.eq.n) go to 9913
 9912 continue
      go to 42
 9913 continue
      gcmin = min(gcmin,global(i),global(i+1))
      gcmax = max(gcmax,global(i),global(i+1))
   42 continue
      izero = pgls + 3*pncoil + 2*pngroup
!
      if(kcycle.le.0) then
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
      do 47 nn=1,ngroupt
      n=nogroupt(nn)
   47 engroup(n) = engrupi(n)
!  47 engroup(n) = engrupi(jnreg(n))
      endif
!
      n = 0
      sum1 = 0._R8
      sum2 = 0._R8
      do 43 i=izero+1,izero+pngroup
      n = n + 1
!......location pgls+3*pncoil+2*pngroup+i -- actual voltage in coil group i
      global(i) = gvsuma(n)
      gvmin = min(gvmin,gvsuma(n))
      gvmax = max(gvmax,gvsuma(n))
!......location pgls+3*pncoil+3*pngroup+i -- prepro voltage in coil group i
      global(i+pngroup) = gvsum0(n)
      gvmin0 = min(gvmin0,gvsum0(n))
      gvmax0 = max(gvmax0,gvsum0(n))
!......location pgls+3*pncoil+4*pngroup+i -- power in coil group i
      global(i+2*pngroup) = gvsuma(n)*grsum(n)
      gpmin = min(gpmin,global(i+2*pngroup))
      gpmax = max(gpmax,global(i+2*pngroup))
      engroup(n) = engroup(n) + global(i+2*pngroup)*dtglobe
!.(old)location pgls+3*pncoil+5*pngroup+i -- energy in coil group i
!     global(i+3*pngroup) = engroup(n)
      global(i+3*pngroup) = sumvsg(n)
      gemin = min(gemin,global(i+3*pngroup))
      gemax = max(gemax,global(i+3*pngroup))
      sum1 = sum1 + gvsuma(n)*grsum(n)
      sum2 = sum2 + engroup(n)
   43 continue
      izero = pgls + 3*pncoil + 6*pngroup
!.....location pgls+3*pncoil+6*pngroup+1 -- total power in coil system
      global(izero+1) = sum1
!.....location pgls+3*pncoil+6*pngroup+2 -- total energy in coil system
      global(izero+2) = sum2
      global(izero+3) = powerr
      global(izero+4) = poweri
      global(izero+5) = energyi
      global(izero+6) = energyr
      gpmint = min(gpmint,sum1,powerr,poweri)
      gpmaxt = max(gpmaxt,sum1,powerr,poweri)
      gemint = min(gemint,sum2,energyi,energyr)
      gemaxt = max(gemaxt,sum2,energyi,energyr)
!
!.....gap resistances, current, voltages
      izero = pgls + 3*pncoil + 6*pngroup+6
      do 451 iii=1,ngroupt
      ii=nogroupt(iii)
      i1 = izero + (ii-1)*3 + 1
      i2 = izero + (ii-1)*3 + 2
      i3 = izero + (ii-1)*3 + 3
      global(i1) = gapr(ii)
      global(i2) = gapi(ii)
      global(i3) = gapv(ii)
  451 continue
!
      izero = pgls + 3*pncoil + 9*pngroup+6
      if(itemp .le. 0) go to 46
      if(icirc.eq.0) call temprise
      if(kcycle.le.0) templmn = 0._R8
      if(kcycle.le.0) templmx = 1._R8
      n=0
      do 44 i=izero+1,izero+numpf
      n = n+1
      global(i) = templ(n)
      templmn = min(templmn,templ(n))
      templmx = max(templmx,templ(n))
   44 continue
   46 continue
!
      return
!
      entry globea
      if(kcycle .eq. 0 .or. (acoef(791).ne.0 .and. ipelav.eq.1))then
      do 98 j=2,npsi
      densum(j)=0._R8
   98 continue
      ipelav = 0
      write(nterm,8811) kcycle,times
 8811 format(" *** densum zeroed at cycle,time= ",i7,1pe12.4)
      endif
!
!
!.....keep a cumulative sum of density source as a function of psi
      do 97 j=2,npsi
      densum(j) = densum(j)+(srave(j)+sraveb(j)+simpe(j)+sraveedg(j)     &  
     &                     + sravejet(j))*dt
   97 continue
!
!.....calculates time averaged quantities between calls to global
!
      if(iplt2.ne.1 .and. iglobea.eq.1) go to 101
      tpergl = 0
      tave1 = 0._R8
      tave2 = 0._R8
      tave3 = 0._R8
      tave4 = 0._R8
      tave5 = 0._R8
      do 48 nn=1,ngroupt
      n=nogroupt(nn)
   48 gvsuma(n) = 0._R8
  101 continue
      iglobea = 1
      dtpergl = times-timedos
      tpergl = tpergl + dtpergl
!
      call groupcur(grsum,grsum0,gvsum,gvsum0,gcur0ka,gcurfka)
      do 50 nn=1,ngroupt
      n=nogroupt(nn)
   50 gvsuma(n) = gvsuma(n)+gvsum(n)*dtpergl
!..............................................
      tave1 = tave1 + hfluxp*dtpergl
!..............................................
      tave2 = tave2 - (paux+palpha)*dtpergl                              &  
     & -(global(119)+global(120)+global(121)+global(123)                 &  
     &               +global(124))*1.E6_R8*dtpergl
!..............................................
      sumw = 0._R8
      do 20 ii=1,nwire
      n = (ncoil-nwire) + ii
      facw = 1._R8
      if(isym.eq.1 .and. zwire(ii).gt.0) facw = 2._R8
      term = -(resave(ii)*tpi-rswire(ii)*ccoil(n))*facw
   20 sumw=sumw+term*udsv*ccoil(n)*udsi
      tave3 = tave3 + sumw*dtpergl
!..............................................
      tave4 = tave4 + enerpb*dtpergl
!..............................................
      sum = 0._R8
      fac = 1._R8
      if(isym.eq.1) fac = 2.0_R8
      do 92 i=iminn,imaxx
      do 92 j=jminn+1,jmaxx
   92 sum = sum - fac*(g(i,j)-gzero/xsqoj(i))
      sumdn = sum
      if(times.eq.timedos .or. sumdo.eq.0) go to 91
      tave5 = tave5 + (sumdn-sumdo)*tpi*gzero*udsp
   91 timedos = times
      sumdo = sumdn
      if(iplt2.le.nskip2) return
      if(tpergl.eq.0) return
      tave1 = tave1/tpergl
      tave2 = tave2/tpergl
      tave3 = tave3/tpergl
      tave4 = tave4/tpergl
      tave5 = tave5/tpergl
      do 49 nn=1,ngroupt
      n = nogroupt(nn)
   49 gvsuma(n) = gvsuma(n)/tpergl
      iglobea = 0
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
