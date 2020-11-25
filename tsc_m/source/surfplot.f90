      subroutine surfplot
!
!.....plots surface functions
!
      USE CLINAM
      USE SAPROP
      USE TCVCOM
      USE WALLCL
      USE EQRUNS
      USE NEWPLOT

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER npoints,imn,imx,i,n,l,ichimx
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 drmax,drmin,fmax,fmin,xfb,zfb,gradsq,dpsidx
      REAL*8 dpsidz,gsval,dlpsi,psval,dx0,dz0,r1,psixx,psixz,psizz
      REAL*8 r2,dx,dz,sgnum,xleft,xright,ps,pval,ppval,eval,epval
      REAL*8 rval,rpval,gval,gpval,gppval,pmax,pmin,tmax,tmin,vmax
      REAL*8 vmin,gmax,qmax,qmin,sigmax,sigmin,chimax,chimin,xmin
      REAL*8 xmax,chismall,ylabset,currfl,plcprnt,xlabset,alamda
      REAL*8 chimx,pos,deneln,a8p,a7p,a9,a10,a5,a17,a13,zpin,a2,a3
      REAL*8 a11,a12,a4,a6,a16,a14,a1,a15,a18,a19,a20,a21,a22,a23
      REAL*8 a24
!============
!     common/eqruns/findex
!     common /newplot/ chienca(ppsi),chiinca(ppsi),
!    1                 chiicopi(ppsi),chiecopi(ppsi),diffary(ppsi)
!
!     dimension gsumg(pngroup),grsumg(pngroup),gzsumg(pngroup),
!    1          grrsumg(pngroup),grho(pngroup),grho2(pngroup)
!     dimension pary(ppsi),gary(ppsi),psiy(ppsi),psih(ppsi),
!    1      eary(ppsi),reary(ppsi),qary(ppsi),teary(ppsi),
!    2      psit(ppsi),tiary(ppsi),div(2),viary(ppsi),rhary(ppsi),
!    3    dr(pnfob),dum(pnfob),dps(pnfob)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gzsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grrsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grho2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psiy
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psih
      REAL*8, ALLOCATABLE, DIMENSION(:) :: eary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: reary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: qary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: teary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psit
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tiary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: div
      REAL*8, ALLOCATABLE, DIMENSION(:) :: viary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rhary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dum
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dps
!============      
      IF(.not.ALLOCATED(gsumg)) ALLOCATE( gsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grsumg)) ALLOCATE( grsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(gzsumg)) ALLOCATE( gzsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grrsumg)) ALLOCATE( grrsumg(pngroup),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(grho)) ALLOCATE( grho(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grho2)) ALLOCATE( grho2(pngroup), STAT=istat)
      IF(.not.ALLOCATED(pary)) ALLOCATE( pary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(gary)) ALLOCATE( gary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(psiy)) ALLOCATE( psiy(ppsi), STAT=istat)
      IF(.not.ALLOCATED(psih)) ALLOCATE( psih(ppsi), STAT=istat)
      IF(.not.ALLOCATED(eary)) ALLOCATE( eary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(reary)) ALLOCATE( reary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(qary)) ALLOCATE( qary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(teary)) ALLOCATE( teary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(psit)) ALLOCATE( psit(ppsi), STAT=istat)
      IF(.not.ALLOCATED(tiary)) ALLOCATE( tiary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(div)) ALLOCATE( div(2), STAT=istat)
      IF(.not.ALLOCATED(viary)) ALLOCATE( viary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(rhary)) ALLOCATE( rhary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(dr)) ALLOCATE( dr(pnfob), STAT=istat)
      IF(.not.ALLOCATED(dum)) ALLOCATE( dum(pnfob), STAT=istat)
      IF(.not.ALLOCATED(dps)) ALLOCATE( dps(pnfob), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : surfplot  ' 
!============      
!
      if(dpsi.le.0.) return
      if(idata.eq.7) then
      npoints=nfob
      imn=1
      imx=nfob
      if(isym.eq.1) then
        npoints=nfob/2+1
        imn=1
        imx=npoints
        endif
        drmax=-100._R8
        drmin=+100._R8
        fmax=-100._R8
        fmin=+100._R8
        do 1 i=imn,imx
          xfb=xfob(i)
          zfb=zfob(i)
          if(isym.eq.1) zfb=-zfb
          call grap(0,zfb,xfb,gradsq,dpsidx,dpsidz,gsval,                    &  
     &    psval,psixz,psixx,psizz,1)
          dlpsi=psval-psilim
          if(gradsq.lt.1.E-12_R8) gradsq=1.E-12_R8
          if(abs(dlpsi).lt.1.E-12_R8) dlpsi=1.E-12_R8
          dx0=-dlpsi*dpsidx/gradsq
          dz0=-dlpsi*dpsidz/gradsq
          r1=-0.5_R8*dx0**2*psixx-dx0*dz0*psixz-0.5_R8*dz0**2*psizz
          r2=dz0**2*psizz-dx0*dz0*(psizz-psixx)-dx0**2*psixz
          dx=dx0+(r1*dpsidx+r2*dpsidz)/gradsq
          dz=dz0+(r1*dpsidz-r2*dpsidx)/gradsq
          sgnum=dlpsi/abs(dlpsi)
          dr(i-imn+1)=-sgnum*sqrt(dx**2+dz**2)
          dps(i-imn+1)=dlpsi/(psilim-psimin)
          dum(i-imn+1)=i
          fmax=max(fmax,dps(i-imn+1))
          fmin=min(fmin,dps(i-imn+1))
          drmax=max(drmax,dr(i-imn+1))
          drmin=min(drmin,dr(i-imn+1))
1       continue
        xleft=imn-1
        xright=imx+1
        call maps(xleft,xright,fmin,fmax,.142_R8,.800_R8,.700_R8,1.000_R8)    
        do 2 i=1,npoints
          call line(dum(i),0._R8,dum(i),dps(i))
2       continue
        call line(dum(1),0._R8,dum(npoints),0._R8)
        call setld(2._R8,50._R8,0,0,2,1)
        write(s100,1043)
        call gtext(s100,80,0)
1043    format("flux errors")
        if(drmax.gt.0.10_R8) drmax=0.10_R8
        if(drmin.lt.-0.10_R8) drmin=-0.10_R8
        if(drmin.gt.-0.02_R8) drmin=-0.0201_R8
        if(drmax.lt.+0.02_R8) drmax=+0.0201_R8
        call maps(xleft,xright,drmin,drmax,.142_R8,.800_R8,.290_R8,        &  
     &   .584_R8)
        do 3 i=1,npoints
          call line(dum(i),0._R8,dum(i),dr(i))
3       continue
        call linep(dum(1),0.01_R8,dum(npoints),0.01_R8,-1)
        call linep(dum(1),0.02_R8,dum(npoints),0.02_R8,-1)
        call linep(dum(1),-0.01_R8,dum(npoints),-0.01_R8,-1)
        call linep(dum(1),-0.02_R8,dum(npoints),-0.02_R8,-1)
        call line(dum(1),0._R8,dum(npoints),0._R8)
        call setld(2._R8,22._R8,0,0,2,1)
        write(s100,1044)
        call gtext(s100,80,0)
1044    format("pos errors")
        call setld(35._R8,4._R8,0,0,2,0)
        write(s100,1045) times
        call gtext(s100,80,0)
1045    format(" time (sec)  ",f10.2)
        call frscj(6)
      endif
!
 
      npts = nx/2
      if(npsit.gt.1) npts = npsit-1
!
      chimin = 1.e-3_R8
      do 100 n=1,npts+1
        ps = xsv2(n)
        call peval(ps,2,pval,ppval,imag,jmag)
        if(acoef(901).eq.0._R8) call eeval(ps,2,eval,epval,imag,jmag)
        if(acoef(901).eq.0._R8) call reval(ps,idens,1,rval,rpval,imag,jmag)
        call geval(ps,2,gval,gpval,gppval,imag,jmag)
        pary(n) = pval*udsp*1.E-3_R8
        if(acoef(901).eq.0._R8) eary(n) = eval*udsp*1.E-3_R8
        if(acoef(901).eq.0._R8) reary(n) = rval*udsd
        if(acoef(901).eq.0._R8) viary(n) = anhy(n)
        if(acoef(901).eq.0._R8) rhary(n) = anhe(n)
        gary(n) = gval
        qary(n) = qprof2(n)
        if(acoef(901).eq.0._R8) teary(n) = .5_R8*(te(n) + te(n+1))
        if(acoef(901).eq.0._R8) tiary(n) = .5_R8*(ti(n) + ti(n+1))
!
! 100 psit(n) = (xsv2(n)-psimin)*tpi
! 100 psit(n) = rminora(n)
  100 psit(n) = sqrt((float(n-1)*dpsi)/(float(npsit-1)*dpsi))

!
      pmax = pary(npts+1)
      pmin = pary(npts+1)
      if(acoef(901).eq.0._R8) tmax = teary(npts+1)
      if(acoef(901).eq.0._R8) tmin = teary(npts+1)
      if(acoef(901).eq.0._R8) vmax = reary(npts+1)
      if(acoef(901).eq.0._R8) vmin = reary(npts+1)
      gmax = gary(npts+1)
      gmin = gary(npts+1)
      qmax = max(qary(1),qary(npts+1))
      qmin = min(qary(1),qary(npts+1))
      if(acoef(901).eq.0._R8) sigmax = -1.E6_R8
      if(acoef(901).eq.0._R8) sigmin=+1.E6_R8
      if(acoef(901).eq.0._R8) chimax=-1.E6_R8
      xmin = 0._R8
      xmax = psit(npts+1)
      if(xmax.le.xmin) xmax = xmin + 1.
      do 200 n=2,npts
      pmax = max(pmax,pary(n))
      pmin = min(pmin,pary(n))
      if(acoef(901).eq.0._R8) tmax = max(tmax,teary(n))
      if(acoef(901).eq.0._R8) tmin = min(tmin,teary(n))
      if(acoef(901).eq.0._R8) tmax = max(tmax,tiary(n))
      if(acoef(901).eq.0._R8) tmin = min(tmin,tiary(n))
!     CMOD
!     qmax = 5.0
!     qmin = 0.5
      qmax = max(qmax,qary(n))
      qmin = min(qmin,qary(n))
      if(acoef(901).eq.0._R8) vmax = max(vmax,reary(n),viary(n),rhary(n)  &  
     & )
      if(acoef(901).eq.0._R8) vmin = min(vmin,reary(n),viary(n),rhary(n)  &  
     & )
      gmax = max(gmax,gary(n))
      gmin = min(gmin,gary(n))
      if(acoef(901).eq.0._R8) sigmax=max(sigmax,sighot(n))
      if(acoef(901).eq.0._R8) sigmin=min(sigmin,sighot(n))
      if(acoef(901).eq.0._R8) chimax=max(chimax,chiinca(n),chiitima(n),chiicopi(n), &
                                                chienca(n),chietema(n),chiecopi(n))
  200 continue
      chismall = 1.E-3_R8
      if(chimin.lt.chismall) chimin = chismall
      do n=1,npts
      if(itrmod.eq.2.or.(itrmod.ge.8.and.itrmod.le.11).or.itrmod.ge.16) then
        if(chiicopi(n).lt. chismall) chiicopi(n) = chismall
        if(chiecopi(n).lt. chismall) chiecopi(n) = chismall
      endif
      if(itrmod.eq.6.or.itrmod.eq.7.or.itrmod.eq.11.or.itrmod.eq.26) then
        if(chienca(n).le.chismall) chienca(n)=chismall
        if(chiinca(n).le.chismall) chiinca(n)=chismall
      endif
      if((itrmod.ge.6.and.itrmod.le.9).or.itrmod.eq.26) then
        if(chiitima(n).lt.chismall) chiitima(n) = chismall
        if(chietema(n).lt.chismall) chietema(n) = chismall
      endif
      enddo
!
      if(pmax.le.pmin) pmax = pmin+1._R8
      if(acoef(901).eq.0._R8.and.tmax.le.tmin) tmax = tmin+1._R8
      if(qmax.le.qmin) qmax = qmin+1._R8
      if(acoef(901).eq.0._R8.and.vmax.le.vmin) vmax = vmin+1._R8
!
      if(acoef(901) .eq.0) then
        call mapg(xmin,xmax,pmin,pmax,.1_R8,.33_R8,.74_R8,.97_R8)
      else
        call mapg(xmin,xmax,pmin,pmax,.1_R8,.45_R8,.55_R8,0.90_R8)
      endif
      call trace(psit,pary,npts+1,-1,-1,0._R8,0._R8)
      if(acoef(901).eq.0._R8) call trace(psit,eary,npts+1,-1,-1,0._R8,   &  
     & 0._R8)
      ylabset = pmax+(pmax-pmin)*.05_R8
      call setold(xmin+(xmax-xmin)*.1_R8,ylabset,1,0,1,0)
      if(acoef(901).eq.0._R8) write(s100,1001)
      if(acoef(901).ge.1._R8) write(s100,2001)
      call gtext(s100,80,0)
 1001 format("p,e (kpa) vs sqrt")
 2001 format("  p[kN/m^2] vs sqrt")
      if(acoef(901).eq.0._R8) then
        call mapg(xmin,xmax, tmin,tmax,.43_R8,.66_R8,.74_R8,.97_R8)
        call tracec(1he,psit,teary,npts+1,-1,-1,0._R8,0._R8)
        call tracec(1hi,psit,tiary,npts+1,-1,-1,0._R8,0._R8)
        call setold(xmin+(xmax-xmin)*.07_R8,tmax+(tmax-tmin)*.05_R8,1,0,1,  &  
     &     0)
        write(s100,1002)
        call gtext(s100,80,0)
 1002   format("te and ti vs sqrt")
      endif
      if(acoef(901).eq.0._R8) then
        call mapg(xmin,xmax,qmin,qmax,.09_R8,.32_R8,.38_R8,.61_R8)
        call tracec(1hq,psit,qary,npts+1,-1,-1,0._R8,0._R8)
      else
        call mapg(xmin,xmax,qmin,qmax,.6_R8,.95_R8,.55_R8,0.90_R8)
        call trace(psit,qary,npts+1,-1,-1,0._R8,0._R8)
      endif
      call setold(xmin+(xmax-xmin)*.07_R8,qmax+(qmax-qmin)*.05_R8,1,0,1,  &  
     & 0)
      if(acoef(901).eq.0._R8) write(s100,1003)
      if(acoef(901).ge.1._R8) write(s100,2003)
      call gtext(s100,80,0)
 1003 format("q and g vs sqrt")
 2003 format("       q vs sqrt")
!
      fluxlink = 0._R8
      call fieldg(xplas,zplas,gsumg,grsumg,gzsumg,                       &  
     &            grrsumg,grho,grho2)
      do 1169 l=1,pngroup
      currfl = 0._R8
      do 297 n=1,ntpts
  297 currfl = currfl + fact(n)*gcur(n,l)*usdi
 1169 fluxlink = fluxlink+gsumg(l)*tpi*(gcurfb(l)+currfl)
      plcprnt = tcurdtp*tpi*udsi
      xlabset = xmin-(xmax-xmin)*.30_R8
      if(acoef(901).eq.0._R8)                                            &  
     &   call setold(xlabset,qmin-.35_R8*(qmax-qmin),1,0,1,0)
      if(acoef(901).ge.1._R8)                                            &  
     &   call setch(1._R8,15.1_R8,0,1,0)
      if(acoef(901).eq.0._R8)                                            &  
     &   write(s100,1005) kcycle,times,global(42),plcprnt
 1005 format("  kcycle=",i7,"   time=",1pe13.5                           &  
     &      ,"  vsec=",1pe12.4,"  ip=",1pe12.4)
      if(acoef(901).ge.1._R8) then
      write(s100,2005) gzero/xplas,plcprnt, fluxlink
      do i=1,5
      write(nterm,*) " "
      enddo
      write(nterm,2005)  gzero/xplas,plcprnt, fluxlink
 2005 format(" Bt[T]=",1pe10.3,"  Ip[A]  =",1pe10.3,                     &  
     &   " Flux Linkage[W]=",1pe10.3)
      endif
      call gtext(s100,80,0)
!
      alamda = ali2+betapol
      if(acoef(901).eq.0._R8)                                            &  
     &   write(s100,167) ali2,betapol,vol,uint,npsit,beta
  167 format(" li2,bp,vol,uint,npsit,beta",1p4e9.2,i4,0pf10.6)
      if(acoef(901).ge.1._R8) then
      write(s100,2167) 2._R8*ali2,betapol,beta*100._R8
 2167 format(" Li(3)=",1pe10.3,"  Betapol=",1pe10.3,                     &  
     &   " Beta[%]        =",1pe10.3)
      write(nterm,2167) 2._R8*ali2,betapol,beta*100._R8
      call gtext(s100,80,0)
      endif
!
      if(acoef(901).eq.0._R8)                                            &  
     &   write(s100,168) rmajor,rminor,shape3,shape5,shape6,             &  
     &   shape7,ain99
  168 format(" r0,a,delt,eps,xsep,zsep",6f7.3,4h(in=,f4.2,1h) )
      if(acoef(901).ge.1._R8) then
      write(s100,2168) rmajor,rminor
 2168 format(" R0[m]=   ",f7.3,"  a[m]   =   ",f7.3)
      write(nterm,2168) rmajor,rminor
      call gtext(s100,80,0)
!
      write(s100,3168) shape5,shape3
      write(nterm,3168) shape5,shape3
 3168 format(" k    =   ",f7.3,"  d      =   ",f7.3)
      endif
      call gtext(s100,80,0)
!
      if(acoef(901).ge.1._R8) then
      write(s100,21680) el95,del95
21680 format(" k95  =   ",f7.3,"  d95    =   ",f7.3)
      write(nterm,21680) el95,del95
      call gtext(s100,80,0)
!
      write(s100,2169) shape6,shape7
      write(nterm,2169) shape6,shape7
 2169 format(" Rx[m]=   ",f7.3,"  Zx[m]  =   ",f7.3)
      call gtext(s100,80,0)
      endif
!
      if(acoef(901).eq.0._R8) then
      write(s100,169) aliga,del95,el95,ain95
  169 format(" li(GA def)=",f5.2,"      del95,eps95    ",2f7.3,          &  
     &         "  (in-95=",f4.2,1h) )
      call gtext(s100,80,0)
      write(s100,170) global(61),global(15),fluxlink
  170 format(" q(95%)  q(100%)  ",2f7.3," Webers=",1pe10.2)
      endif
!
      if(acoef(901).ge.1._R8) then
      write(s100,2170) qprof2(1),global(15)
      write(nterm,2170) qprof2(1),global(15)
 2170 format(" q0   =   ",f7.3,"  q100   =   ",f7.3)
      call gtext(s100,80,0)
      write(s100,3170) global(61),global(62)
      write(nterm,3170) global(61),global(62)
 3170 format(" q95  =   ",f7.3,"  qstar  =   ",f7.3)
      call gtext(s100,80,0)
!      write(s100,3171) findex
! 3171 format(" nindx=   ",f7.3)
!      call gtext(s100,80,0)
 
      write(nsc1,4180)
 4180 format(" p and q vs sqrt")
      call frscj(6)
      endif
 
!
      if(acoef(901).eq.0._R8) then
        div(1) = gmin
        div(2) = gmax
        call map(xmin,xmax,gmin,gmax,.09_R8,.32_R8,.38_R8,.61_R8)
        call tracec(1hg,psit,gary,npts+1,-1,-1,0._R8,0._R8)
        call gaxisf(xmax,gmin,xmax,gmax,0,1,0,'f5.2',4,div)
        write(nsc1,1088) kcycle
 1088   format(" surface profiles,cycle=",i7)
        call mapg(xmin,xmax,vmin,vmax,.77_R8,1._R8,.38_R8,.61_R8)
        call tracec(1he,psit,reary,npts+1,-1,-1,0._R8,0._R8)
        call tracec(1hi,psit,viary,npts+1,-1,-1,0._R8,0._R8)
        call tracec(1hh,psit,rhary,npts+1,-1,-1,0._R8,0._R8)
        call setold(xmin+(xmax-xmin)*.1_R8,vmax+(vmax-vmin)*.05_R8,1,0,1,0)
        write(s100,1004)
        call gtext(s100,80,0)
 1004   format("n vs sqrt")
!
        if(kcycle.gt.0.and.ipres.eq.0.and.chimax.gt.0) then
          ichimx=log10(chimax)
          chimx=10._R8**(1._R8+ichimx)
          call mapgsl(xmin,xmax,chimin,chimax,.77_R8,1.0_R8,.74_R8,.97_R8)
          if(itrmod .ne. 6 .and. itrmod .ne. 26) call tracec(1hT,psit,chiesec,npts+1,-1,-1,0._R8,0._R8)
          if(itrmod.eq.6.or.itrmod.eq.7.or.itrmod.eq.9 .or.itrmod.eq.11.or.itrmod.eq.26)  &
            call tracec(1hN,psit,chienca,npts+1,-1,-1,0._R8,0._R8)
          if(itrmod.eq.6.or.(itrmod.ge.8.and.itrmod.le.11).or.itrmod.eq.26)               &
          call tracec(1hG,psit,chietema,npts+1,-1,-1,0._R8,0._R8)
!         if(itrmod.eq.2 .or. (itrmod.ge.8.and.itrmod.le.11))                  &
          if(itrmod.eq.2 .or. (itrmod.ge.8.and.itrmod.le.11).or.itrmod.ge.16)                  &
          call tracec(1hC,psit,chiecopi,npts+1,-1,-1,0._R8,0._R8)
          call setold(xmin+(xmax-xmin)*.07_R8,chimx+(chimx-chimin)*.10_R8,1,  &  
          0,1,0)
          write(s100,1007)
          call gtext(s100,80,0)
1007      format("chiE(m**2/sec)")
        endif
!
        if(kcycle.gt.0.and.ipres.eq.0.and.chimax.gt.0) then
          ichimx=log10(chimax)
          chimx=10**(1._R8+ichimx)
          call mapgsl(xmin,xmax,chimin,chimax,.45_R8,.68_R8,.38_R8,.61_R8)
          if(itrmod .ne. 6 .and. itrmod .ne. 26) call tracec(1hT,psit,chiisec,npts+1,-1,-1,0._R8,0._R8)
          if(itrmod.eq.6.or.itrmod.eq.7.or.itrmod.eq.9 .or.itrmod.eq.11.or.itrmod.eq.26)  &
!                                                                        &  
          call tracec(1hN,psit,chiinca,npts+1,-1,-1,0._R8,0._R8)
          if(itrmod.eq.6.or.(itrmod.ge.8.and.itrmod.le.11).or.itrmod.eq.26)               &
          call tracec(1hG,psit,chiitima,npts+1,-1,-1,0._R8,0._R8)
!         if(itrmod.eq.2 .or. (itrmod.ge.8.and.itrmod.le.11))                                 &
          if(itrmod.eq.2 .or. (itrmod.ge.8.and.itrmod.le.11).or.itrmod.ge.16)              &
          call tracec(1hC,psit,chiicopi,npts+1,-1,-1,0._R8,0._R8)
          call setold(xmin+(xmax-xmin)*.07_R8,chimx+(chimx-chimin)*.10_R8,1,  &  
             0,1,0)
          write(s100,1008)
          call gtext(s100,80,0)
1008      format("chiI(m**2/sec)")
        endif
        call frscj(6)
        if(acoef(3101).gt.0) call tempplot
!>>>>>
!>>>>> DEBUG
!       do i=1,npts+1
!       write(73,1073) i,chiesec(i),chienca(i),chietema(i),chiisec(i),chiinca(i),chiitima(i)
!1073   format(i5,1p6e12.4)
!       enddo
!       stop '73'
      endif
 
        pos=global(25)*global(37)*tcurdtp*global(77)*global(34)*global(35)     
      if(pos .le. 0) return
!
!..set up call to popcon if ialpha>0
!
      if(global(78).eq.0.or.global(14).eq.0.or.global(29).eq.0)          &  
     &  go to 201
      if(ialpha.eq.1.and.kcycle.gt.nskipsf) then
      deneln=global(77)*1.E20_R8
      a8p=global(35)
      a7p=global(34)
      a9=global(37)
      if(iplim.lt.0) a9=global(65)
      a10=global(36)
      if(iplim.lt.0) a10=global(64)
      a5=tcurdtp*tpi*udsi
      a17=acoef(55)
      a13=global(135)
      zpin=global(45)*1.E6_R8
      a2=global(79)/global(78)-1._R8
      a3=global(29)-1._R8
      a11=global(78)*1.E20_R8
!     a12=global(27)/global(29)*(1.+a3)*(1.+a2)/(1.+a2+a3)
!     a12=a12/1.e3
      a12=global(27)/1.E3_R8
      a4=global(63)/global(14)-1._R8
      a6=gzero/a7p
      a16=0.5_R8
      a14=2.5_R8
!
      if(a14 .le. 0.or.a9 .le. 0.or.a5 .le. 0.or.deneln.le.0             &  
     &             .or.a6 .le. 0.or.a8p.le. 0 .or.a7p.le.0)go to 201
      if(pohmic+paux+prad+palpha .le. 0) go to 201
!
      a1=(global(25)*1.E-3_R8)/(.0521_R8*sqrt(a14)*a9**0.25_R8*(a5/      &  
     & 1.E6_R8)**.85_R8                                                  &  
     &   *(deneln/1.e19_R8)**0.1_R8*a6**0.3_R8*a8p**0.3_R8               &  
     &   *a7p**0.85_R8*sqrt(1.e6_R8/(pohmic+paux+prad+palpha)))
      a15=zimp
!
      a18=kcycle
      a19=times
      a20=pohmic
      a21=palpha
      a22=-global(119)*1.E6_R8
      a23=-global(120)*1.E6_R8
      a24=global(80)*1.E6_R8
      call popcon(a1,a2,a3,a4,a5,a6,a7p,a8p,a9,a10,a11,a12,a13,a14,a15,  &  
     &   a16,a17,a18,a19,a20,a21,a22,a23,a24)
      endif
  201 continue
!
 1786 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
