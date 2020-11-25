      subroutine surfplot3
!.....6.95 shape
!.....6.95 shape
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
      INTEGER n,l,i,ichimx
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ps,pval,ppval,eval,epval,rval,rpval,gval,gpval
      REAL*8 gppval,pmax,pmin,tmax,tmin,vmax,vmin,gmax,qmax,qmin
      REAL*8 sigmax,sigmin,chimax,chimin,xmin,xmax,chismall,ylabset
      REAL*8 currfl,plcprnt,xlabset,alamda,chimx,pos
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
      if (istat .ne. 0) stop 'Allocation Error : surfplot3  ' 
!============      
!
!
!
!......disabled 1/18/10 as it is redundant
       return
 
      if(dpsi.le.0) return
      npts = nx/2
      if(npsit.gt.1) npts = npsit-1
!
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
!     CMOD
!     qmax = 5.0
!     qmin = 0.5
      if(acoef(901).eq.0._R8) sigmax = -1.E6_R8
      if(acoef(901).eq.0._R8) sigmin=+1.E6_R8
      if(acoef(901).eq.0._R8) chimax=-1.E6_R8
      if(acoef(901).eq.0._R8) chimin=+1.E6_R8
      xmin = 0._R8
      xmax = psit(npts+1)
      do 200 n=2,npts
      pmax = max(pmax,pary(n))
      pmin = min(pmin,pary(n))
      if(acoef(901).eq.0._R8) tmax = max(tmax,teary(n))
      if(acoef(901).eq.0._R8) tmin = min(tmin,teary(n))
      if(acoef(901).eq.0._R8) tmax = max(tmax,tiary(n))
      if(acoef(901).eq.0._R8) tmin = min(tmin,tiary(n))
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
      if(acoef(901).eq.0._R8) chimax=max(chimax,chiesec(n),chiisec(n))
      if(acoef(901).eq.0._R8) chimin=min(chimin,chiesec(n),chiisec(n))
  200 continue
      chismall = 1.E-8_R8
      if(chimin.lt.chismall) chimin = chismall
      do n=1,npts
      if(itrmode.eq.2.or.(itrmode.ge.8.and.itrmode.le.11)) then
        if(chiicopi(n).lt. chismall) chiicopi(n) = chismall
        if(chiecopi(n).lt. chismall) chiecopi(n) = chismall
      endif
      if(itrmode.eq.6.or.itrmode.eq.7.or.itrmode.eq.11.or.itrmode.eq.26) then
        if(chienca(n).le.chismall) chienca(n)=chismall
        if(chiinca(n).le.chismall) chiinca(n)=chismall
      endif
      if((itrmode.ge.6.and.itrmode.le.9).or.itrmode.eq.26) then
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
 1001 format("p,e (kpa) vs <r>")
 2001 format("  p[kN/m^2] vs <r>")
      if(acoef(901).eq.0._R8) then
      call mapg(xmin,xmax, tmin,tmax,.43_R8,.66_R8,.74_R8,.97_R8)
      call tracec(1he,psit,teary,npts+1,-1,-1,0._R8,0._R8)
      call tracec(1hi,psit,tiary,npts+1,-1,-1,0._R8,0._R8)
      call setold(xmin+(xmax-xmin)*.07_R8,tmax+(tmax-tmin)*.05_R8,1,0,1,  &  
     & 0)
      write(s100,1002)
      call gtext(s100,80,0)
 1002 format("te and ti vs <r>")
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
 1003 format("q and g vs <r>")
 2003 format("       q vs <r>")
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
 4180 format(" p and q vs <r>")
      call frscj(6)
      endif
 
!
      if(acoef(901).eq.0._R8) then
      div(1) = gmin
      div(2) = gmax
      call map(xmin,xmax,gmin,gmax,.09_R8,.32_R8,.38_R8,.61_R8)
      call tracec(1hg,psit,gary,npts+1,-1,-1,0._R8,0._R8)
      call gaxisf(xmax,gmin,xmax,gmax,0,1,0,'f5.2',4,div)
      write(nsc1,1088) times
 1088 format(" profiles vs <r>,time=",f8.3)
      call mapg(xmin,xmax,vmin,vmax,.77_R8,1._R8,.38_R8,.61_R8)
      call tracec(1he,psit,reary,npts+1,-1,-1,0._R8,0._R8)
      call tracec(1hi,psit,viary,npts+1,-1,-1,0._R8,0._R8)
      call tracec(1hh,psit,rhary,npts+1,-1,-1,0._R8,0._R8)
      call setold(xmin+(xmax-xmin)*.1_R8,vmax+(vmax-vmin)*.05_R8,1,0,1,  &  
     & 0)
      write(s100,1004)
      call gtext(s100,80,0)
 1004 format("n vs <r>")
!
      do n=2, npts
      chimin = min(chimin,chienca(n),chiinca(n))
      enddo
      chimin = max(chismall, chimin)
!
      if(kcycle.gt.0.and.ipres.eq.0.and.chimax.gt.0) then
      ichimx=log10(chimax)
      chimx=10**(1._R8+ichimx)
      call mapgsl(xmin,xmax,chimin,chimax,.77_R8,1.0_R8,.74_R8,.97_R8)
      call tracec(1hT,psit,chiesec,npts+1,-1,-1,0._R8,0._R8)
      if(itrmode.eq.6.or.itrmode.eq.7.or.itrmode.eq.9 .or.itrmode.eq.11.or.itrmode.eq.26)  &  
!    &                                                                   &  
     &call tracec(1hN,psit,chienca,npts+1,-1,-1,0._R8,0._R8)
!     if(itrmode.eq.6.or.(itrmode.ge.8.and.itrmode.le.11))               &  
      if(itrmode.eq.6.or. itrmod.eq.14 .or.                              &
     & (itrmode.ge.8.and.itrmode.le.11).or.itrmode.eq.26)                                 &
     &call tracec(1hG,psit,chietema,npts+1,-1,-1,0._R8,0._R8)
!     if(itrmode.ge.8.and.itrmode.le.11)                                 &  
      if((itrmode.ge.8.and.itrmode.le.11).or. itrmod .eq. 14)            &
     &call tracec(1hC,psit,chiecopi,npts+1,-1,-1,0._R8,0._R8)
      call setold(xmin+(xmax-xmin)*.07_R8,chimx+(chimx-chimin)*.10_R8,1,  &  
     & 0,1,0)
      write(s100,1007)
      call gtext(s100,80,0)
1007  format("chiE(m**2/sec)")
      endif
!
      if(kcycle.gt.0.and.ipres.eq.0.and.chimax.gt.0) then
      ichimx=log10(chimax)
      chimx=10**(1._R8+ichimx)
      chimin=min(chimin,chiinca(n),chiisec(n))
      call mapgsl(xmin,xmax,chimin,chimax,.45_R8,.68_R8,.38_R8,.61_R8)
      call tracec(1hT,psit,chiisec,npts+1,-1,-1,0._R8,0._R8)
      if(itrmode.eq.6.or.itrmode.eq.7.or.itrmode.eq.9 .or.itrmode.eq.11.or.itrmode.eq.26)  &  
!    &                                                                   &  
     &call tracec(1hN,psit,chiinca,npts+1,-1,-1,0._R8,0._R8)
      if(itrmode.eq.6.or.itrmod.eq.14.or.                                &
     & (itrmode.ge.8.and.itrmode.le.11).or.itrmode.eq.26)                                 &
!     if(itrmode.eq.6.or.(itrmode.ge.8.and.itrmode.le.11))               &  
     &call tracec(1hG,psit,chiitima,npts+1,-1,-1,0._R8,0._R8)
!     if(itrmode.ge.8 .and. itrmode.le.11)                               &  
      if((itrmode.ge.8 .and. itrmode.le.11) .or. itrmod .eq. 14)         &
     &call tracec(1hC,psit,chiicopi,npts+1,-1,-1,0._R8,0._R8)
      call setold(xmin+(xmax-xmin)*.07_R8,chimx+(chimx-chimin)*.10_R8,1,  &  
     & 0,1,0)
      write(s100,1008)
      call gtext(s100,80,0)
1008  format("chiI(m**2/sec)")
      endif
      call frscj(6)
      endif
 
      pos=global(25)*global(37)*tcurdtp*global(77)*global(34)*global(35)     
      if(pos .le. 0) return
!
!..set up call to popcon if ialpha>0
!
      if(global(78).eq.0 .or. global(14).eq.0 .or. global(29).eq.0)      &  
     &  go to 201
  201 continue
!
 1786 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
