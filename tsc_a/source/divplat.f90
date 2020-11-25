      subroutine divplat
!
!.....compute strike point locations on divertor plate and
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ifrsdp,ios23,l,n,i,iaxis,icount,iv,ia,istrike,ilook
      INTEGER itry,jtry,ir,jr,ii,iseg1,jseg1,iseg2,jseg2,isw,irecdv
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 formf,ais,xsegah,zsegah,asegah,psega,dpsega,xinter
!     REAL*8 zinter,dpsdxh,dpsdzh,disi,xv,psiv,xa,psia,pinterp
      REAL*8 disi,xv,psiv,xa,psia,pinterp
      REAL*8 pgrad,pgradin,pscrapin,dissav,xtry,ztry,ptry,aisegt,xr
      REAL*8 zr,psir,aisegr,xave,zave,disave,aisegtt,aisegrt,aigt
      REAL*8 ailt,dtry,xc,zc,disc,sum1,sum2,xseg1,zseg1,pseg1,xseg2
      REAL*8 zseg2,pseg2,pseg,aih,fpln1,fpln2,area,term1,term2
      REAL*8 gradsq,gsval,psval,psixz,psixx,psizz,plossmw
!============
      data  ifrsdp / 1 /
!
!.....heat deposition profile
!
!     dimension formf(pnseg),ais(2),xsegah(pnplat,pnseg),
!    1          zsegah(pnplat,pnseg),asegah(pnplat,pnseg),
!    2          psega(pnplat,pnseg),dpsega(pnplat,pnseg),
!    3          xinter(2,pnplat),zinter(2,pnplat),
!    4          dpsdxh(pnplat,pnseg),dpsdzh(pnplat,pnseg)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: formf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ais
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xsegah
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: zsegah
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: asegah
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: psega
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: dpsega
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xinter
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: zinter
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: dpsdxh
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: dpsdzh
!============      
      IF(.not.ALLOCATED(formf)) ALLOCATE( formf(pnseg), STAT=istat)
      IF(.not.ALLOCATED(ais)) ALLOCATE( ais(2), STAT=istat)
      IF(.not.ALLOCATED(xsegah)) ALLOCATE( xsegah(pnplat,pnseg),         &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(zsegah)) ALLOCATE( zsegah(pnplat,pnseg),         &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(asegah)) ALLOCATE( asegah(pnplat,pnseg),         &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(psega)) ALLOCATE( psega(pnplat,pnseg),           &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(dpsega)) ALLOCATE( dpsega(pnplat,pnseg),         &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(xinter)) ALLOCATE( xinter(2,pnplat), STAT=istat)
      IF(.not.ALLOCATED(zinter)) ALLOCATE( zinter(2,pnplat), STAT=istat)
      IF(.not.ALLOCATED(dpsdxh)) ALLOCATE( dpsdxh(pnplat,pnseg),         &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(dpsdzh)) ALLOCATE( dpsdzh(pnplat,pnseg),         &  
     &                                 STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : divplat  ' 
!============      
!
      if(iplate .eq. 0) return
      if(nplate .eq. 0) return
      if(lrswtch.ne.0) return
      if(ineg .ne. 0) return
!
      if(ifrsdp .ne. 1) go to 501
      ifrsdp = 0
!
!.....divertor plate initialization goes here
!
!     ifiledi(1:6) = 'divhis'
!     ifiledi(7:7) = isuffix(1:1)
      if( numargs .lt. 1 ) then
         ifiledi = 'divhis' // isuffix(1:1)
      else
         ifiledi = 'divhis' // '.' // trim(suffix)
      end if
      open(ndivlu,file=trim(ifiledi),status='unknown',iostat=ios23)
      call tversion(ndivlu)
      write(ndivlu,6000) (name(l),l=1,8)
 6000 format(10a8)
      write(ndivlu,6002) nplate,(nseg(n),fplate(n,1),fplate(n,2),        &  
     &                                 n=1,nplate)
 6002 format(1x,i5,8(/1x,i5,1p2e12.4))
      do 601 n=1,nplate
      do 601 i=1,nseg(n)
      xsegah(n,i) = .5_R8*(xsega(n,i)+xsega(n,i+1))
      zsegah(n,i) = .5_R8*(zsega(n,i)+zsega(n,i+1))
      disi =                                                             &  
     & sqrt((xsega(n,i+1)-xsega(n,i))**2+(zsega(n,i+1)-zsega(n,i))**2)
      asegah(n,i) = tpi*xsegah(n,i)*disi
  601 continue
      write(ndivlu,6001) ((xsegah(n,i),i=1,nseg(n)),n=1,nplate)
      write(ndivlu,6001) ((zsegah(n,i),i=1,nseg(n)),n=1,nplate)
      write(ndivlu,6001) ((asegah(n,i),i=1,nseg(n)),n=1,nplate)
 6001 format(1x,1p9e14.6)
!
!
  501 continue
!
!...flux gradient at outer midplane
      iaxis = int((xmag-ccon)/deex) + 2
      icount = 0
      xv = xmag
      iv = iaxis
      psiv = psimin
  705 xa = xv
      ia = iv
      psia = psiv
      xv = xv + deex
      iv = iv + 1
      psiv = pinterp(xv,zmag,iv,jmag)
      icount = icount + 1
      if(icount.ge.nx) return
      if(psiv.lt.psilim) go to 705
      pgrad = abs((psiv-psia)/(xv-xa))
!
!...flux gradient at outer midplane
      iaxis = int((xmag-ccon)/deex) + 2
      icount = 0
      xv = xmag
      iv = iaxis
      psiv = psimin
  706 xa = xv
      ia = iv
      psia = psiv
      xv = xv - deex
      iv = iv - 1
      psiv = pinterp(xv,zmag,iv,jmag)
      icount = icount + 1
      if(icount.ge.nx) return
      if(psiv.lt.psilim) go to 706
      pgradin = abs((psiv-psia)/(xv-xa))
!
!...flux increment associated with 0.6cm on outer midplane
      pscrape = .006_R8*pgrad
!...flux increment associated with 0.6cm on inner midplane
      pscrapin = .006_R8*pgradin
!
!.....first calculate strike points using binary search
!
      do 101 n=1,nplate
!
!.....search for two strike points on each plate
      do 100 istrike=1,2
      ilook = 0
      dissav = 1.E8_R8
!
      if(iplim.gt.0) go to 200
      go to (201,202),istrike
  201 xtry = xsega(n,1)
      ztry = zsega(n,1)
      itry = (xtry-ccon)/deex + 2.0_R8
      jtry = (ztry-zzero)/deez + 2.0_R8
      ptry = pinterp(xtry,ztry,itry,jtry)
      aisegt = 1._R8
      do 210 i=2,nseg(n)+1
      xr = xsega(n,i)
      zr = zsega(n,i)
      ir = (xr-ccon)/deex + 2.0_R8
      jr = (zr-zzero)/deez + 2.0_R8
      psir = pinterp(xr,zr,ir,jr)
      aisegr = i
      if(ptry.ge.psisep .and. psir.lt.psisep) go to 209
      ptry = psir
      xtry = xr
      ztry = zr
      aisegt = aisegr
      go to 210
  209 ilook = 1
      xave = .5_R8*(xr+xtry)
      zave = .5_R8*(zr+ztry)
      disave = sqrt((xave-xsep(1))**2+(zave-zsep(1))**2)
      if(    (zave*zsep(1).lt.0  .or. xsep(1).le.0) .and. xsep(2).gt.0)  &  
     &    disave=sqrt((xave-xsep(2))**2+(zave-zsep(2))**2)
      if(disave.gt.dissav) go to 210
      aisegtt = aisegt
      aisegrt = aisegr
      dissav = disave
  210 continue
      go to(200,150),ilook+1
  202 xr = xsega(n,nseg(n)+1)
      zr = zsega(n,nseg(n)+1)
      ir = (xr-ccon)/deex + 2.0_R8
      jr = (zr-zzero)/deez + 2.0_R8
      psir = pinterp(xr,zr,ir,jr)
      aisegr = nseg(n)+1
      do 220 ii=1,nseg(n)
      i = 1+nseg(n)-ii
      xtry = xsega(n,i)
      ztry = zsega(n,i)
      itry = (xtry-ccon)/deex + 2.0_R8
      jtry = (ztry-zzero)/deez + 2.0_R8
      ptry = pinterp(xtry,ztry,itry,jtry)
      aisegt = i
      if(psir.ge.psisep .and. ptry.lt.psisep) go to 219
      psir = ptry
      xr = xtry
      zr = ztry
      aisegr = aisegt
      go to 220
  219 ilook = 1
      xave = .5_R8*(xr+xtry)
      zave = .5_R8*(zr+ztry)
      disave = sqrt((xave-xsep(1))**2+(zave-zsep(1))**2)
      if(    (zave*zsep(1).lt.0  .or. xsep(1).le.0) .and. xsep(2).gt.0)  &  
     &    disave=sqrt((xave-xsep(2))**2+(zave-zsep(2))**2)
      if(disave.gt.dissav) go to 220
      aisegtt = aisegt
      aisegrt = aisegr
      dissav = disave
  220 continue
      go to(200,149),ilook+1
  149 continue
!
      aigt = aisegrt
      ailt = aisegtt
      go to 160
  150 aigt = aisegtt
      ailt = aisegrt
  160 continue
!
!.....binary search to find strike point
      do 180 l=1,12
      aisegt = .5_R8*(aigt+ailt)
      call divdis(n,aisegt,xtry,ztry,dtry)
      itry = (xtry-ccon)/deex + 2.0_R8
      jtry = (ztry -zzero)/deez + 2.0_R8
      ptry = pinterp(xtry,ztry,itry,jtry)
      if(ptry .gt. psisep) go to 185
      ailt = aisegt
      go to 180
  185 continue
      aigt = aisegt
  180 continue
      xinter(istrike,n) = xtry
      zinter(istrike,n) = ztry
      aisegt = .5_R8*(aigt+ailt)
      call divdis(n,aisegt,xc,zc,disc)
      strike(n,istrike) = disc
      ais(istrike) = aisegt
!
!.....calculate distance from strike point to x-point dsep
      if(zc*zsep(1) .lt. 0) go to 181
      dsep(n,istrike) = sqrt((xc-xsep(1))**2 + (zc-zsep(1))**2)
      go to 100
  181 continue
      if(zc*zsep(2) .lt. 0) go to 182
      dsep(n,istrike) = sqrt((xc-xsep(2))**2 + (zc-zsep(2))**2)
      go to 100
  182 dsep(n,istrike) = 0._R8
      go to 100
  199 continue
!
!
!.....separatrix surface does not intersect plate
  200 continue
      strike(n,istrike) = 0._R8
      ais(istrike) = 0._R8
      dsep(n,istrike) = 0._R8
  100 continue
!
!.....power per meter**2 in each zone segment
!      note:  one sided exponentials with no power deposited between
!              strike points
      sum1 = 0._R8
      sum2 = 0._R8
      do 95 i=1,nseg(n)
      xseg1 = xsega(n,i)
      zseg1 = zsega(n,i)
      iseg1 = (xseg1-ccon)/deex + 2.0_R8
      jseg1 = (zseg1-zzero)/deez + 2.0_R8
      pseg1 = pinterp(xseg1,zseg1,iseg1,jseg1)
!
      xseg2 = xsega(n,i+1)
      zseg2 = zsega(n,i+1)
      iseg2 = (xseg2-ccon)/deex + 2.0_R8
      jseg2 = (zseg2-zzero)/deez + 2.0_R8
      pseg2 = pinterp(xseg2,zseg2,iseg2,jseg2)
!
      psega(n,i) = .5_R8*(pseg1+pseg2)
      dpsega(n,i) = abs(pseg1-pseg2)
      pseg = psega(n,i)
      if(iplim.gt.0) go to 97
!
!.....divertor form factor
      formf(i) = exp(-abs((pseg-psisep)/pscrape))*dpsega(n,i)
      aih = i+0.5_R8
      if(ais(1).gt.0 .and. aih.le.ais(1)) sum1 = sum1 + formf(i)
      if(ais(2).gt.0 .and. aih.ge.ais(2)) sum2 = sum2 + formf(i)
      go to 95
   97 continue
!
!.....limiter form factor
      formf(i) = exp(-abs((pseg-psilim)/pscrape))*dpsega(n,i)
!
!.....check for wrong side of x-point
      if(iexs(iseg1,jseg1).eq.1.or.iexs(iseg2,jseg2).eq.1) formf(i)      &  
     & =0._R8
      sum1 = sum1 + formf(i)
   95 continue
      fpln1=0.5_R8
      if(iplim.lt.0) fpln1=fplate(n,1)
      fpln2=0.0_R8
! note fplate(n,1) and fplate(n,2) have been hardwired in as 0.5, 0.0
! for limiter limited configuration.
      if(iplim.lt.0) fpln2=fplate(n,2)
      do 96 i=1,nseg(n)
      area = asegah(n,i)
      if(iplim.gt.0) go to 98
      aih = i+0.5_R8
!
!.....mwatts per m**2 on divertor plate
!
!.....calculate heat flux based on one sided exponentials
      term1 = 0._R8
      if(sum1.gt.0 .and. ais(1).gt.0 .and. aih.le.ais(1))                &  
     &term1 = hfluxav*(formf(i)/(sum1*area))*fpln1*1.E-6_R8
!
      term2 = 0._R8
      if(sum2.gt.0 .and. ais(2).gt.0 .and. aih.ge.ais(2))                &  
     &term2 = hfluxav*(formf(i)/(sum2*area))*fpln2*1.E-6_R8
!
      hplate(n,i) = term1 + term2
      go to 96
   98 continue
      hplate(n,i) = 0._R8
      if(sum1.gt.0._R8)                                                  &  
     &hplate(n,i) = hfluxav*(formf(i)/(sum1*area))*fpln1*1.E-6_R8
   96 continue
!
      isw=0
      do 91 i=1,nseg(n)
      call grap(3,zsegah(n,i),xsegah(n,i),gradsq,dpsdxh(n,i),            &  
     &   dpsdzh(n,i),gsval,psval,psixz,psixx,psizz,isw)
91    continue
  101 continue
!
!
      plossmw=(hfluxav-prad)*1.E-6_R8
      irecdv = irecdv + 1
      write(ndivlu,6001) times,psilim,plimiter,psisep,pgrad,pgradin,     &  
     &                   pscrape,plossmw,gzero,                          &  
     & ((xinter(istrike,n),zinter(istrike,n),istrike=1,2),n=1,nplate),   &  
     &                   (dsep(n,1),dsep(n,2),n=1,nplate),               &  
     &                   (strike(n,1),strike(n,2),n=1,nplate)
      write(ndivlu,6001) ((psega(n,i),i=1,nseg(n)),n=1,nplate)
      write(ndivlu,6001) ((dpsega(n,i),i=1,nseg(n)),n=1,nplate)
      write(ndivlu,6001) ((hplate(n,i),i=1,nseg(n)),n=1,nplate)
      write(ndivlu,6001) ((dpsdxh(n,i),i=1,nseg(n)),n=1,nplate)
      write(ndivlu,6001) ((dpsdzh(n,i),i=1,nseg(n)),n=1,nplate)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
