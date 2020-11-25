      subroutine spdd3d
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE CLINAM
      USE SAPROP
      USE SCR3
      USE FEED
      USE DIAGCOM
      USE LOOP1
      USE PROBE1
      USE ARRAYS

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER kchop,ksign,nchop,icount,ios7,ios22,ishot,ktime
      INTEGER limitr,kfcoil,ksilop,kmgprbe,kecoil,kvesel,kesum
      INTEGER kslref,kco2r,kco2v,kgrean,i,k,i1,j1,i2,j2,ii,jw
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 vps
      REAL*8 vpull,d,d2,v,t2,hv1,hv2,delt,rzer,zzer,ang,r1,z1,r2,z2
      REAL*8 p1,pinterp,p2,pone,tay,tt,zp00,zp2,dfzp0,dfzp2,zp,dfzp
      REAL*8 tvde,zpp,rpp,wvspaip,zxp,shapep,gapinp,dc,v01,csum
      REAL*8 plcur,tekev,rlin
!============
!     common/diagcom/verta,vertb,dfzpd
!       COMMON /loop1/
!    *PSF2A,PSF4A,PSF5A,PSF6NA,
!    *PSF8A,PSF1B,PSF2B,PSF3B,PSF4B,PSF5B,PSF6NB,
!    *PSF8B,PSF9B,PSI7A,PSF6FA,PSI7B,PSF6FB
!     COMMON /probe1/
!    *AMPI9A ,AMPI67A ,AMPI6FA ,AMPI66M ,AMPI6FB ,
!    *AMPI67B ,AMPI7FB ,AMPI9B ,AMPI89B ,
!    *AMPI8B ,AMPI4B ,AMPI3B ,AMPI2B ,AMPI1B ,
!    *AMPI11M ,AMPI1A ,AMPI2A ,AMPI5A
!
!     common / arrays / PF(18),PFE(18),FCOM(18),vchopper(18),FCOMSV(18)
!     common/feed/ viabturn(18),viabturno(18)
!
!      probe identification
!
! no.      name       no      name
!
! 1.   mpi11m067      16.   mpi1b067
! 2.   mpi1a067       17.   mpi2b067
! 3.   mpi2a067(157)  18.   mpi3b067
! 4.   mpi3a067       19.   mpi4b067
! 5.   mpi4a067       20.   mpi5b067
! 6.   mpi5a067       21.   mpi8b067
! 7.   mpi8a067       22.   mpi89b067
! 8.   mpi9a067       23.   mpi9b067
! 9.   mpi79a067      24.   mpi79b067
! 10.  mpi7fa067      25.   mpi7fb067
! 11.  mpi7na067      26.   mpi7nb067
! 12.  mpi67a067      27.   mpi67b067
! 13.  mpi6fa067      28.   mpi6fb067
! 14.  mpi6na067      29.   mpi6nb067
! 15.  mpi66m067
!
!
!     psi loop identification
!
! no.   name     no.   name
!
! 1.   psf1a     22.  psi34a
! 2.   PSF2A     23.  psii45a
! 3.   psf3a     24.  psi58a
! 4.   PSF4A     25.  psi9a
! 5.   PSF5A     26.  psf7fa
! 6.   PSF6NA    27.  PSI7A
! 7.   psf7na    28.  PSF6FA
! 8.   PSF8A     29.  psi6a
! 9.   psf9a     30.  psi12b
! 10.  PSF1B     31.  psi23b
! 11.  PSF2B     32.  psi34b
! 12.  PSF3B     33.  psi45b
! 13.  PSF4B     34.  psi58b
! 14.  PSF5B     35.  psi9b
! 15.  PSF6NB    36.  psf7fb
! 16.  psf7nb    37.  PSI7B
! 17.  PSF8B     38.  PSF6FB
! 18.  PSF9B     39.  psi6b
! 19.  psi11m    40.  psi89fb
! 20.  psi12a    41.  psi89nb
! 21.  psi23a
!
!
      integer plimitr,pkco2r,pkco2v,pkfcoil,pksilop,pmgprbe,pkecoil,     &  
     &        pkvesel,pktime,pkesum,pkgrean,pkslref
      parameter (plimitr=50   ,pkco2r=1     ,pkco2v=3     ,pkfcoil=18 ,  &  
     &           pksilop=41   ,pmgprbe=64   ,pkecoil=122  ,pkvesel=24 ,  &  
     &           pktime=501   ,pkesum=2     ,pkgrean=10   ,pkslref=1 )
!
!     dimension  sxlim(plimitr)  ,sylim(plimitr)  ,
!    1           schordr(pkco2r) ,schordv(pkco2v)  ,
!    2           srf(pkfcoil)    ,szf(pkfcoil)    ,swf(pkfcoil)    ,
!    3           shf(pkfcoil)    ,saf(pkfcoil)    ,saf2(pkfcoil)  ,
!    4           sturnfc(pkfcoil),srsisfc(pkfcoil),
!    5           srsi(pksilop)   ,szsi(pksilop)   ,
!    6           sxmp2(pmgprbe)  ,symp2(pmgprbe)  ,samp2(pmgprbe)  ,
!    7           sre(pkecoil)    ,sze(pkecoil)    ,swe(pkecoil)    ,
!    8           she(pkecoil)    ,secid(pkecoil)  ,
!    9           srvs(pkvesel)   ,szvs(pkvesel)   ,swvs(pkvesel)   ,
!    .           shvs(pkvesel)   ,savs(pkvesel)   ,savs2(pkvesel)  ,
!    1           srsisvs(pkvesel),
!    2           stime(pktime)   ,sfwtfc(pkfcoil) ,
!    3           sfccurt(pktime,pkfcoil) , seccurt(pktime,pkesum) ,
!    4           spasmat(pktime) ,svloopt(pktime) ,sfwtsi(pksilop) ,
!    5           ssilopt(pktime,pksilop) , sfwtmp2(pmgprbe)        ,
!    6           sexpmpi(pktime,pmgprbe) , sdenrt(pktime,pkco2r)   ,
!    7           sdenvt(pktime,pkco2v)   , srgrean(pktime,pkgrean) ,
!    8           stegrean(pktime,pkgrean), sbcentr(pktime)        ,
!    9           srpp(pktime)            , szpp(pktime)           ,
!    1           selp(pktime)            , shv1volt(pktime)       ,
!    2           shv2volt(pktime)        , pasmap(pktime)         ,
!    3           rppe(pktime)            , v7off(pktime)          ,
!    4          sd1(pktime),sd2(pktime),sv1(pktime),swvsp2p(pktime)
!
!     dimension bfield(pmgprbe),pflux(pksilop),fcoil(pkfcoil),
!    1         bfieldo(pmgprbe),pfluxo(pksilop)
!
!
!
      dimension VPS(18),KCHOP(18),KSIGN(18),NCHOP(18),VPULL(18)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sxlim
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sylim
      REAL*8, ALLOCATABLE, DIMENSION(:) :: schordr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: schordv
      REAL*8, ALLOCATABLE, DIMENSION(:) :: srf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: szf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: swf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: shf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: saf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: saf2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sturnfc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: srsisfc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: srsi
      REAL*8, ALLOCATABLE, DIMENSION(:) :: szsi
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sxmp2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: symp2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: samp2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sre
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sze
      REAL*8, ALLOCATABLE, DIMENSION(:) :: swe
      REAL*8, ALLOCATABLE, DIMENSION(:) :: she
      REAL*8, ALLOCATABLE, DIMENSION(:) :: secid
      REAL*8, ALLOCATABLE, DIMENSION(:) :: srvs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: szvs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: swvs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: shvs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: savs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: savs2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: srsisvs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: stime
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sfwtfc
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sfccurt
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: seccurt
      REAL*8, ALLOCATABLE, DIMENSION(:) :: spasmat
      REAL*8, ALLOCATABLE, DIMENSION(:) :: svloopt
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sfwtsi
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ssilopt
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sfwtmp2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sexpmpi
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sdenrt
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sdenvt
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: srgrean
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: stegrean
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sbcentr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: srpp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: szpp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: selp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: shv1volt
      REAL*8, ALLOCATABLE, DIMENSION(:) :: shv2volt
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pasmap
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rppe
      REAL*8, ALLOCATABLE, DIMENSION(:) :: v7off
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sd1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sd2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sv1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: swvsp2p
      REAL*8, ALLOCATABLE, DIMENSION(:) :: bfield
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pflux
      REAL*8, ALLOCATABLE, DIMENSION(:) :: fcoil
      REAL*8, ALLOCATABLE, DIMENSION(:) :: bfieldo
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pfluxo
!============      
      IF(.not.ALLOCATED(sxlim)) ALLOCATE( sxlim(plimitr), STAT=istat)
      IF(.not.ALLOCATED(sylim)) ALLOCATE( sylim(plimitr), STAT=istat)
      IF(.not.ALLOCATED(schordr)) ALLOCATE( schordr(pkco2r), STAT=istat)
      IF(.not.ALLOCATED(schordv)) ALLOCATE( schordv(pkco2v), STAT=istat)
      IF(.not.ALLOCATED(srf)) ALLOCATE( srf(pkfcoil), STAT=istat)
      IF(.not.ALLOCATED(szf)) ALLOCATE( szf(pkfcoil), STAT=istat)
      IF(.not.ALLOCATED(swf)) ALLOCATE( swf(pkfcoil), STAT=istat)
      IF(.not.ALLOCATED(shf)) ALLOCATE( shf(pkfcoil), STAT=istat)
      IF(.not.ALLOCATED(saf)) ALLOCATE( saf(pkfcoil), STAT=istat)
      IF(.not.ALLOCATED(saf2)) ALLOCATE( saf2(pkfcoil), STAT=istat)
      IF(.not.ALLOCATED(sturnfc)) ALLOCATE( sturnfc(pkfcoil),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(srsisfc)) ALLOCATE( srsisfc(pkfcoil),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(srsi)) ALLOCATE( srsi(pksilop), STAT=istat)
      IF(.not.ALLOCATED(szsi)) ALLOCATE( szsi(pksilop), STAT=istat)
      IF(.not.ALLOCATED(sxmp2)) ALLOCATE( sxmp2(pmgprbe), STAT=istat)
      IF(.not.ALLOCATED(symp2)) ALLOCATE( symp2(pmgprbe), STAT=istat)
      IF(.not.ALLOCATED(samp2)) ALLOCATE( samp2(pmgprbe), STAT=istat)
      IF(.not.ALLOCATED(sre)) ALLOCATE( sre(pkecoil), STAT=istat)
      IF(.not.ALLOCATED(sze)) ALLOCATE( sze(pkecoil), STAT=istat)
      IF(.not.ALLOCATED(swe)) ALLOCATE( swe(pkecoil), STAT=istat)
      IF(.not.ALLOCATED(she)) ALLOCATE( she(pkecoil), STAT=istat)
      IF(.not.ALLOCATED(secid)) ALLOCATE( secid(pkecoil), STAT=istat)
      IF(.not.ALLOCATED(srvs)) ALLOCATE( srvs(pkvesel), STAT=istat)
      IF(.not.ALLOCATED(szvs)) ALLOCATE( szvs(pkvesel), STAT=istat)
      IF(.not.ALLOCATED(swvs)) ALLOCATE( swvs(pkvesel), STAT=istat)
      IF(.not.ALLOCATED(shvs)) ALLOCATE( shvs(pkvesel), STAT=istat)
      IF(.not.ALLOCATED(savs)) ALLOCATE( savs(pkvesel), STAT=istat)
      IF(.not.ALLOCATED(savs2)) ALLOCATE( savs2(pkvesel), STAT=istat)
      IF(.not.ALLOCATED(srsisvs)) ALLOCATE( srsisvs(pkvesel),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(stime)) ALLOCATE( stime(pktime), STAT=istat)
      IF(.not.ALLOCATED(sfwtfc)) ALLOCATE( sfwtfc(pkfcoil), STAT=istat)
      IF(.not.ALLOCATED(sfccurt)) ALLOCATE( sfccurt(pktime,pkfcoil),     &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(seccurt)) ALLOCATE( seccurt(pktime,pkesum),      &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(spasmat)) ALLOCATE( spasmat(pktime), STAT=istat)
      IF(.not.ALLOCATED(svloopt)) ALLOCATE( svloopt(pktime), STAT=istat)
      IF(.not.ALLOCATED(sfwtsi)) ALLOCATE( sfwtsi(pksilop), STAT=istat)
      IF(.not.ALLOCATED(ssilopt)) ALLOCATE( ssilopt(pktime,pksilop),     &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(sfwtmp2)) ALLOCATE( sfwtmp2(pmgprbe),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(sexpmpi)) ALLOCATE( sexpmpi(pktime,pmgprbe),     &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(sdenrt)) ALLOCATE( sdenrt(pktime,pkco2r),        &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(sdenvt)) ALLOCATE( sdenvt(pktime,pkco2v),        &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(srgrean)) ALLOCATE( srgrean(pktime,pkgrean),     &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(stegrean)) ALLOCATE( stegrean(pktime,pkgrean),   &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(sbcentr)) ALLOCATE( sbcentr(pktime), STAT=istat)
      IF(.not.ALLOCATED(srpp)) ALLOCATE( srpp(pktime), STAT=istat)
      IF(.not.ALLOCATED(szpp)) ALLOCATE( szpp(pktime), STAT=istat)
      IF(.not.ALLOCATED(selp)) ALLOCATE( selp(pktime), STAT=istat)
      IF(.not.ALLOCATED(shv1volt)) ALLOCATE( shv1volt(pktime),           &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(shv2volt)) ALLOCATE( shv2volt(pktime),           &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(pasmap)) ALLOCATE( pasmap(pktime), STAT=istat)
      IF(.not.ALLOCATED(rppe)) ALLOCATE( rppe(pktime), STAT=istat)
      IF(.not.ALLOCATED(v7off)) ALLOCATE( v7off(pktime), STAT=istat)
      IF(.not.ALLOCATED(sd1)) ALLOCATE( sd1(pktime), STAT=istat)
      IF(.not.ALLOCATED(sd2)) ALLOCATE( sd2(pktime), STAT=istat)
      IF(.not.ALLOCATED(sv1)) ALLOCATE( sv1(pktime), STAT=istat)
      IF(.not.ALLOCATED(swvsp2p)) ALLOCATE( swvsp2p(pktime), STAT=istat)
      IF(.not.ALLOCATED(bfield)) ALLOCATE( bfield(pmgprbe), STAT=istat)
      IF(.not.ALLOCATED(pflux)) ALLOCATE( pflux(pksilop), STAT=istat)
      IF(.not.ALLOCATED(fcoil)) ALLOCATE( fcoil(pkfcoil), STAT=istat)
      IF(.not.ALLOCATED(bfieldo)) ALLOCATE( bfieldo(pmgprbe),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(pfluxo)) ALLOCATE( pfluxo(pksilop), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : spdd3d  ' 
!============      
 
!---------------
      DATA VPS /                                                         &  
     &  225._R8, 436._R8,   0._R8, 230._R8, 230._R8, 518._R8, 585._R8,   &  
     & 230._R8,   0._R8,                                                 &  
     &  259._R8, 307._R8, 436._R8, 225._R8, 288._R8, 518._R8, 585._R8,   &  
     & 288._R8, 230._R8/
!
        DATA KCHOP /                                                     &  
     &  1, 1, 0, 1, 1, 2, 2, 2, 0,                                       &  
     &  1, 1, 1, 1, 1, 2, 2, 1, 1 /
 
        DATA KSIGN /                                                     &  
     & -1,-1,-1,+1,+1,-1,-1,+1,+0,                                       &  
     & -1,-1,-1,+1,+1,-1,-1,+1,+1 /
!
        DATA NCHOP /                                                     &  
     &  2, 2, 1, 1, 1, 3, 3, 1, 0,                                       &  
     &  2, 1, 1, 1, 2, 4, 4, 2, 2 /
!
!------------------
!
!.......check if initialization is necessary, otherwise proceed to 100
      if(ifrst(9).ne.1) go to 100
      ifrst(9) = 0
      icount = 0
!
      pfbval = 0._R8
      D=369._R8
      D2=465._R8
      V=360._R8
      T2=398._R8
      HV1=864._R8
      HV2=825._R8
      VPS(1)=d2
      VPS(2)=d2
      VPS(4)=d
      VPS(5)=d
      VPS(6)=hv1
      VPS(7)=hv2
      VPS(8)=v
      VPS(1+9)=d2
      VPS(2+9)=d2
      VPS(3+9)=d2
      VPS(4+9)=v
      VPS(5+9)=v
      VPS(6+9)=hv1
      VPS(7+9)=hv2
      VPS(8+9)=v
      VPS(9+9)=t2
!
!.......................................................................
!
!.....read data file
!
!.......................................................................
!     ifileen(1:4) = 'enin'
!     ifileen(5:5) = isuffix(1:1)
!     ifile7a(1:5) = 'o167a'
!     ifile7a(6:6) = isuffix(1:1)
      if( numargs .lt. 1 ) then
        ifileen = 'enin' // isuffix(1:1)
        ifile7a = 'o167a' // isuffix(1:1)
      else
        ifileen = 'enin' // '.' // trim(suffix)
        ifile7a = 'o167a' // '.' // trim(suffix)
      end if
      open(nenin,file=trim(ifileen),status='old',iostat=ios7)
      open(no167a,file=trim(ifile7a),status='unknown',iostat=ios22)
!
!
!
      read(nenin,7000) ishot,ktime,limitr,kfcoil,ksilop,                 &  
     &                  kmgprbe,kecoil,kvesel,kesum,kslref
      read(nenin,7000) kco2r,kco2v,kgrean
!
      write(nout,1001) ishot,ktime,limitr,kfcoil,ksilop,                 &  
     &              kmgprbe,kecoil,kvesel,kesum,kslref
 1001 format(1h1," ishot= ",i6,"    ktime= ",i6,"    limitr=",i6,        &  
     &   /,1x,"    kfcoil=",i6,"   ksilop= ",i6,"   kmgprbe=",i6,        &  
     &   /,1x,"    kecoil=",i6,"   kvesel= ",i6,"    kesum =",i6,        &  
     &   /,1x,"    kgrean=",i6)
!
      if(ktime .gt. pktime  .or. limitr .gt. plimitr .or.                &  
     &   kfcoil .gt. pkfcoil .or. ksilop .gt. pksilop .or.               &  
     &   kmgprbe .gt. pmgprbe .or. kecoil .gt. pkecoil .or.              &  
     &   kvesel .gt. pkvesel .or. kesum .gt. pkesum .or.                 &  
     &   kslref .gt. pkslref .or. kco2r .gt. pkco2r .or.                 &  
     &   kco2v .gt. pkco2v .or. kgrean .gt. pkgrean) ineg=6
      if(ineg.ne.0) return
!
      read(nenin,7020) (sxlim(i),sylim(i),i=1,limitr)
      read(nenin,7020) (schordr(i),i=1,kco2r)
      read(nenin,7020) (schordv(i),i=1,kco2v)
      read(nenin,7020) (srf(k),szf(k),swf(k),shf(k),saf(k),              &  
     &  saf2(k),sturnfc(k),srsisfc(k),k=1,kfcoil)
      read(nenin,7020) (srsi(i),szsi(i),i=1,ksilop)
      read(nenin,7020) (sxmp2(i),symp2(i),samp2(i),i=1,kmgprbe)
      write(nout,7021) (i,srsi(i),szsi(i),i=1,ksilop)
      write(nout,7022) (i,sxmp2(i),symp2(i),samp2(i),i=1,kmgprbe)
 7021 format(i6,1p2e12.4)
 7022 format(i6,1p3e12.4)
 7023 format(" i      rloop       zloop")
 7024 format(" i      sxmp2       symp2       samp2")
      close(nenin)
 7000 format(1x,11i7)
 7020 format(1x,1p6e13.6)
!
!
!
!.....end of initialization section
!
!
  100 continue
!
      if(kcycle.lt.0) return
!
!
!
!.....field from simulation at probe points
      delt = .01_R8
!     write(nout,8800)
!8800 format(////////,"   i      sxmp2(i)    symp2(i)    samp2(i)")
!8801 format(1x,i5,1p3e12.4)
      do 200 i=1,kmgprbe
      rzer = sxmp2(i)
      zzer = symp2(i)
      ang = samp2(i)*tpi/360._R8
!     write(nout,8801)  i, rzer, zzer, samp2(i)
      r1 = rzer - delt*sin(ang)
      z1 = zzer + delt*cos(ang)
      r2 = rzer + delt*sin(ang)
      z2 = zzer - delt*cos(ang)
      i1 = (r1-ccon)/deex + 2
      j1 = (z1-zzero)/deez + 2
      i2 = (r2-ccon)/deex + 2
      j2 = (z2-zzero)/deez + 2
      p1 = pinterp(r1,z1,i1,j1)
      p2 = pinterp(r2,z2,i2,j2)
      bfieldo(i) = bfield(i)
      bfield(i) = (p1-p2)/(2._R8*rzer*delt)
      if(kcycle.le.1) go to 200
      bfield(i) = 0.5_R8*(bfield(i)+bfieldo(i))
  200 continue
!     stop
!
!.....flux from simulation
      pone = 0._R8
      do 300 i=1,ksilop
      r1 = srsi(i)
      z1 = szsi(i)
      i1 = (r1-ccon)/deex + 2
      j1 = (z1-zzero)/deez + 2
      p1 = pinterp(r1,z1,i1,j1)
      pfluxo(i) = pflux(i)
      pflux(i) =  -(p1+pone)
      pone = pflux(1)
      if(kcycle.le.1) go to 300
      pflux(i) = 0.5_R8*(pflux(i) + pfluxo(i))
  300 continue
      do 305 k=1,18
  305 fcoil(k) = 0._R8
!
!.....actual current in f-coils
      do 307 k=1,18
      do 306 ii=1,nwire
      jw = jwire(ii)
      if(igroupw(ii).ne.k) go to 306
      if(zary(jw).gt.0 .or. k.gt.9)                                      &  
     &  fcoil(k)=fcoil(k)+ccoil(ncoil-nwire+ii)*udsi
      if(zary(jw).lt.0 .and. k.le.9)                                     &  
     & fcoil(k+9)=fcoil(k+9)+ccoil(ncoil-nwire+ii)*udsi
  306 continue
  307 continue
!
!
!......fluxes,fields,and currents for use in feedback systems
!
      PSF1B =  pflux(10)
      PSI7B = pflux(37)
      PSI7A = pflux(27)
      PSF6FB = pflux(38)
      PSF6FA = pflux(28)
      PSF2B  = pflux(11)
      PSF2A  = pflux(2)
      PSF3B  = pflux(12)
      PSF4B  = pflux(13)
      PSF4A  = pflux(4)
      PSF5B  = pflux(14)
      PSF5A  = pflux(5)
      PSF8B  = pflux(17)
      PSF8A  = pflux(8)
      PSF9B = pflux(18)
      PSF6NA = pflux(6)
      PSF6NB = pflux(15)
!
      AMPI2A  = bfield(3)
      AMPI2B  = bfield(17)
      AMPI67A = bfield(12)
      AMPI67B = bfield(27)
      AMPI66M = bfield(15)
      AMPI6FA = bfield(13)
      AMPI6FB = bfield(28)
      AMPI1B  = bfield(16)
      AMPI9A  = bfield(8)
      AMPI89B = bfield(22)
      AMPI5A  = bfield(6)
      AMPI8B = bfield(21)
      AMPI9B = bfield(23)
      AMPI7FB = bfield(25)
      AMPI4B = bfield(19)
      AMPI3B = bfield(18)
      AMPI11M = bfield(1)
      AMPI1A = bfield(2)
!
!....calculate current per turn...
      do 171 i=1,18
      PF(i) = (fcoil(i)/sturnfc(i))*1.E-3_R8
  171 continue
      TAY = 4.0_R8*(dt+dtold)*udst*1000._R8
!     TAY = 1.
      if(TAY .le. 0) TAY = 1._R8
!
!
!
!
      AIP = .5_R8*(apl+aplo)
      TT = times
      ZP00 = zp2
      DFZP0 = dfzp2
      zp2 = zp
      dfzp2 = dfzp
!
!.....Define TVDE,ZPP,RPP,WVSPAIP,ZXP,SHAPEP,GAPINP,ZP00,DFZP0 for shot 83964 time 1.800
!
      TVDE = 2.000_R8
      ZPP = 1.0_R8
      RPP = 0.56_R8
      WVSPAIP = 0.99792_R8
      ZXP = -5.9509_R8
      SHAPEP = -1.9501_R8
      GAPINP = 1.8768_R8
!
!
!.....special for shot 88806
      ZPP = -.0081912_R8
      RPP = .56415_R8
      WVSPAIP = .99792_R8
      ZXP = -6.2012_R8
      SHAPEP = -2.4887_R8
      GAPINP =  1.8724_R8
!
      call control(AIP,TT,TVDE,ZPP,RPP,WVSPAIP,ZXP,SHAPEP,GAPINP,        &  
     &                  ZP00,zp,DFZP0,dfzp,TAY)
!
!
      do 278 i=1,18
      if(FCOM(i) .gt. 8._R8) FCOM(i) = 8._R8
      if(FCOM(i) .lt.-8._R8) FCOM(i) =-8._R8
      DC = 0.0375_R8*FCOM(i)+0.5_R8
      if(KCHOP(i).eq.0) then
      vchopper(i) = 0._R8
      go to 278
      endif
!
      call pulln (VPULL(i),KCHOP(i),PF(i),NCHOP(i),FCOM(i),              &  
     &            DC,VPS(i),V01)
      vchopper(i) = ksign(i)*(DC*V01 + (1._R8-DC)*VPULL(i))
      viabturn(i) = vchopper(i)/sturnfc(i)*usdv/tpi
  278 continue
!
!...diag
!     write(nout,1278) kcycle,viabturn(15),vertb,dfzp,dfzpd
!1278 format(" kcycle,viabturn(15),vertb,dfzp,dfzpd",i5,1p4e12.4)
!
!
!     do 60 l=1,ntpts
!  60 fact(l) = 0.
!     fact(istart) = 1.
!
!.....only write every nskip2 cycles
      if(kcycle.le.0) return
      if(iplt2.le.nskip2) return
      icount = 0
!
!.....calculate total current inside vessel and plasma
      csum = 0._R8
      do 400 ii=1,nwire
      if(igroupw(ii).le.18 .or. igroupw(ii).ge.43) go to 400
      csum = csum + ccoil(ncoil-nwire+ii)*udsi
  400 continue
      plcur = tcurdtp*tpi*udsi + csum
!
      write(no167a,1009) ishot,kcycle,ksilop,kmgprbe,times,csum,plcur
       tekev = te(1)*1.E-3_R8
      rlin = 0._R8
      write(no167a,1010) (fcoil(k),k=1,18),tekev,rlin,                   &  
     &                   (bfield(i),i=1,kmgprbe),                        &  
     &                   (pflux(i),i=1,ksilop),                          &  
     &                   (FCOM(i),i=1,18),                               &  
     &                   (VPULL(i),i=1,18),                              &  
     &                   (PF(i),i=1,18)
 1009 format(4i10,1p3e20.12)
 1010 format(1x,1p8e13.5)
      return
  101 continue
!
!.....error exit
      write(nout,1777)
 1777 format(" error , ier ne 0 after call to chop in spdata")
      ineg=15
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
