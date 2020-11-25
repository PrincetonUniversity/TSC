      subroutine lhwrite(nlhcdp,iopen)
!
      USE CLINAM
      USE RUNAWAY
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iopen,nlhcdp,imnpl,imxpl,jmnpl,jmxpl,i,j,nxpl,nzpl
      INTEGER npsitp,l,n,ios11,kcycmod
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 amassh,sumhe,anum,ps,pval,ppval,gval
      REAL*8 gpval,gppval,gval2j,gpj,gppj,gval2jm,gpjm,gppjm
      REAL*8 factori,rgzero,bgzero
!============
!     dimension anecc(ppsi),tekev(ppsi),tikev(pimp+2,ppsi),
!    1          anicc(pimp+2,ppsi),
!    1          amass(pimp+2),achrg(pimp+2),
!    2          rho(ppsi),pary(ppsi),ppary(ppsi),gary(ppsi),
!    3          gpary(ppsi),vptemp(ppsi)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: anecc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tekev
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tikev
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: anicc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: amass
      REAL*8, ALLOCATABLE, DIMENSION(:) :: achrg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ppary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gpary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: vptemp
!============      
      IF(.not.ALLOCATED(anecc)) ALLOCATE( anecc(ppsi), STAT=istat)
      IF(.not.ALLOCATED(tekev)) ALLOCATE( tekev(ppsi), STAT=istat)
      IF(.not.ALLOCATED(tikev)) ALLOCATE( tikev(pimp+2,ppsi),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(anicc)) ALLOCATE( anicc(pimp+2,ppsi),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(amass)) ALLOCATE( amass(pimp+2), STAT=istat)
      IF(.not.ALLOCATED(achrg)) ALLOCATE( achrg(pimp+2), STAT=istat)
      IF(.not.ALLOCATED(rho)) ALLOCATE( rho(ppsi), STAT=istat)
      IF(.not.ALLOCATED(pary)) ALLOCATE( pary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ppary)) ALLOCATE( ppary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(gary)) ALLOCATE( gary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(gpary)) ALLOCATE( gpary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(vptemp)) ALLOCATE( vptemp(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : lhwrite  ' 
!============      
!
!.....first define special quantities needed to write disk file
      call sprop
      imnpl = imag
      imxpl = imag
      jmnpl = jmag
      jmxpl = jmag
!
!
      do 710 i=iminn,imaxx
      do 710 j=jminn,jmaxx
      if(iexv(i,j).eq.1) go to 710
      if(iexs(i,j).eq.1) go to 710
      if(psi(i,j) .gt. psilim) go to 710
      imnpl = min(i,imnpl)
      imxpl = max(i,imxpl)
      jmnpl = min(j,jmnpl)
      jmxpl = max(j,jmxpl)
  710 continue
      imnpl = imnpl-1
      imxpl = imxpl+1
      jmxpl = jmxpl+1
      if(isym.eq.0) jmnpl = jmnpl-1
      nxpl = imxpl-imnpl+1
      nzpl = jmxpl-jmnpl+1
!
      nspc = 1
      amassh = 1.6726E-24_R8
      amass(1) = amassh*amgas
      achrg(1) = 1._R8
      npsitp = npsit+1
      do 10 l=2,npsitp
   10 anicc(1,l) = anhy(l)*1.E-6_R8
!
!.....check if helium is a species
      sumhe = 0._R8
      do 20 l=2,npsitp
   20 sumhe = sumhe + anhe(l)
      if(sumhe.le.0) go to 31
      nspc = nspc + 1
      amass(nspc) = 4._R8*amassh
      achrg(nspc) = 2._R8
      do 30 l=2,npsitp
   30 anicc(nspc,l) = anhe(l)*1.E-6_R8
   31 continue
      if(iimp.le.0) go to 61
!
!.....check if oxygen is a species
      if(fraci(1).le.0) go to 41
      nspc = nspc + 1
      amass(nspc) = 16._R8*amassh
      achrg(nspc) = 8._R8
      do 40 l=2,npsitp
   40 anicc(nspc,l) = anox(l)*1.E-6_R8
   41 continue
!
!.....check if carbon is a species
      if(fraci(2).le.0) go to 51
      nspc = nspc + 1
      amass(nspc) = 12._R8*amassh
      achrg(nspc) = 6._R8
      do 50 l=2,npsitp
   50 anicc(nspc,l) = anca(l)*1.E-6_R8
   51 continue
!
!.....check if iron is a species
      if(fraci(3).le.0) go to 61
      nspc = nspc + 1
      amass(nspc) = 52._R8*amassh
      achrg(nspc) = 26._R8
      do 60 l=2,npsitp
   60 anicc(nspc,l) = anca(l)*1.E-6_R8
   61 continue
      do 100 l=2,npsitp
      anecc(l) = ane(l)*1.E-6_R8
      tekev(l) =  te(l)*1.E-3_R8
      do 90 n=1,nspc
      tikev(n,l) =  ti(l)*1.E-3_R8
   90 continue
!     voltlp(l) = .5*(as(l)+as(l-1))*udsv
      anum = l-1.5_R8
      rho(l) = sqrt(anum/(npsit-1+fraclst))
      ps = xsv(l)
      call peval(ps,2,pval,ppval,imag,jmag)
      call geval(ps,2,gval,gpval,gppval,imag,jmag)
      pary(l) = pval*udsp
      ppary(l)= ppval*udsp
      gary(l) = gval
      gpary(l)= gpval
      vptemp(l)=.5_R8*tpi*(qprof2(l)*vp2(l) + qprof2(l+1)*vp2(l+1))
!
!.....new calculation of loop voltage...cell centered
      call geval(xsv2(l  ),2,gval2j ,gpj ,gppj ,imag,jmag)
      call geval(xsv2(l-1),2,gval2jm,gpjm,gppjm,imag,jmag)
!
      factori= (vp2(j)*bsqar(j))/ (1.E8_R8*tpi)
      voltlp(l) = 0.5_R8*(etpara(l)*etafac(l)+etpara(l-1)*etafac(l-1))*  &  
     & udsv*(gval**2*rdpsi*(gxmja2(l)/gval2j - gxmja2(l-1)/gval2jm)      &  
     &   - (tpi*factori)*(ajavcd(l)+ajavfw(l)+ajavlh(l)+ajavbs(l)+ajavec(l)))
  100 continue
      rgzero = xplas
      bgzero = gzero/xplas
      if(iopen.eq.1) go to 103
      if(ifk.eq.2) go to 101
!
!
!
!     Create disk file with name LHCDOUa if running on channel a
!                                LHCDOUb "     "     "    "    b  etc.
!            assign to logical unit nlhcd (set in main program)
!
!     ilhcdou(1:6) = 'lhcdou'
!     ilhcdou(7:7) = isuffix(1:1)
      if( numargs .lt. 1 ) then
         ilhcdou = 'lhcdou' // isuffix(1:1)
      else
         ilhcdou = 'lhcdou' // '.' // trim(suffix)
      end if
      open(nlhcdp,file=trim(ilhcdou),status='unknown',iostat=ios11)
  101 continue
      rewind nlhcdp
  103 continue
!
!
!      Description of variables:
!            (note--everything in MKS if not specified)
!
!      name:   80 character description of run from TSC title card
!      npsit: number of flux surfaces for surface averaged quantities
!              (note npsit = npsitp-1)
!      nspc:   total number of ion species
!      kcycle: TSC time cycle number (starts at 0)
!      times:  TSC problem time in sec
!      anecc:  electron density in particles/cc
!      anicc:       ion density in particles/cc
!      tekev:  electron temp in kev
!      tikev:  ion temp in kev
!      amass:  mass of each ion species in grams
!      achrg:  charge of eqch ion species (hydrogen is 1)
!      voltlp: toroidal loop voltage in volts
!      rho:    sqrt of the normalized toroidal flux
!      xsv:    poloidal flux per radian
!      pary:   pressure (MKS units)
!      ppary:  derivative of pressure wrt xsv
!      gary:   toroidal field function R*Bt (MKS)
!      gpary:  derivative of gary wrt xsv
!      vptemp: derivative of volume wrt xsv
!      nxpl:   number of cartesian mesh points in x direction (for psi)
!      nzpl:   number of cartesian mesh points in z direction (for psi)
!      isym:   symmetry option,  0-no symmetry    1- up/down symmetry
!      iplim:  limiter switch   pos-plasma rests on limiter   neg-diverted
!      xary(imnpl):leftmost boundary of cartesian mesh
!      xary(imxpl): rightmost boundary of cartesian mesh
!      zary(jmnpl):   bottom boundary of cartesian mesh
!      zary(jmxpl):    top boundary of cartesian mesh
!      psimin,psilim:  pol flux per radian at mag axis and P/V boundary
!      xmag,zmag:      (x,z) coordinates of magnetic axis
!      rgzero:    nominal major radius of machine
!      bgzero:    vacuum field strength at rgzero
!      apl:       plasma current in amperes
!      psep,xsep,zsep:  flux value,x and z coordinates of separatrix (2)
!      psi:       poloidal flux per radian
!
      write(nlhcdp,1001) (name(i),i=1,8)
      kcycmod = mod(kcycle,100000)
      write(nlhcdp,1002) npsit,nspc,kcycmod,times
      write(nlhcdp,1003) (anecc(l),l=2,npsitp)
      write(nlhcdp,1003) ( tekev(l),l=2,npsitp)
      write(nlhcdp,1003) ((anicc(n,l),l=2,npsitp),n=1,nspc)
      write(nlhcdp,1003) ((tikev(n,l),l=2,npsitp),n=1,nspc)
      write(nlhcdp,1003) (amass(n),n=1,nspc)
      write(nlhcdp,1003) (achrg(n),n=1,nspc)
      write(nlhcdp,1003) (voltlp(l),l=2,npsitp)
      write(nlhcdp,1003) (rho(l),l=2,npsitp)
      write(nlhcdp,1003) (vptemp(l),l=2,npsitp)
!
      write(nlhcdp,1003) (xsv(l),l=2,npsitp)
      write(nlhcdp,1003) (pary(l),l=2,npsitp)
      write(nlhcdp,1003) (ppary(l),l=2,npsitp)
      write(nlhcdp,1003) (gary(l),l=2,npsitp)
      write(nlhcdp,1003) (gpary(l),l=2,npsitp)
      write(nlhcdp,1004) nxpl,nzpl,isym,iplim
      write(nlhcdp,1003) xary(imnpl),xary(imxpl),zary(jmnpl),zary(jmxpl)    
      write(nlhcdp,1003) psimin,psilim,xmag,zmag
      write(nlhcdp,1003) rgzero,bgzero,apl
      write(nlhcdp,1003) psep(1),xsep(1),zsep(1)
      write(nlhcdp,1003) psep(2),xsep(2),zsep(2)
      write(nlhcdp,1003) ((psi(i,j),i=imnpl,imxpl),j=jmnpl,jmxpl)
      if(ifk.eq.2.or.iopen.eq.1) go to 102
      close(nlhcdp)
  102 continue
!
      return
 1001 format(8a10)
 1002 format(3i5,1pe16.6)
 1003 format(1p5e16.6)
 1004 format(4i5)
!.....DESCRIPTION OF PARAMETERS PASSED TO LSC:
!
!     bpowsmw   -  total LH power in MW
!     NLHCD     -  Logical Unit number for writing to LSC
!     NLHCD2    -  Logical Unit number for reading from LSC
!     NLSCGNU   -  LU number LSC to write GNU plot information
!     NOUT      -  LU for standard output
!     NLHCDIN   -  LU for LSC input
!     NTERM     -  LU for writing to terminal
!     IRAYT     -  =1 for LSC to trace rays,  2 otherwise
!     IPLOT     -  =1 to make plots,  0 otherwise
!     IERR      -  =0 on normal exit
!
!     ARRAYS WRITTEN BY LSC ON DISK FILE:
!
!     powtsc   -  power density in W/cc
!     currtsc  -  current density in amps/cm**2
!     djdetsc  -  derivative of current density wrt E (cell centered)
!     djde2tsc -  derivative of current density wrt E (boundary centered)
!
!     call lsc(bpowsmw,nlhcd,nlhcd2,nlscgnu,nout,nlhcdin,nterm,irayt
!    1        ,iplot,ierr)
!     iraytd = irayt
!     if(ierr.gt.0) ineg=45
!     if(ierr.lt.0) write(nout,1991) ierr
!1991 format("LSC returned ierr=",i5)
!     if(ierr.lt.0) write(nterm,1991) ierr
!     if(ineg.ne.0) return
!     rewind nlhcd2
!     read(nlhcd2,1001) (powtsc(l),l=2,npsit)
!     read(nlhcd2,1001) (currtsc(l),l=2,npsit)
!     read(nlhcd2,1001) (djdetsc(l),l=2,npsit)
!     read(nlhcd2,1001) (djdets2(l),l=2,npsit)
!1001 format(1p5e16.6)
! 651 continue
!
!     powsum = 0.
!     cursum = 0.
!     do 653 l=2,npsi
!     powsum = powsum + powtsc(l)*vp(l)*dpsi
!     cursum = cursum + currtsc(l)*1.e4*vp(l)*dpsi/(tpi*xplas)
! 653 continue
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
