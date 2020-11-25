      subroutine wrpest
!......6.75 wrpest
!
!.....writes a pest formated output file
!
      USE CLINAM
      USE SAPROP
      USE SCR10
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER, PARAMETER :: nwallnstx=32
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nxpest,nzpest,n,i,jpest,j,jj,k,idum,nzdim,nwalliter
      INTEGER nlimiter,kk
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rwalliter,zwalliter,rlimiter,zlimiter,alrpest
      REAL*8 alzpest,rzpest,rpest,pmpest,plpest,dps,ps,pval,ppval
      REAL*8 gval,gpval,gppval,t1,t3,t5,rdim,zdim,zmid,bcentr,xdum
      REAL*8 zip,psival,pprimesav2,pprimesav1,ffprimsav2,ffprimsav1
      REAL*8 rwallnstx,zwallnstx,ajphi2
!============
!     dimension te2(ppsi),ane2(ppsi),ti2(ppsi),avez2(ppsi),zeff2(ppsi)
!     dimension fpol(pnx),pres(pnx),ffprim(pnx),pprime(pnx),qpsi(pnx)
      dimension rwalliter(41),zwalliter(41),rlimiter(14),zlimiter(14)
      dimension rwallnstx(nwallnstx),zwallnstx(nwallnstx)
!     dimension psiaux(pnx,2*pnz)
!
!...first wall R, m
      data rwalliter /4.047_R8,4.047_R8,4.047_R8,4.047_R8,4.047_R8,      &  
     & 4.047_R8,                                                         &  
     &                4.047_R8,4.314_R8,4.855_R8,5.705_R8,6.494_R8,      &  
     & 7.413_R8,                                                         &  
     &                7.907_R8,8.267_R8,8.478_R8,8.312_R8,7.902_R8,      &  
     & 7.254_R8,                                                         &  
     &                6.348_R8,6.311_R8,6.377_R8,6.263_R8,6.062_R8,      &  
     & 5.877_R8,                                                         &  
     &                5.725_R8,5.618_R8,5.561_R8,5.562_R8,5.298_R8,      &  
     & 5.164_R8,                                                         &  
     &                5.030_R8,4.871_R8,4.703_R8,4.475_R8,4.080_R8,      &  
     & 4.407_R8,                                                         &  
     &                4.481_R8,4.459_R8,4.372_R8,4.231_R8,4.047_R8/
!...first wall Z, m
      data zwalliter /-2.474_R8,-1.523_R8,-0.507_R8,0.510_R8,1.527_R8,   &  
     & 2.544_R8,                                                         &  
     &                3.560_R8,4.300_R8,4.703_R8,4.539_R8,3.943_R8,      &  
     & 3.140_R8,                                                         &  
     &                2.436_R8,1.663_R8,0.667_R8,-0.442_R8,-1.372_R8,-   &  
     & 2.309_R8,                                                         &  
     &                -3.016_R8,-3.095_R8,-3.173_R8,-3.203_R8,-3.227_R8,  &  
     & -3.309_R8,                                                        &  
     &                -3.441_R8,-3.613_R8,-3.902_R8,-4.609_R8,-3.979_R8,  &  
     & -3.765_R8,                                                        &  
     &                -3.664_R8,-3.608_R8,-3.620_R8,-3.675_R8,-3.872_R8,  &  
     & -3.340_R8,                                                        &  
     &                -3.087_R8,-2.903_R8,-2.738_R8,-2.616_R8,-2.474_R8/     
!...first wall R, m
      data rwallnstx / 0.17_R8, 0.17_R8, 0.25_R8, 0.25_R8, 0.28_R8, &
     & 0.28_R8, 0.55_R8, 0.55_R8, 0.59_R8, 0.63_R8, 0.64_R8, 1.12_R8, &
     & 1.12_R8, 1.35_R8, 1.51_R8, 1.70_R8, 1.70_R8, 1.51_R8, 1.35_R8, &
     & 1.12_R8, 1.12_R8, 0.64_R8, 0.63_R8, 0.60_R8, 0.60_R8, 0.55_R8, &
     & 0.55_R8, 0.28_R8, 0.28_R8, 0.25_R8, 0.25_R8, 0.17_R8/
!...first wall Z, m
      data zwallnstx /-1.02_R8, 1.02_R8, 1.14_R8, 1.62_R8, 1.62_R8, &
     & 1.66_R8, 1.66_R8, 1.74_R8, 1.74_R8, 1.72_R8, 1.66_R8, 1.46_R8, &
     & 1.34_R8, 1.04_R8, 0.56_R8, 0.56_R8,-0.56_R8,-0.56_R8,-1.04_R8, &
     & -1.34_R8,-1.46_R8,-1.66_R8,-1.72_R8,-1.74_R8,-1.86_R8,-1.86_R8, &
     & -1.66_R8,-1.66_R8,-1.62_R8,-1.62_R8,-1.14_R8,-1.02_R8/
!...limiter R, m
      data rlimiter /8.271_R8,8.258_R8,8.269_R8,8.275_R8,                &  
     &               8.282_R8,8.283_R8,8.281_R8,8.276_R8,                &  
     &               8.267_R8,8.255_R8,8.240_R8,8.222_R8,                &  
     &               8.201_R8,8.205_R8/
!...limiter Z, m
      data zlimiter /-0.436_R8,-0.274_R8,-0.112_R8,0.000_R8,             &  
     &               0.213_R8,0.375_R8,0.537_R8,0.700_R8,0.862_R8,       &  
     &               1.024_R8,1.185_R8,1.347_R8,1.508_R8,1.670_R8/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: te2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ane2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ti2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: avez2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zeff2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: fpol
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pres
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ffprim
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pprime
      REAL*8, ALLOCATABLE, DIMENSION(:) :: qpsi
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: psiaux
!============      
      IF(.not.ALLOCATED(te2)) ALLOCATE( te2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ane2)) ALLOCATE( ane2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ti2)) ALLOCATE( ti2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(avez2)) ALLOCATE( avez2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(zeff2)) ALLOCATE( zeff2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(fpol)) ALLOCATE( fpol(pnx), STAT=istat)
      IF(.not.ALLOCATED(pres)) ALLOCATE( pres(pnx), STAT=istat)
      IF(.not.ALLOCATED(ffprim)) ALLOCATE( ffprim(pnx), STAT=istat)
      IF(.not.ALLOCATED(pprime)) ALLOCATE( pprime(pnx), STAT=istat)
      IF(.not.ALLOCATED(qpsi)) ALLOCATE( qpsi(pnx), STAT=istat)
      IF(.not.ALLOCATED(psiaux)) ALLOCATE( psiaux(pnx,2*pnz),            &  
     &                                 STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : wrpest  ' 
!============      
!
      if(ineg.ne.0) return
!
      write(nout,1000) kcycle,times,ipest
 1000 format("**** pest file written , cycle,times,ipest",               &  
     &        i7,1pe12.4,i3,"  ****")
!
      if(ipest.eq.2) then
      call lhwrite(neqplt,1)
      return
      endif
!
      if(ipest.eq.1) go to 200
!
      nxpest = nx
      nzpest = nz
      if(isym.eq.1) nzpest = 2*nz-1
      alrpest = alx-ccon
      alzpest = 2*alz
      rzpest = ccon
      rpest = gzero
      pmpest = psimin*tpi
      plpest = psilim*tpi
!
      dps = (psilim-psimin)/(nxpest-1)
!
      do 39 n=1,nxpest
      ps = psimin + dps*(n-1)
      call peval(ps,1+isurf,pval,ppval,imag,jmag)
      call geval(ps,1+isurf,gval,gpval,gppval,imag,jmag)
      ppest(n) = pval
   39 gpest(n) = gval/gzero
!
!
!     leave off title code,  implies running mapper with nocase=1
! ==> write (neqplt) (name(i),i=1,8),(name(i),i=1,8),idate
      write (neqplt) nxpest,nzpest,alrpest,alzpest,rpest,rzpest
!
      do 50 i=1,nxpest
      if(isym.eq.0) go to 45
      jpest = 0
      do 30 j=3,nzp
      jj = nzp+3-j
      jpest = jpest + 1
      opest(jpest,i) = psi(i+1,jj)*tpi
   30 continue
      go to 46
   45 jpest = 0
   46 continue
      do 40 j=nh,nzp
      jpest = jpest + 1
   40 opest(jpest,i) = psi(i+1,j)*tpi
   50 continue
      write(neqplt) ((opest(j,i),j=1,nzpest),i=1,nx)
!
      write (neqplt) xmag,zmag,pmpest,plpest
!     write (neqplt) alphap,alphag,p0,gp1
      write (neqplt) times,apl,p0,gp1
      write (neqplt) (gpest(i),i=1,nxpest)
      write (neqplt) (ppest(i),i=1,nxpest)
!
!.....special output to be processed by j-solver equilibrium code
!
!
  200 continue
      if(isurf.ne.1) return
 
      write(neqdsk) kcycle,isym,ipest,npsit,kmax
      write(neqdsk) times,xmag,zmag,gzero,apl,beta,betapol,ali2,         &  
     &              qsaw,psimin,psilim
      do 100 j=1,npsit
      ps = xsv2(j)
      call peval(ps,1+isurf,pval,ppval,imag,jmag)
      prpest2(j) = pval
      pppest2(j) = ppval
      te2(j)    = .5_R8*(te(j)   +te(j+1))
      ane2(j)   = .5_R8*(ane(j)  +ane(j+1))
      ti2(j)    = .5_R8*(ti(j)   +ti(j+1))
      avez2(j)  = .5_R8*(avez(j) +avez(j+1))
      zeff2(j) = .5_R8*(zeffa(j)+zeffa(j+1))
      if(j.eq.1 .or. j.eq.npsit) go to 100
      t1 = rdpsi*tpi*(qprof2(j))**2/(xmja2(j))**2
      t3 = gxmja(j+1)*xmja(j+1)/(.5_R8*(qprof2(j)+qprof2(j+1)))
      t5 = gxmja(j)*xmja(j)/(.5_R8*(qprof2(j)+qprof2(j-1)))
      ajpest2(j) = t1*(t3-t5)
  100 continue
      ajpest2(1) = 2._R8*ajpest2(2) - ajpest2(3)
      ajpest2(npsit) = 0._R8
      te2(1)    = 1.5_R8*te(2)    - 0.5_R8*te(3)
      ane2(1)   = 1.5_R8*ane(2)   - 0.5_R8*ane(3)
      ti2(1)    = 1.5_R8*ti(2)    - 0.5_R8*ti(3)
      avez2(1)  = 1.5_R8*avez(2)  - 0.5_R8*avez(3)
      zeff2(1) = 1.5_R8*zeffa(2) - 0.5_R8*zeffa(3)
      write(neqdsk) (prpest2(j),j=1,npsit)
      write(neqdsk) (pppest2(j),j=1,npsit)
      write(neqdsk) (ajpest2(j),j=1,npsit)
      write(neqdsk) (xsv2(j),j=1,npsit)
      write(neqdsk) (xplot(1,k),k=1,kmax+1)
      write(neqdsk) (zplot(1,k),k=1,kmax+1)
!
!......write equilibrium data in ASCII form (changed 9/22/02)
      write(neqdska,2001) (name(i),i=1,8)
 2001 format(20x,8a10)
      write(neqdska,6100) kcycle,isym,ipest,npsit,kmax
      write(neqdska,6101) times,xmag,zmag,gzero,apl,beta,                &  
     & betapol,ali2,qsaw,psimin,psilim
      write(neqdska,6101) (prpest2(j),j=1,npsit)
      write(neqdska,6101) (pppest2(j),j=1,npsit)
      write(neqdska,6101) (ajpest2(j),j=1,npsit)
      write(neqdska,6101) (xsv2(j),   j=1,npsit)
      write(neqdska,6101) (xplot(1,k),k=1,kmax+1)
      write(neqdska,6101) (zplot(1,k),k=1,kmax+1)
      write(neqdska,6101) (te2(j),    j=1,npsit)
      write(neqdska,6101) (ane2(j),   j=1,npsit)
      write(neqdska,6101) (ti2(j),    j=1,npsit)
      write(neqdska,6101) (avez2(j),  j=1,npsit)
      write(neqdska,6101) (zeff2(j), j=1,npsit)
 6100 format(5i10)
 6101 format(1p5e20.12)
      if(ncycle.eq.0) then
      if( numargs .lt. 1 ) then
         filename = 'boundary'
      else
         filename = 'boundary' // '.' // trim(suffix)
      end if
      open(unit=46,file=trim(filename),status='unknown')
      do k=2,kmax
      write(46,460) xplot(1,k),zplot(1,k)
 460  format(2f20.4)
      enddo
      do k=3,kmax-1
      write(46,460) xplot(1,kmax+2-k),-zplot(1,kmax+2-k)
      enddo
      endif
!.....write GEQDSK ASCII file in the format used by EFIT (added 6/04/03)
 2000 format(6a8,3i4)
 2020 format(5e16.9)
 2022 format(2i5)
      rdim = alx - ccon
      zdim = 2*alz
      zmid = 0._R8
      bcentr = gzero/xplas
      xdum = 0._R8
      idum = 13
      nzdim = nz + isym*(nz-1)
        do i=2,nx+1
        psiaux(i,j) = psi(i,j)
          if(isym.ne.0) then
             do j=3,nz+1
               psiaux(i,j+nz-1) = psi(i,j)
               psiaux(i,nz+3-j) = psi(i,j)
              enddo
               psiaux(i,nz+1) = psi(i,2)
          else
             do j=2,nz+1
               psiaux(i,j) = psi(i,j)
             enddo
          endif
        enddo
      write(neqdskg,2000) (name(i),i=1,6),idum,nx,nzdim
      write(neqdskg,2020) rdim,zdim,xplas,ccon,zmid
      write(neqdskg,2020) xmag,zmag,psimin,psilim,bcentr
!
!     global(75) is beta based on vacuum toroidal field
!     global(76) is central beta x 0.5
!     global(174) is troyon coefficient
!
      zip = tcurdtp*tpi*udsi
      write(neqdskg,2020) zip,psimin,global(75),xmag,global(76)*2._R8
      write(neqdskg,2020) zmag,global(174),psilim,xdum,xdum
      do 250 i=1,nx
      psival = psimin + (i-1)*(psilim-psimin)/(nx-1)
      call geval(psival,2,fpol(i),gpval,ffprim(i),imag,jmag)
      call peval(psival,2,pres(i),pprime(i),imag,jmag)
      do 249 jj=2,npsit
      j = jj
      if(psival .ge. xsv2(j-1) .and. psival .le. xsv2(j)) go to 251
 249  continue
 251  continue
      qpsi(i) = qprof2(j-1) + (psival-xsv2(j-1))                         &  
     &        /(xsv2(j)-xsv2(j-1))*(qprof2(j)-qprof2(j-1))
 250  continue
!
!.....temporary patch
      pprimesav2 = pprime(2)
      pprimesav1 = pprime(1)
      ffprimsav2 = ffprim(2)
      ffprimsav1 = ffprim(1)
!
      pprime(2) = 2._R8*pprime(3) - pprime(4)
      pprime(1) = 2._R8*pprime(2) - pprime(3)
      ffprim(2) = 2._R8*ffprim(3) - ffprim(4)
      ffprim(1) = 2._R8*ffprim(2) - ffprim(3)
!     write(nterm,6666) pprimesav1,pprimesav2,ffprimsav1,ffprimsav2
!     write(nterm,6666) pprime(1) ,pprime(2) ,ffprim(1) ,ffprim(2)
!6666 format("***diag***  ",1p4e12.4)
!
!.....end of temporary patch
      write(neqdskg,2020) (fpol(i),i=1,nx)
      write(neqdskg,2020) (udsp*pres(i),i=1,nx)
      write(neqdskg,2020) (ffprim(i),i=1,nx)
      write(neqdskg,2020) (udsp*pprime(i),i=1,nx)
      write(neqdskg,2020) ((psiaux(i,j),i=2,nxp),j=2,nzdim+1)
      write(neqdskg,2020) (qpsi(i),i=1,nx)
      nwalliter = 41
      nlimiter = kmax
      if(isym.eq.1) then
        nlimiter = 2*(kmax-2)
          do kk=kmax+2,nlimiter
          k = kk-kmax
          xplot(1,kmax+k) = xplot(1,kmax-k)
          zplot(1,kmax+k) =-zplot(1,kmax-k)
          enddo
        endif
!

      if (acoef(1).eq.2 .or. acoef(1).eq.6) then !NSTX
         write(neqdskg,2022) nlimiter+1,nwallnstx
         write(neqdskg,2020) (xplot(1,k),zplot(1,k),k=1,nlimiter+1)
         write(neqdskg,2020) (rwallnstx(i),zwallnstx(i),i=1,nwallnstx)
      else !ITER
         write(neqdskg,2022) nlimiter,nwalliter
         write(neqdskg,2020) (xplot(1,k),zplot(1,k),k=1,nlimiter)
         write(neqdskg,2020) (rwalliter(i),zwalliter(i),i=1,nwalliter)
      end if
!.....added 10/13/2011 for M3D to import profiles
      write(neqdskg,2020) ((ajphi(i,j),i=2,nxp),j=2,nzdim+1)
!
!.....special debug output for Josh and Hank
       do i=2,nxp
       do j=2,nx
         psival = psiaux(i,j)
         call geval(psival,2,fpol(i),gpval,ffprim(i),imag,jmag)
         call peval(psival,2,pres(i),pprime(i),imag,jmag)
         ajphi2 = -(xary(i)*pprime(i) + ffprim(i)/xary(i))
         write(201,2021) i,j,psival,ffprim(i),pprime(i),ajphi2,ajphi(i,j)
  2021   format(2i5,1p5e12.4)
       enddo
       enddo


      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
