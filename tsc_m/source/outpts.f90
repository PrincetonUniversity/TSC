      subroutine outpts
!......6.70 outpts
!......6.70 outpts
!                              R. O. S.     29 Jul 1987
!
!.....special output called every cycle
!
      USE CLINAM
      USE SAPROP
      USE SCR3
      USE NEWPLOT

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ntrite,iprs,lparm,j,ii,i,l,i2,n,jj,jmin,jmax
      INTEGER imax
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 adii,preskpa,t1,t3,t4
      REAL*8 t2,dum1,atuv,t6,t8,t9,t7,dum2,atu2,ajmax,ellip,delta
      REAL*8 x1,x2,z1,z2,currfl,qzer,qedge,aloopv,zmore,residd
      REAL*8 residmax
      REAL*8 rho,hnbe,hnbi,hice,hici,hece,heci,hlhe,hlhi
      REAL*8 hale,hali,ajnb,ajic,ajec,ajlh,ajbs,ajto,sab
      REAL*8 ajt,abs,hbre,hcyc,hlin
!============
!     common /newplot/ chienca(ppsi),chiinca(ppsi),
!    1                 chiicopi(ppsi),chiecopi(ppsi),diffary(ppsi)
!     dimension gsumg(pngroup),grrsumg(pngroup),grho(pngroup),
!    1  grsumg(pngroup),gzsumg(pngroup)
!    1          ,grho2(pngroup)
!     dimension icurpr(pnx)
!
!
      character*2 ibl,ichc,idot,inum(0:10),ichar(pnx)
!
!
      data ibl/"  "/
      data inum/" 0"," 1"," 2"," 3"," 4"," 5"," 6"," 7"," 8",            &  
     &          " 9","10"/
      data ichc/" c"/
      data idot/" ."/
!
      data ntrite/0/
      data iprs/100/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grrsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gzsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grho2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: icurpr
!============      
      IF(.not.ALLOCATED(gsumg)) ALLOCATE( gsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grrsumg)) ALLOCATE( grrsumg(pngroup),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(grho)) ALLOCATE( grho(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grsumg)) ALLOCATE( grsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(gzsumg)) ALLOCATE( gzsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grho2)) ALLOCATE( grho2(pngroup), STAT=istat)
      IF(.not.ALLOCATED(icurpr)) ALLOCATE( icurpr(pnx), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : outpts  ' 
!============      
!
!.....print header line every 50 cycles
!
      nxm = nx-1
      if(iprnt.eq.2) go to 60
      if(kcycle.gt.ncycle) go to 60
      iprs = iprs + 1
      lparm = 50
      if(ncycle.le.300) lparm = 10
      if(iprs.lt.lparm) go to 70
   60 continue
      iprs = 0
!
      if(isurf.eq.0 .or. lrswtch.ne.0) go to 401
!------------------------------------------------------------------------
!
!.....surface averaged transport variables
      write(nout,1059) kcycle,times,dts
      do 300 j=1,npsit
      if(adi(j).ne.0) adii = 1._R8/adi(j)
      if(j .eq. 1) then 
      ajt = (gzero/xplas)*((gxmja2(2)-gxmja2(1))*rdpsi*udsi/tpi)
      abs = (gzero/xplas)*(.5*(ajavbs(2)+ajavbs(1))*udsi)
      else
      ajt = (gzero/xplas)*((gxmja2(j)-gxmja2(j-1))*rdpsi*udsi/tpi)
      abs = (gzero/xplas)*(.5*(ajavbs(j)+ajavbs(j-1))*udsi)
      endif
      preskpa = adp(j)/vpg(j)*udsp*1.E-3_R8
  300 write(nout,2059) j,ti(j),te(j),anhy(j),ane(j),chiisec(j),          &  
     &  ratioe(j),ratioi(j),preskpa,xsv2(j),qprof2(j),                   &  
     &       ibaloon(j),ajt,abs,rminora(j)
 1059 format(1h1," cycle",i7,"    time",1pe12.4,"     dt",1pe12.4,//,4x,  &  
     &                                                                   &  
     & "j      ti(ev)      te(ev)     ni(mks)     ne(mks)     chii",     &  
     &"     ratioe     ratioi     press(kPa)     xsv2      qprof2",      &  
     &"    ibal     jdbt      jdbbs    rminor")
      write(nout,1069) kcycle,times,dts
      do 900 j=1,npsit
  900 write(nout,2069) j,diffnema(j),diffary(j),ane(j), chiitima(j),     &  
     &  chiinca(j),chiicopi(j), ti(j), chietema(j), chienca(j),          &  
     &  chiecopi(j), te(j)
 1069 format(1h1," cycle",i7,"    time",1pe12.4,"     dt",1pe12.4,//,4x,  &  
     &                                                                   &  
     & "j   diffnema    D ~ 1/n    ne(mks)  chiitima   chiinca  ",       &  
     & "  chi ~ 1/n    ti(ev)   chietema    chienca   che ~ 1/n",        &  
     & "     te" )
!2059 format(i5,1p10e12.4,i3)
 2059 format(i5,1p10e12.4,i3,1p3e12.4)
 2069 format(i5,1p11e11.3)
!
!...output for ITER simulations
!
!     do 1919 jj=3,npsit-1
!     te(jj) = (te(jj-2)+te(jj-1)+te(jj)+te(jj+1)+te(jj+2))/5.0
!     ti(jj) = (ti(jj-2)+ti(jj-1)+ti(jj)+ti(jj+1)+ti(jj+2))/5.0
!     savea(jj) = (savea(jj-2)+savea(jj-1)+savea(jj)+savea(jj+1)         &
!    & +savea(jj+2))/5.0
!     savia(jj) = (savia(jj-2)+savia(jj-1)+savia(jj)+savia(jj+1)         &
!    & +savia(jj+2))/5.0
!     ajavbs(jj) = (ajavbs(jj-2)+ajavbs(jj-1)+ajavbs(jj)+ajavbs(jj+1)         &
!    & +ajavbs(jj+2))/5.0
!     chiisec(jj)=(chiisec(jj-2)+chiisec(jj-1)+chiisec(jj)+chiisec(jj+1)         &
!    & +chiisec(jj+2))/5.0
!     chiesec(jj)=(chiesec(jj-2)+chiesec(jj-1)+chiesec(jj)+chiesec(jj+1)         &
!    & +chiesec(jj+2))/5.0
!1919 continue
      write(nout,9069) kcycle,times,dts
      write(nout,9271)
      write(nout,9272)
      write(nout,9169)
      do 904 j=2,npsit
      rho = sqrt(float(j-1)*dpsi/(float(npsit-1)*dpsi))
      hnbe = savee(j)*udsp/udst
      hnbi = savei(j)*udsp/udst
      hice = savefw(j)*udsp/udst
      hici = savifw(j)*udsp/udst
      hece = savebm(j)*udsp/udst
      heci = savibm(j)*udsp/udst
      hlhe = savelh(j)*udsp/udst
      hlhi = savilh(j)*udsp/udst
  904 write(nout,9269) rho,te(j),ti(j),hnbe,hnbi,hice,hici,hece,heci,   &
     &hlhe,hlhi
      write(nout,9271)
      write(nout,9272)
      write(nout,9170)
      do 901 j = 2, npsit
      rho = sqrt(float(j-1)*dpsi/(float(npsit-1)*dpsi))
      hale = savea(j)*udsp/udst
      hali = savia(j)*udsp/udst
      ajnb = 0.5*(ajavcd(j)+ajavcd(j-1))*udsi
      ajic = 0.5*(ajavfw(j)+ajavfw(j-1))*udsi
      ajec = 0.5*(ajavec(j) + ajavec(j-1))*udsi
      ajlh = 0.5*(ajavlh(j)+ajavlh(j-1))*udsi
      ajbs = 0.5*(ajavbs(j)+ajavbs(j-1))*udsi
      ajto = (gxmja2(j)-gxmja2(j-1))*rdpsi*udsi/tpi
      sab = gzero/rmajora(j)
  901 write(nout,9269) rho,hale,hali,ajnb,ajic,ajec,ajlh,ajbs,ajto,     &
     &qprof2(j),sab
      write(nout,9271)
      write(nout,9272)
      write(nout,9173)
      do 905 j = 2, npsit
      rho = sqrt(float(j-1)*dpsi/(float(npsit-1)*dpsi))
      hbre = savebre(j)*udsp/udst
      hcyc = savecyc(j)*udsp/udst
      hlin = (sradion(j) + saveimp(j)*udsp/udst)
  905 write(nout,9269) rho,chiisec(j),chiesec(j),hbre,hcyc,hlin,        &
     & ane(j),anhe(j),fraci(4)*ane(j),fraci(7)*ane(j),zeffa(j)
 9069 format(1h1," cycle",i7,"    time",1pe12.4,"     dt",1pe12.4)
 9169 format(                                                           &
     & "rho  Te        Ti        HNBe      HNBi      HICe      ",       &
     & "HICe      HECe      HECi      HLHe      HLHi")
 9170 format(                                                           &
     & "rho  HALe      HALi      JNB       JIC       JEC       ",       &
     & "JLH       Jbs       Jtot      q         <Bphi>")
 9173 format(                                                           &
     & "rho  ChiI      ChiE      Hbre      Hcyc      Hlin      ",       &
     & "ne        nHe       nBe       nAr       Zeff  ")
 9269 format(f5.3,1p10e10.3)
 9271 format(//)
 9272 format("rho = sqrt(phi/phi_b), H___ is W/m3, J__ is A/m2-T",      &
     &"T is eV",/,                                                      &
     &"AL-alpha,NB-neutral beam,IC-ion cyc,EC-electron cyc,"            &
     &"LH-lower hybrid,bs-bootstrap,tot-total,bre-bremsstrahlung,"      &
     &"cyc-cyclotron,lin-line,He-helium,Be-beryllium,Ar-argon")
!
!...end of ITER output
!                        skip coil print if coil plot skipped (ROS, Mar 88)
!.....coil currents
  401 if(nwire.eq.0 .or. noplot(5).gt.0)   go to 403
      write(nout,3500)
 3500 format(//'    i xwire zwire cwire(ka) cwire0(ka) diff',            &  
     &       '  volts  t-f1     ',1x,                                    &  
     &       '    j xwire zwire cwire(ka) cwire0(ka) diff',              &  
     &       '  volts  t-f1  ')
      do 85 ii=nc0,nwire,2
      i = ncoil-nwire+ii
      t1 = ccoil(i)*udsi*1.E-3_R8
      t3 = cwire0(ii)*udsi*1.E-3_R8
      t4 = tpi*resave(ii)*udsv
      t2 = t3-t1
      dum1 = 0._R8
      atuv = 0._R8
      do 186 l=1,ntpts
  186 atuv = atuv + atnvw1(l,ii)*fact(l)
      i2 = ii + 1
      j = ncoil-nwire+i2
      t6 = ccoil(j)*udsi*1.E-3_R8
      t8 = cwire0(i2)*udsi*1.E-3_R8
      t9 = tpi*resave(i2)*udsv
      t7 = t8-t6
      dum2 = 0._R8
      atu2 = 0._R8
      do 187 l=1,ntpts
  187 atu2 = atu2 + atnvw1(l,i2)*fact(l)
      write(nout,3600) ii,xcoil(i),zcoil(i),t1,t3,t2,t4,atuv,            &  
     &              i2,xcoil(j),zcoil(j),t6,t8,t7,t9,atu2
 3600 format(1x,i5,2f6.2,2f9.2,f8.2,f8.2,f6.2,5x,                        &  
     &          i5,2f6.2,2f9.2,f8.2,f8.2,f6.2)
   85 continue
  403 continue
!------------------------------------------------------------------------
      if(ncoil-nwire .eq. 0) go to 407
      write(nout,3601)
 3601 format(//"    i xcoil zcoil ccoil(ka) ")
      do 86 n=nc0,ncoil-nwire
      t1 = ccoil(n)*udsi*1.E-3_R8
   86 write(nout,3600) n,xcoil(n),zcoil(n),t1
!  skip over printer graphs except last cycle
      if(kcycle.le.ncycle.or.nx.gt.68) go to 407
      ajmax = 0._R8
!
!.....printer graph of current ,boundary, and coils
      do 402 i=2,nxp
      do 402 j=2,nzp
  402 ajmax = max(ajmax,ajphi(i,j))
      if(ajmax .le. 0) go to 407
      write(nout,1403) kcycle,ajmax
 1403 format('1 max current at cycle',i7,1pe12.4)
      do 405 jj=2,nzp
      j = 2 + nzp - jj
      do 406 i=3,nx
  406 icurpr(i) = int(10._R8*ajphi(i,j)/ajmax)
      do 409 i=3,nx
  409 ichar(i) = ibl
      if(j.lt.jminn.or.j.gt.jmaxx) go to 419
      do 410 i=iminn,imaxx
  410 ichar(i) = inum(icurpr(i))
  419 continue
      do 420 ii=nc0,nwire
      if(jwire(ii).ne.j) go to 420
      ichar(iwire(ii)) = ichc
  420 continue
      do 425 i=3,nx-1
      if(psi(i,j).ge.psilim .and. psi(i+1,j).le.psilim) go to 426
      if(psi(i,j).le.psilim .and. psi(i+1,j).ge.psilim) go to 427
      go to 425
  426 continue
      ichar(i) = idot
      go to 425
  427 ichar(i+1) = idot
  425 continue
      write(nout,1404) j,(ichar(i),i=4,nx-1)
 1404 format(1x,0p,i2,64a2)
  405 continue
!
!------------------------------------------------------------------------
  407 continue
!
!.....calculate and printout shape,leftmost,rightmost, top and bottom
      if(isurf.le.0) then
      call shape(ellip,delta,x1,x2,z1,z2)
      write( nout,3086) kcycle,times,ellip,x2,x1,z1,z2
!     write(nterm,3086) kcycle,times,ellip,x2,x1,z1,z2
 3086 format(/" cycl       time     ellip      left     right       top  &  
     &   bottom"/   i5,1pe11.3,1p5e10.2)
      endif
      if(irfp.ne.1) then
      write( nout,168) rmajor,rminor,shape3, shape5,shape6,shape7
      if(acoef(901).eq.0._R8) write(nterm,168) rmajor,rminor,shape3,     &  
     &shape5,shape6,shape7
  168 format(/"         r0        a     delt      eps     xsep     zsep"  &  
     &       /2x,6f9.5)
      if(numsep.ge.2) then
      do 7168 n=2,numsep
      write( nout,6168) xsep(n),zsep(n)
      write(nterm,6168) xsep(n),zsep(n)
 7168 continue
 6168 format(38x,2f9.5)
                      endif
      else
      write( nout,1168) rmajor,rminor,shape3, shape5,shape6,shape7
      if(acoef(901).eq.0._R8) write(nterm,1168) rmajor,rminor,shape3,    &  
     &shape5,shape6,shape7
 1168 format(/"         r0        a     psi1     psi2     psi3           &  
     & psi4"                                                             &  
     &  /2x,6f9.5)
      endif
      fluxlink = 0._R8
      call fieldg(xplas,zplas,gsumg,grsumg,gzsumg,                       &  
     &            grrsumg,grho,grho2)
      do 1169 l=1,pngroup
      currfl = 0._R8
      do 297 n=1,ntpts
  297 currfl = currfl +fact(n)*gcur(n,l)*usdi
 1169 fluxlink = fluxlink+gsumg(l)*tpi*(gcurfb(l)+currfl)
      if(acoef(901).eq.0._R8) then
      if(isurf.eq.1) write( nout,169) del95,el95,fluxlink
      if(isurf.eq.1) write(nterm,169) del95,el95,fluxlink
  169 format( "     95% flux surf: ",2f9.5,"  flux linkage ",            &  
     &          1pe12.4," Webers")
!
      qzer = qprof(1)
      qedge = qprof(npts)
      if(isurf.ne.0) qzer = qprof2p
      if(isurf.ne.0) qedge = qprof2(npts)
      write(nout,3087) qzer,qedge,betapol,ali2,aliga,vol,uint,npsit
 3087 format(/'          q1          qe     betapol        li/2  li(GA d  &  
     &ef)        vol        uint  npsit'/1p7e12.4,i4)
      endif
      if(isurf.eq.0 .or. ibalsw .le.0) go to 670
      write(nout,1670)
 1670 format(" balloon stability: (0-stable,1-bal unstable,2-mercier"    &  
     &," unstable 3-both")
      do 660 j=1,npsit,20
      jmin = j
      jmax = j+19
      if(jmax.gt.npsit) jmax = npsit
      write(nout,1680) (ibaloon(jj),jj=jmin,jmax)
 1680 format(10x,20i3)
  660 continue
  670 continue
      write(nout,4000)
      if(acoef(901).eq.0._R8)  write (nterm, 4020)
      ntrite = 0
!------------------------------------------------------------------------
   70 continue
      aloopv = tpi*reboun*udsv
!
!.....write 1 line cycle printout
!
      zmore = zmag
      if(isym.eq.1) zmore = gzero
      call residcal(residd,imax,jmax,residmax)
      write ( nout,3000) kcycle,times,dts,tauems ,idtmin,jdtmin,tauekg,  &  
     &                   zmore,ekin,apl,amach,xmag,betat0,aloopv,iplim
 
!
!                                       Output to terminal every other print out
       ntrite = ntrite + 1
       if (ntrite.ne.1)   return
       ntrite = - 4
      if(acoef(901).eq.0._R8) write (nterm,3020)  kcycle,times, idtmin,  &  
     &jdtmin,ekin,apl,xmag,zmag,betat0,iplim,resid,ali2,global(141),     &  
     & global(213)
      return
 3000 format(i7,1pe11.4,1pe9.2, 0pf6.1,2i3,1p6e11.3,1p2e10.2,i4)
 3020 format(i7,1pe11.4, 2i3, 1p2e11.3, 0pf8.3,f7.3, 1pe9.2, i4,         &  
     &      1p4e9.2)
!
 4000 format(1h1,"  cyc...  time      dt taue(ms)  i  j     tauekg    ",  &  
     &                                                                   &  
     &"  zmag(gzero)  ekin     pl cur      amach      xmag       beta ",  &  
     &                                                                   &  
     &"loop v  iplim")
 4020 format(   "   cyc       time  i  j       ekin     pl cur    xmag   &  
     &  zmag    beta iplim  resid    ali2     boot     tot")
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
