      subroutine runaway_sub
!
!
!.....advance runaway electron currents
!
      USE CLINAM
      USE RUNAWAY
      USE SAPROP
      USE SCADVAN
      USE SCR1
      USE SCR4
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!--------
!       ajpre   - cell vertice runaway electron current density (advan23)
!       ajprecc - cell center runaway electron current density
!       ajphisf - tsc surface average current density
!       ajpresf - surface average runaway current density
!       anre    - surface average runaway density
!       adnre   - anre*vp(j)
!       recur   - runaway current
!       sresf   - surface average source term
!       sreav   - avalanche production
!       sumre   - total runaway production
!--------
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,kloop,i,jplas,k1,k,k2,m,n
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cntr,epso,ameo,cfjev,cfevj,qe,bump,aratio,c,dtre
      REAL*8 trmo,gameps,sumre1,anreo,adnreo,areatsc,etamid,ef,trma
      REAL*8 amvc2,fnew,dfda,dela,denom,efc,zfac,denomo,ecre,taure
      REAL*8 avseed,sre,extrps,srel,exptrm,trmld,eovrec,avtrm,trmla
      REAL*8 sreh,trmhd,trmha,trmh,trml,dfds,patrm,dsre,sresgo
      REAL*8 sresgn,avajp2,avajre,etamlt,sumre2,xmin,xmax,zmin,zmax
      REAL*8 dlzmx,dlxmx,vzmn,vxmn,vzmx,vxmx,ajpremx,ajpremn,psimx
      REAL*8 psimn,dcntr,ecmin1,ecmax1,ecmin2,ecmax2,xset,yset
      REAL*8 temax,ajtscmax,ajremax,etamax
!============
!     dimension ajphisf(ppsi), yplot(pnx)  , zarym(pnz), efs(ppsi)
!    .,         efcd(ppsi)   , efca(ppsi)  , xs(ppsi)  , etamidy(ppsi)
!    .,         etsc(ppsi)   , ajtsc(ppsi)
      dimension cntr(30)
!
      REAL*8 muo, nuo, Lambda
!
      data epso / 8.854E-12_R8/
      data ameo /  0.91E-30_R8/
      data cfjev / 6.242E18_R8/
      data cfevj / 1.6021E-19_R8/
      data qe / 1.6021E-19_R8/
      data muo / 1.2566E-06_R8/
      data bump / 1.E-5_R8/
      data aratio / 0.25_R8/
      data c / 3.E8_R8/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ajphisf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: yplot
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zarym
      REAL*8, ALLOCATABLE, DIMENSION(:) :: efs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: efcd
      REAL*8, ALLOCATABLE, DIMENSION(:) :: efca
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: etamidy
      REAL*8, ALLOCATABLE, DIMENSION(:) :: etsc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ajtsc
!============      
      IF(.not.ALLOCATED(ajphisf)) ALLOCATE( ajphisf(ppsi), STAT=istat)
      IF(.not.ALLOCATED(yplot)) ALLOCATE( yplot(pnx), STAT=istat)
      IF(.not.ALLOCATED(zarym)) ALLOCATE( zarym(pnz), STAT=istat)
      IF(.not.ALLOCATED(efs)) ALLOCATE( efs(ppsi), STAT=istat)
      IF(.not.ALLOCATED(efcd)) ALLOCATE( efcd(ppsi), STAT=istat)
      IF(.not.ALLOCATED(efca)) ALLOCATE( efca(ppsi), STAT=istat)
      IF(.not.ALLOCATED(xs)) ALLOCATE( xs(ppsi), STAT=istat)
      IF(.not.ALLOCATED(etamidy)) ALLOCATE( etamidy(ppsi), STAT=istat)
      IF(.not.ALLOCATED(etsc)) ALLOCATE( etsc(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ajtsc)) ALLOCATE( ajtsc(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : runaway  ' 
!============      
!
      dtre = dts
!
! * * estimate runaway production
!
      if(recur.lt.apl .and. igone.eq.0)  then
!
! * * map toroidal current density on to TSC 1D grid
!
!cc      call map21(ajphi, ajphisf)
!
      trmo   = qe**3/(2._R8*tpi*epso**2)
      gameps = 1._R8/(1._R8+1.46_R8*sqrt(aratio)+1.72_R8*aratio)
      sumre1 = 0._R8
!
      do 50 j = 2,npsit
      anreo      = anre(j)
      adnreo     = adnre(j)
      areatsc = (vary(j) - vary(j-1)) / xplas
      ajtsc(j) = (gxmja2(j) - gxmja2(j-1)) * udsi / areatsc
      ajphisf(j) = ajtsc(j)
!
! * * * calculate critical velocity given plasma electric field
! * * *     (Equ. 2.9 - H. Knoepfel & D. A. Spong, Nuclear Fusion
!                 Vol. 19, No. 6, pp 785-829 (1979) )
!
      etamid  = .5_R8*(etpara(j)+etpara(j-1))*udsr
      etamidy(j) = etamid
      ef = etamid*(ajphisf(j)-qe*anreo*3.E8_R8)
      ef = max(1.E-2_R8,ef)
      trma  = 2._R8*trmo*ane(j)*cfjev/ef
      amvc2 = 24._R8*trma
      kloop = 0
   60 kloop = kloop + 1
      Lambda = 24._R8-log(sqrt(ane(j)*1.E-6_R8)/amvc2)
      fnew   = trma*Lambda - amvc2
      dfda   = trma*(1._R8/amvc2) - 1._R8
      dela   = - fnew/dfda
      if(abs(dela).gt.0.5_R8*amvc2) dela = 0.5_R8*sign(1.0_R8,dela)*     &  
     & amvc2
      amvc2 = amvc2 + dela
!
      if(amvc2.ne.0.0_R8) dela = dela/amvc2
      if(abs(dela).gt.1.E-6_R8.and. kloop.lt.30) go to 60
      if(amvc2.le.3._R8*te(j)) amvc2 = 3._R8*te(j)
      Lambda = 24._R8-log(sqrt(ane(j)*1.E-6_R8)/amvc2)
!
      if(kloop.ge.30) then
      write(nout,*) '****** runaway did not converge amvc2 ******'
      write(nout,238) times,i,j,amvc2,dela
  238 format(10x,'time=',1pe11.4,' i=',i3,' j=',i3,' amvc2=',1pe11.4     &  
     &,   ' dela=',1pe11.4)
               endif
!
! * * * calculate critical field
! * * *     (Equ. 2.12 - H. Knoepfel & D. A. Spong )
!
      denom   = te(j)/cfjev
      efc     = max(ef,trmo*ane(j)*Lambda/denom)
      efcd(j) = efc
!
! * * calculate runaway source term
! * * *     (Equ. 2.16 - H. Knoepfel & D. A. Spong )
!
      if(te(j).le.10.0_R8) then
      Lambda = 23._R8-log(sqrt(ane(j)*1.E-6_R8)/te(j)**1.5_R8)
            else
      Lambda = 24._R8-log(sqrt(ane(j)*1.E-6_R8)/te(j)     )
            endif
      nuo    =  trmo*qe*(ane(j)*Lambda)                                  &  
     &           /(te(j)/cfjev)**1.5_R8/sqrt(ameo)
      zfac   = 0.556_R8+.36_R8*zeffa(j)-.099_R8*zeffa(j)**2+.0064_R8*    &  
     & zeffa(j)**3
      denomo = sqrt(0.19_R8+4._R8*pi*(zeffa(j)+1._R8)**2                 &  
     &       / (3._R8*gameps*(zeffa(j)+5._R8)*(0.51_R8+4._R8/gameps**2))  &  
     & )
      ecre   = ane(j)*qe**3*Lambda/(4._R8*pi*epso**2*ameo*c**2)
      efca(j)= ecre
      taure  = ameo*c/(qe*ecre)
      avseed = 8.6E-23_R8*ane(j)/(qe*c)
      kloop  = 0
      sre    = sresf(j)
      extrps = 1.0_R8
   70 kloop  = kloop + 1
      srel   = sre * (1._R8- bump)
      if(srel.eq.0._R8) srel = -1._R8
      ef = etamid*(ajphisf(j)-qe*c*(anreo + srel*dtre))
      ef = max(1.E-2_R8,ef)
      exptrm = min(100.0_R8,efc/(4._R8*ef)+sqrt((zeffa(j)+1.0_R8)*efc/   &  
     & ef))
      trmld  = zfac*(efc/ef)**(3._R8*(zeffa(j)+1.0_R8)/16._R8)           &  
     &                            *ane(j)*nuo*exp(-exptrm)
      eovrec = max(1.0_R8,ef/ecre)
      denom  = sqrt(1._R8-1._R8/eovrec+4._R8*pi*(zeffa(j)+1._R8)**2      &  
     &       / (3._R8*gameps*(zeffa(j)+5._R8)*(eovrec**2+4._R8/gameps**  &  
     & 2-1._R8)))
      avtrm  = (1._R8/(taure*Lambda))*sqrt(pi*gameps/(3._R8*(zeffa(j)+5)  &  
     & ))                                                                &  
     &       *((eovrec-1._R8)/denom-0.23_R8/denomo)
      trmla  =  (anreo + srel*dtre)*avtrm
!
      sreh    = sre * (1._R8+ bump)
      if(srel.eq.0._R8) srel = 1._R8
      ef = etamid*(ajphisf(j)-qe*c*(anreo + sreh*dtre))
      ef = max(1.E-2_R8,ef)
      exptrm = min(100.0_R8,efc/(4._R8*ef)+sqrt((zeffa(j)+1.0_R8)*efc/   &  
     & ef))
      trmhd = zfac*(efc/ef)**(3._R8*(zeffa(j)+1.0_R8)/16._R8)            &  
     &                            *ane(j)*nuo*exp(-exptrm)
      eovrec = max(1.0_R8,ef/ecre)
      denom  = sqrt(1._R8-1._R8/eovrec+4._R8*pi*(zeffa(j)+1._R8)**2      &  
     &       / (3._R8*gameps*(zeffa(j)+5._R8)*(eovrec**2+4._R8/gameps**  &  
     & 2-1._R8)))
      avtrm  = (1._R8/(taure*Lambda))*sqrt(pi*gameps/(3._R8*(zeffa(j)+5)  &  
     & ))                                                                &  
     &       *((eovrec-1._R8)/denom-0.23_R8/denomo)
      trmha  =  (anreo + sreh*dtre)*avtrm
!
      trmh = trmhd + trmha
      trml = trmld + trmla
      dfds = (trmh - trml)/(sreh - srel) - 1.0_R8
!
      ef = etamid*(ajphisf(j)-qe*c*(anreo + sre*dtre))
      ef = max(1.E-2_R8,ef)
      efs(j) = ef
      eovrec = max(1.0_R8,ef/ecre)
      denom  = sqrt(1._R8-1._R8/eovrec+4._R8*pi*(zeffa(j)+1._R8)**2      &  
     &       / (3._R8*gameps*(zeffa(j)+5._R8)*(eovrec**2+4._R8/gameps**  &  
     & 2-1._R8)))
      avtrm  = (1._R8/(taure*Lambda))*sqrt(pi*gameps/(3._R8*(zeffa(j)+5)  &  
     & ))                                                                &  
     &       *((eovrec-1._R8)/denom)
      patrm  = (1._R8/(taure*Lambda))*sqrt(pi*gameps/(3._R8*(zeffa(j)+5)  &  
     & ))                                                                &  
     &       *(-0.23_R8/denomo)
      exptrm = min(100.0_R8,efc/(4._R8*ef)+sqrt((zeffa(j)+1.0_R8)*efc/   &  
     & ef))
!          Dreicer contribution
      fnew = zfac*(efc/ef)**(3._R8*(zeffa(j)+1.0_R8)/16._R8)             &  
     &                            *ane(j)*nuo*exp(-exptrm) - sre         &  
!          Avalance contribution
     &     + (anreo + sre*dtre)*avtrm                                    &  
!          Additional pitch angle scattering loss
     &     + (anreo + sre*dtre)*patrm                                    &  
!          Putvinski avalanche source
     &     + avseed
!
      dsre = - fnew/dfds
      sresgo = sresgn
      sresgn = sign(1.0_R8,dsre)
      if(kloop.gt.5 .and. sresgn.ne.sresgo) extrps = 0.5_R8*extrps
      sre = sre + extrps*dsre
      if(sre.ne.0.0_R8) dsre = dsre/sre
      if(abs(dsre).gt.1.E-6_R8.and. kloop.lt.40) go to 70
!
      if(kloop.ge.40 .and. sre.gt.1.E10_R8) then
      write(nout,*) '****** runaway did not converge sre ******'
      write(nout,239) times,j,sre,dsre
  239 format(10x,'time=',1pe11.4,' j=',i3,' sre=',1pe11.4                &  
     &,   ' dsre=',1pe11.4)
      sre = 0._R8
                   endif
!
      sresf(j)   = sre
      adnre(j)   = max(0._R8,adnreo + dtre*sresf(j)*vp(j))
      anre(j)    = adnre(j)/vp(j)
!
!.....Change suggested by Merrill 9/4/98
      ajpresf(j) = min(ajtsc(j),c*qe*anre(j))
!.....added 09/10/98
      ajpresf(j) = max(0._R8,ajpresf(j))
!
      sumre      = sumre + dpsi*vp(j)*(sre-anre(j)*patrm)*dtre
      sreav      = sreav + dpsi*vp(j)*anre(j)*avtrm*dtre
      sumre1     = sumre1+ dpsi*vp(j)*anre(j)
!     write(nterm,1998) j, sresf(j), adnre(j), anre(j)
!1998 format("j=",I3," sre=",1pe12.4," adn=",1pe12.4," anr=",1pe12.4)
!     write(nterm,1999) ajpresf(j), sumre, sreav, sumre1
!1999 format("aj=",1pe12.4," su=",1pe12.4," sr=",1pe12.4," s1=",1pe12.4)
!     write(nterm,2000) kcycles, times
!2000 format("cycle=",I3," time=",1pe11.4)
!
   50 continue
!
!     write(nterm,*)'****** Comparison E and J with TSC ********'
!     do j=2,npsit
!     etsc(j) = as(j) * udsv / (tpi*xplas)
!     areatsc = (vary(j) - vary(j-1))/xplas
!     ajtsc(j) = (gxmja2(j) - gxmja2(j-1)) * udsi / areatsc
!     write(nterm,2050)j, efs(j), ajphisf(j)
!     write(nterm,2100)etamidy(j), etsc(j), ajtsc(j)
!2050 format("j=",I3," efs(j)=",1pe12.4," ajphisf(j)=",1pe12.4)
!2100 format("etamidy=",1pe12.4," etsc=",1pe12.4
!    ., " ajtsc(j)=",1pe12.4)
!     end do
                       endif
      do 55 j=2,npsit
      avajp2 = max(1.E-6_R8,0.5_R8*(ajphisf(j)+ajphisf(j+1)))
      avajre = 0.5_R8*(ajpresf(j)+ajpresf(j+1))
      etamlt = max(1.E-2_R8, (avajp2-avajre)/avajp2)
      etamlt = min(1.0_R8,etamlt)
      etafac(j) = 0.1_R8*(etamlt-etafac(j))+etafac(j)
   55 continue
      etafac(1) = etafac(2)
      efcd(1)   = efcd(2)
      efca(1)   = efca(2)
      efs(1)    = efs(2)
!
      do 65 j = npsit+1, npsi
      anre(j)    = 0._R8
      sresf(j)   = 0._R8
      ajpresf(j) = 0._R8
      etafac(j)  = 1.0_R8
   65 continue
!
      call map12(ajpresf, ajprecc)
!
      do 75 i = 2,nxp
      do 75 j = 2,nzp
   75 ajpre(i,j) = 0._R8
!
      recur  = 0._R8
      sumre2 = 0._R8
      do 80 i = iminn, imaxx
      do 85 j = jminn, jmaxx
!
      ajpre(i,j) = 0.25_R8*usdi*( ajprecc(i  ,j  ) + ajprecc(i+1,j  )    &  
     &           +             ajprecc(i+1,j+1) + ajprecc(i  ,j+1) )
!
      if(j.gt.jminn) then
      recur  = recur  + (1._R8+isym)*ajprecc(i,j)*deex*deez
      sumre2 = sumre2                                                    &  
     &       + (1._R8+isym)*ajprecc(i,j)*tpi*xarh(i)*deex*deez/(qe*c)
              endif
   85 continue
      if(isym.eq.1) ajpre(i,1) = ajpre(i,3)
   80 continue
!
      if(mod(kcycle,nskipr/2).eq.0) then
      write(nterm,200) times, recur, sreav/max(1.E-6_R8,sumre)
                   endif
  200 format(10x,"time=",1pe11.4," runaway current=",1pe11.4             &  
     &,          " avalanche fraction=",1pe11.4)
!
! * * plot contors
!
!---------
      if(mod(kcycle,nskipl).eq.0) then
!
      call colora("white")
      xmin = xary(iminn)
      xmax = xary(imaxx)
      zmin =-zary(jmaxx)
      zmax = zary(jmaxx)
      dlzmx = zmax - zmin
      dlxmx = xmax - xmin
      vzmn = .35_R8
      vxmn = .2_R8
      if(dlxmx.gt.dlzmx) then
      vzmx = .6_R8*(dlzmx/dlxmx) + vzmn
      vxmx = .8_R8
             else
      vxmx = .6_R8*(dlxmx/dlzmx) + vxmn
      vzmx = .95_R8
             endif
!
      ajpremx = -1.E10_R8
      ajpremn =  1.E10_R8
      do 192 i = 1,nx
      do 193 j = 1,nz
      if(ajpre(i,j).gt.0.0_R8) ajpremn  = min(ajpremn,ajpre(i,j))
      ajpremx  = max(ajpremx,ajpre(i,j))
  193 continue
      jplas    = (zplas-zzero)/deez + 2
      yplot(i) = ajpre(i,jplas)
  192 continue
      ajpremn = min(100._R8*ajpremn,ajpremx/100._R8)
      if(ajpremn.eq.0 .and. ajpremx.eq.0._R8) go to 195
!
      call maps(xmin, xmax, zmin, zmax, vxmn, vxmx, vzmn, vzmx)
!
      psimx = -1.E10_R8
      psimn =  1.E10_R8
      do 180 j = 1,nz
      zarym(j) = -zary(j)
  180 continue
!
      write(s100,181) times, cntr(2)
  181 format(" time(s)=",1pe11.4                                         &  
     &,      " contr interval=",1pe11.4)
      call setlch(xmin, 1.1_R8*zmax, 1, 2, 2, 0)
      call gtext(s100,80,0)
!
      k1 = 20
      dcntr = (psimx-psimn)/(k1-1)
      cntr(1) = psimn
      do 182 k = 2, k1-1
      cntr(k) = cntr(k-1) + dcntr
  182 continue
      k2 = k1
      cntr(k2) = cntr(k2-1) + 1.5_R8*dcntr
      call rcontr(k1,cntr,k2,psi,pnx,xary,1,nx,1,zary ,1,nz,1)
      if(isym.eq.1) then
      call rcontr(k1,cntr,k2,psi,pnx,xary,1,nx,1,zarym,1,nz,1)
      endif
!
      call setch(5._R8,35._R8,2,2,1,-1)
      write(s100,2007)
 2007 format('Height (m)')
      call gtext(s100,10,0)
      call setch(25._R8,4._R8,2,2,0,-1)
      write(s100,2017)
 2017 format('Major Radius (m)')
      call gtext(s100,17,0)
 
      call colora("red")
!
      cntr(1) = ajpremn
      cntr(2) = abs(ajpremx-ajpremn)/20._R8
!
      call rcontr(0,cntr,0,ajpre,pnx,xary,1,nx,1,zary ,2,nz,1)
      if(isym.eq.1) then
      call rcontr(0,cntr,0,ajpre,pnx,xary,1,nx,1,zarym,2,nz,1)
      endif
 
      call colora("white")
      call maps(xmin, xmax, ajpremn, ajpremx, vxmn, vxmx, .2_R8, vzmn)
!
      call colora("red")
      do 194 m = 2, nx
      call line(xary(m-1), yplot(m-1),xary(m), yplot(m))
  194 continue
 
      call colora("white")
      call frscj(1)
!
      ecmin1 = 1.E30_R8
      ecmax1 =-1.E30_R8
      ecmin2 = 1.E30_R8
      ecmax2 =-1.E30_R8
      xmin   = 1.E30_R8
      xmax   =-1.E30_R8
      do 300 j=2,npsit
      ecmax1 = max(ecmax1, efcd(j), efs(j))
      ecmax2 = max(ecmax2, efca(j), efs(j))
      ecmin1 = min(ecmin1, efcd(j), efs(j))
      ecmin2 = min(ecmin2, efca(j), efs(j))
      if(j.le.2) go to 300
      xs(j)  = sqrt((j-1.5_R8)/(npsit-1))
  300 continue
      xmax   = 1.0_R8
      xmin   = 0.0_R8
      ecmin1 = max(ecmin1,1.E-4_R8*ecmax1)
      ecmin2 = max(ecmin2,1.E-4_R8*ecmax2)
!
      call map(xmin, xmax, ecmin1, ecmax1, 0.15_R8, 0.95_R8, 0.20_R8,    &  
     & 0.7_R8)
      write(s100,325)
  325 format("Electric Fields (V/M)")
      xset = xmin - 0.10_R8*(xmax-xmin)
      yset = ecmin1 + 0.15_R8*(ecmax1-ecmin1)
      call setold(xset, yset, 1, 0, 1, 1)
      call gtext(s100,80,0)
      write(s100,326)
  326 format("r/a")
      xset = xmin + 0.4_R8*(xmax-xmin)
      yset = ecmin1 - 0.25_R8*(ecmax1-ecmin1)
      call setold(xset, yset, 1, 0, 1, 0)
      call gtext(s100,80,0)
!
      call setold(xmin, 1.3_R8*ecmax1, 1, 0, 1, 0)
      call colora("red")
      write(s100,428) times, kcycle
  428 format(" time = ",1pe10.2," sec,    cycle =",i7)
      call gtext(s100,80,0)
      call colora("green")
      write(s100,328)
  328 format(" D - Critical electric field for Dreicer source")
      call gtext(s100,80,0)
      call colora("magenta")
      write(s100,327)
  327 format(" A - Critical electric field for Avalanche source")
      call gtext(s100,80,0)
      call colora("cyan")
      write(s100,329)
  329 format(" E - Electric field")
      call gtext(s100,80,0)
!
      call colora("white")
      call mapgsl(xmin, xmax, ecmin1, ecmax1, 0.15_R8, 0.50_R8, 0.20_R8,  &  
     &  0.7_R8)
      call colora("green")
      call tracec("D",xs(2), efcd(2), npsit-1,-1,-1,0._R8,0._R8)
      call colora("cyan")
      call tracec("E",xs(2), efs(2) , npsit-1,-1,-1,0._R8,0._R8)
      call colora("white")
      call mapgsl(xmin, xmax, ecmin2, ecmax2, 0.60_R8, 0.95_R8, 0.20_R8,  &  
     &  0.7_R8)
      call colora("magenta")
      call tracec("A",xs(2), efca(2), npsit-1,-1,-1,0._R8,0._R8)
      call colora("cyan")
      call tracec("E",xs(2), efs(2) , npsit-1,-1,-1,0._R8,0._R8)
!
      call colora("white")
      call frscj(1)
!
      call mapgsl(0._R8,1.0_R8,0.0001_R8,1.0_R8,.200_R8,.800_R8,.200_R8,  &  
     & .700_R8)
      call scalea(te(2),yplot(2),npsit-1,1._R8,temax,1._R8)
      do 330 n=2,npsit-1
  330 yplot(n) = max(yplot(n),0.0001_R8)
      call colora("magenta")
      call tracec(1hT,xs(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(ajtsc(2),yplot(2),npsit-1,1._R8,ajtscmax,1._R8)
      do 335 n=2,npsit-1
  335 yplot(n) = max(yplot(n),0.0001_R8)
      call colora("cyan")
      call tracec(1hJ,xs(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(ajpresf(2),yplot(2),npsit-1,1._R8,ajremax,1._R8)
      do 340 n=2,npsit-1
  340 yplot(n) = max(yplot(n)*ajremax/ajtscmax,0.0001_R8)
      call colora("green")
      call tracec(1hI,xs(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(etamidy(2),yplot(2),npsit-1,1._R8,etamax,1._R8)
      do 345 n=2,npsit-1
  345 yplot(n) = max(yplot(n),0.0001_R8)
      call colora("yellow")
      call tracec(1hR,xs(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call setch(10._R8,31._R8,1,2,0,-1)
!
      call colora("magenta")
      write(s100,9001)temax
      call gtextm(s100,80,0,1,1)
!
      call colora("green")
      write(s100,9002)
      call gtextm(s100,80,0,1,1)
!
      call colora("cyan")
      write(s100,9003) ajtscmax
      call gtextm(s100,80,0,1,1)
!
      call colora("yellow")
      write(s100,9004) etamax
      call gtextm(s100,80,0,1,1)
!
      call setch(10._R8,4.0_R8,1,2,0,-1)
      call colora("yellow")
      write(s100,9007) kcycle,times
      call gtextm(s100,80,0,1,1)
 9007 format(" cycle=",i6,"  times=",1pe12.4)
!
 9001 format("T...Electron temp scaled to ",1pe12.4)
 9002 format("I...RE current density &")
 9003 format("J...TSC current density scaled to ",1pe12.4)
 9004 format("R...Resistivity scaled to ",1pe12.4)
!
      call colora("white")
      call frscj(1)
!
                     endif
!
  195 return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
