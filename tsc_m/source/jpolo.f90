        subroutine  jpolo (ajphimax, polmax)
!*************************************************************************
!
!               Compute arrays VECX containing poloidal current in amperes
!                  VECZ containing temperature in eV at (i,j).
!                         ROS  2 Apr 1990.
!
!        21 APR 1993.  Compute R & Z forces on wire groups due to POLOIDAL
!               currents and write to file "wfpol"
!
!        21 APR 1993.  Compute R & Z forces on wire groups due to TOROIDAL
!               currents and write to file "wftor"
!
!        21 Apr 1993.  Revised to read wire numbers for force calculations
!               so that one version works for all machines.  Six
!               groups of wires are treated.
!        10 Apr 1994.  Vary THALO according to acoef(941-6).
!
      USE CLINAM
      USE SCRATCH
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ncall,nwf,nwt,nvv1,nvv2,nsiup1,nsiup2,nsilo1,nsilo2
      INTEGER nsoup1,nsoup2,nsolo1,nsolo2,nallw1,nallw2,ios22,i,j
      INTEGER ii,jj
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 polmax,ajphimax,dum2,dum3,tcur,tcurmain,tcurhalo,dip
      REAL*8 dipmain,diphalo,speedm,tempmax,g5,g6,g8,g9,ai,aj,speed
      REAL*8 tempval,ptestpl,tcurrat,plascold,timeold,tcurmainold
      REAL*8 tcurhalold,frvvtor,fzvvtor,frsiuptor,fzsiuptor
      REAL*8 frsilotor,fzsilotor,frsouptor,fzsouptor,frsolotor
      REAL*8 fzsolotor,frallwtor,fzallwtor,dpsir,dpsiz
!============
      logical  exforfil
      character*10   dum1(13)
      dimension  dum2(11), dum3(14)
!
      data  ncall, nwf, nwt /0, 97, 98/
!============      
!
      ncall = ncall + 1
!---------------
      if (ncall.eq.1)   then
!                 Set indices for wire groups
      nvv1   = acoef(951)
      nvv2   = acoef(952)
      nsiup1 = acoef(953)
      nsiup2 = acoef(954)
      nsilo1 = acoef(955)
      nsilo2 = acoef(956)
      nsoup1 = acoef(957)
      nsoup2 = acoef(958)
      nsolo1 = acoef(959)
      nsolo2 = acoef(960)
      nallw1 = acoef(961)
      nallw2 = acoef(962)
!
!           open file for writing plasma main & halo currents
!
      if (numargs .lt. 1 ) then
         filename = 'pcurdat'
      else 
         filename = 'pcurdat' // '.' // trim(suffix)
      end if
        inquire (file=trim(filename), exist=exforfil, err=9000)
        if (.not.exforfil)   then
         open (n13, file=trim(filename), status='new', iostat=ios22)
         write (nterm, 10)
   10       format('pcurdat',' non-existent: creating pcurdat.')
               call  tversion (n13)
                  write (n13,20)   name
   20             format (10a8)
               write (n13,30)
   30 format(/'      time        apl  tcuricalc       tcur   tcurmain    &  
     &tcurhalo    tcurrat      Thalo      dipdt    dipmain    diphalo    &  
     &   .....')
      else
!                         Append to pcurdat
        open (n13, file=trim(filename), status='old', err=9000)
      print *, ' Appending to pcurdat:'
   12    format (13a10)
   13    read (n13, 12, err=9500, end=16)  dum1
      write (nterm,12)   dum1
      go to 13
   16    continue
                 endif
!
!
!           open file for writing R & Z forces on wire groups
!           due to TOROIDAL currents
!
      if( numargs .lt. 1 ) then
         filename = 'wftor'
      else
         filename = 'wftor' // '.' // trim(suffix)
      end if
        inquire (file=trim(filename), exist=exforfil, err=9000)
        if (.not.exforfil)   then
         open (nwf, file=trim(filename), status='new', iostat=ios22)
         write (nterm, 33)
   33       format('wftortemp',' non-existent: creating new file.')
               call  tversion (nwf)
                  write (nwf,20)   name
               write (nwf,53)
   53   format (/' Forces on wire groups due to toroidal currents')
         write (nwf,31)  nvv1,nvv2, nsiup1,nsiup2, nsilo1,nsilo2,        &  
     &                  nsoup1,nsoup2, nsolo1,nsolo2, nallw1,nallw2
   31    format (/'    VV:  wires',i4,'-',i4                             &  
     &          /' SI_UP:  wires',i4,'-',i4                              &  
     &          /' SI_LO:  wires',i4,'-',i4                              &  
     &          /' SO_UP:  wires',i4,'-',i4                              &  
     &          /' SO_LO:  wires',i4,'-',i4                              &  
     &          /'  ALLW:  wires',i4,'-',i4)
            write (nwf,32)
   32 format(/'      time    fzvvtor  fzsiuptor fzsilotor fzsouptor      &  
     & fzso                                                              &  
     &lotor  fzallwtor    frvvtor  frsiuptor frsilotor frsouptor         &  
     & frsolot                                                           &  
     &or  frallwtor')
      else
!                       Append to wftor
        open (nwf, file=trim(filename), status='old', err=9000)
      print *, ' Appending to file wftor:'
   43    read (nwf, 12, err=9500, end=46)  dum1
      write (nterm, 12)   dum1
      go to 43
   46    continue
                 endif
!...
      if (acoef(904).eq.0.0_R8)  write (nout, 47)
      if (acoef(904).lt.0.0_R8)  write (nout, 48)
      if (acoef(904).gt.0.0_R8)      write (nout, 50)
      if (acoef(904).eq.0.0_R8)  write (nterm, 47)
      if (acoef(904).lt.0.0_R8)  write (nterm, 48)
      if (acoef(904).gt.0.0_R8)      write (nterm, 50)
   47    format (' JPOLO: JAYPH contains AJPHI, JPOLX, JPOLZ')
   48    format (' JPOLO: JAYPH contains AJPHI, JPOL, Temp in eV')
   50    format (' JPOLO: JAYPH contains AJPHI, BpolR, BpolZ')
!...
            endif
!---------------
      tcur = 0.0_R8
      tcurmain = 0.0_R8
      tcurhalo = 0.0_R8
      dip = 0.0_R8
      dipmain = 0.0_R8
      diphalo = 0.0_R8
!                       ******* 13 Jun 92
!
      speedm = 0
        tempmax = 0.0_R8
        ajphimax = 0.0_R8
!
        if (acoef(904).gt.0.0_R8)   go to 250
!
!                                       Start of loop
      do 200 i=iminn,imaxx
      do 100 j=jminn,jmaxx
      vecx(i,j) = 0._R8
      vecz(i,j) = 0._R8
      g5 = g(i,j)*xsqoj(i)
      g6 = g(i+1,j)*xsqoj(i+1)
      g8 = g(i,j+1)*xsqoj(i)
      g9 = g(i+1,j+1)*xsqoj(i+1)
      ai = .5_R8*(g9+g6-g8-g5)
      aj = .5_R8*(g8+g9-g5-g6)
      vecx(i,j) = - tpi*udsi*aj
      vecz(i,j) =   tpi*udsi*ai
      speed = vecx(i,j)**2 + vecz(i,j)**2
      if(speed.gt.speedm) speedm = speed
        if (ajphi(i,j).gt.ajphimax)   ajphimax = ajphi(i,j)
      ii = i
      jj = j
      call  getemp (ii, jj, tempval)
        if (tempval .gt. tempmax)   tempmax = tempval
      if (acoef(904).lt.0.0_R8)  then
                    vecx(i,j) = sqrt(speed)
                    vecz(i,j) = tempval
               endif
        if (iexv(i,j).eq.1)   go to 100
!                    point inside vessel
      ptestpl = phalo
      if (acoef(97).le.0.0_R8)   ptestpl = psilim
        if (psi(i,j).gt.ptestpl) go to 100
!                    point inside plasma+halo
      tcur = tcur + ajphi(i,j)
!                    point inside main plasma
      if (psi(i,j).lt.psilim)   tcurmain = tcurmain + ajphi(i,j)
!
  100 continue
  200 continue
!
        polmax = sqrt(speedm)
        ajphimax = ajphimax * udsi
      write (nout, 1004)   kcycle, times, ajphimax, polmax, tempmax
 1004 format (' JPOLO: kcycle, t, MAX ajphi, Ipol temp =' i7, 1p4e12.4)
      tcur     = tcur     * dxdz * udsi
      tcurmain = tcurmain * dxdz * udsi
      tcurhalo = tcur - tcurmain
      if (tcur.gt.1.E-4_R8)   tcurrat = tcurhalo/tcur
      if (ncall.gt.1)   then
         dip     = (apl - plascold)        /(times - timeold)
         dipmain = (tcurmain - tcurmainold)/(times - timeold)
         diphalo = (tcurhalo - tcurhalold) /(times - timeold)
               endif
      write (n13, 1006)  1000._R8*times, apl, tcurdtp*tpi*udsi,          &  
     &            tcur, tcurmain, tcurhalo, tcurrat,                     &  
     &            acoef(98), dip, dipmain, diphalo
 1006    format (0pf10.2, 1p13e11.3)
      plascold = apl
      tcurmainold = tcurmain
      tcurhalold  = tcurhalo
      timeold = times
!                       ******* 6 Jan 92
!ccccccccccccccc  if (idata.ne.3)      return
         call  sumwft (nvv1, nvv2, frvvtor, fzvvtor)
         call  sumwft (nsiup1, nsiup2, frsiuptor, fzsiuptor)
         call  sumwft (nsilo1, nsilo2, frsilotor, fzsilotor)
         call  sumwft (nsoup1, nsoup2, frsouptor, fzsouptor)
         call  sumwft (nsolo1, nsolo2, frsolotor, fzsolotor)
         call  sumwft (nallw1, nallw2, frallwtor, fzallwtor)
!ccccccc       call  sumwfp (nvv1, nvv2, frvvpol, fzvvpol)
!ccccccc       fzallwtot = fzallwtor
!ccccccc       fzallwtot = fzvvpol + fzallwtor
      write (nwf, 2006)  1000._R8*times, fzvvtor, fzsiuptor, fzsilotor,  &  
     & fzsouptor, fzsolotor, fzallwtor, frvvtor, frsiuptor, frsilotor,   &  
     & frsouptor, frsolotor, frallwtor
 2006    format (0pf10.2, 1p2e11.3, 3e10.3, 3e11.3, 3e10.3, 2e11.3)
!
!              Mod 10 Apr 1994.  Programmed thalo variation.
!     941   942   943   944   945   946
!     Time1    Time2    THALO1   THALO2   THALOMIN    THALOMAX
!
      if (acoef(944).eq.0.0_R8)  return
!
      if (times.lt.acoef(941) .or. times.gt.acoef(942))   return
      if (acoef(941) .ge. acoef(942))   return
      acoef(98) = acoef(943) + (acoef(944)-acoef(943))                   &  
     &         * (times-acoef(941)) / (acoef(942)-acoef(941))
      if (acoef(98).lt.acoef(945))  acoef(98) = acoef(945)
      if (acoef(98).gt.acoef(946))  acoef(98) = acoef(946)
!ccccccc    write (nterm, 240) 1000.*times,  acoef(98)
  240    format (' HALO TEMP: t,  th =', f10.3, 1p4e11.3)
        return
!..............................................................
  250   continue
!
!   acoef(904) .gt. 0       Put        Bpolr in VECX,          Bpolz in VECZ
!
      do 400 i=iminn,imaxx
      do 300 j=jminn,jmaxx
      vecx(i,j) = 0._R8
      vecz(i,j) = 0._R8
        dpsir = 0.2_R8* (psi(i+1,j) - psi(i-1,j))                        &  
     &        + 0.4_R8* (psi(i+2,j) - psi(i-2,j))
        if (j.lt.3)   then
                      dpsiz = psi(i,j+1) - psi(i,j-1)
        else
        dpsiz = 0.2_R8* (psi(i,j+1) - psi(i,j-1))                        &  
     &        + 0.4_R8* (psi(i,j+2) - psi(i,j-2))
                      endif
!              Bpolr and Bpolz
      vecx(i,j) = + 0.5_R8* dpsiz / (deez * xary(i))
      vecz(i,j) = - 0.5_R8* dpsir / (deex * xary(i))
        if (abs(vecx(i,j)).gt.abs(speedm))    speedm  = vecx(i,j)
        if (abs(vecz(i,j)).gt.abs(tempmax))   tempmax = vecz(i,j)
        if (ajphi(i,j).gt.ajphimax)   ajphimax = ajphi(i,j)
  300 continue
  400 continue
!
        ajphimax = ajphimax * udsi
        polmax = speedm
      write ( nout, 3003)   times, ajphimax, polmax, tempmax
      write (nterm, 3003)   times, ajphimax, polmax, tempmax
 3003    format (' JPOLO:  t, MAX ajphi, Br, Bz =' 1p4e12.4)
      return
 9000   write (nterm, 9010)
 9010    format (' Error in opening pcurdat or wftor')
      return
 9500   write (nterm, 9510)
 9510    format (' Error in reading pcurtemp or wftortemp')
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
