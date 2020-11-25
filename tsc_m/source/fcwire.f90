      subroutine fcwire
!......6.99 fcwire
!
!         calculate forces at wires (type 10 elements)
!
!         binary write R & z poloidal field (units = wb/m**2)
!                  and R & z forces (units = nt/rad)  at wires
!         into file wiref\suffix\
!
!         For a restart run, read from file wirer\suffix\
!                            write to  file wiref\suffix\
!
!         R. O. Sayer          14 June 88
!
!         binary write plasma toroidal current (units = amps/m**2) into
!         file jayph\suffix\
!
!         10 Feb 1988      5 point difference for psi derivatives
!                          New wire file format!! Includes r0, a, delt, eps
!                              xsep, zsep, xmag, zmag, + dumdum(1), ...,
!                              dumdum(20) for future output.
!
!         27 Nov 1989      modified to include poloidal current paths
!
!          2 Apr 1990      Format of JAYPH file modified to include poloidal
!                                       currents.
!          8 Jul 1993      polsegpl = current segments between plasma and wall
!
!          2 Sep 1993      Bpol R,Z in vecx,vecz  if acoef(904) = 55.
!
!         16 Jun 1994     vpolsegpl = poloidal voltages for polsegpl segments
!              are written to the plwall file.
!
      USE CLINAM
      USE FVV1
      USE SCR16
      USE SCR21
      USE SCR22
      USE SCR23
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!ccccccccc        common /scr23/  vpolsegpl(2*pncoil), vpolsegmax
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER namdum,iwopen,ibump,ierropn,ierrin,ierrout,luplwall,i
      INTEGER j,n,iw,jw,m,kowir,ii,iwdum,indxdum,ind2dum
      INTEGER jmindum,jmaxdum,iwrd,jwrd,iww
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dumdum,dumold,diamag,timdum,ajphimax,polmax,plcurdum
      REAL*8 timeold,plascold,xmagold,zmagold,r0old,aold,deltold
      REAL*8 epsold,xsepold,zsepold,dipdum,diadum,dpsir,dpsiz
!============
!       character*8  ijayph, iwirer, iwiref
        character*8  cijayph, ciwirer, ciwiref, ciplwall
!       equivalence (cijayph,ijayph), (ciwirer,iwirer), (ciwiref,iwiref)
        dimension  namdum(10), dumdum(20), dumold(20)
!       equivalence  (diamag, global(12))
!       equivalence  (vvl, dumdum(1))
        data   iwopen, ibump, ierropn, ierrin, ierrout / 0, 0, 0, 0, 0/
        data  cijayph,ciwirer,ciwiref /'jayph   ','wirer   ','wiref   '/    
        data  luplwall, ciplwall /48, 'plwall '/
!============      
!------------------------------------------------------------------------
!
        diamag = global(12)
        dumdum(1) = vvl
        if (iwopen.gt.0)   go to 80
!
!......initialize poloidal current array calculation
      do 109 i=1,nxp
      do 109 j=1,nzp
  109 itap1(i,j) = 0
!
      do 111 n=1,nwire
      iw = iwire(n)
      jw = jwire(n)
  111 itap1(iw,jw) = 1
!
      m = 0
      do 112 n=1,nwire
      iw = iwire(n)
      jw = jwire(n)
      if(itap1(iw+1,jw).eq.1) then
        m = m + 1
          xpc1(m) = xary(iw)
          zpc1(m) = zary(jw)
          xpc2(m) = xary(iw+1)
          zpc2(m) = zary(jw)
          orient(m) =  - tpi*udsi
      endif
!
      if(itap1(iw,jw+1).eq.1) then
        m = m+1
          xpc1(m) = xary(iw)
          zpc1(m) = zary(jw)
          xpc2(m) = xary(iw)
          zpc2(m) = zary(jw+1)
          orient(m) = + tpi*udsi
      endif
  112 continue
      mmax = m
      call  fcplwall (iwopen)
!------------------------------------------------------------------------
        if (irst1 .ne.1  .or. iwayne.le.1)   then
        if( numargs .lt. 1 ) then
           filename = trim(cijayph)
        else
           filename = trim(cijayph) // '.' // trim(suffix)
        end if
        open (lujayp, file=trim(filename), status='old', err=20)
        close(lujayp, status='delete')
   20   open (lujayp, file=trim(filename), status='new',                 &  
     &                form='unformatted',                                &  
     &                              err=104)
                           write(lujayp) nx,nz,alx,alz,isym,ccon
!
        if( numargs .lt. 1 ) then
           filename = trim(ciplwall)
        else
           filename = trim(ciplwall) // '.' // trim(suffix)
        end if
        open (luplwall, file=trim(filename), status='old', err=21)
        close(luplwall, status='delete')
   21   open (luplwall, file=trim(filename), status='new',               &  
     &                                       form='unformatted',         &  
     &                                       err=105)
        write (luplwall)  name
        write (luplwall)  nx,nz,alx,alz,isym,ccon, nwire, mplmax
        write (luplwall) (jplexv(m), kgrouplw(m), xplw1(m), zplw1(m),    &  
     &                    xplw2(m), zplw2(m), m=1,mplmax)
!
        if( numargs .lt. 1 ) then
           filename = trim(ciwiref)
        else
           filename = trim(ciwiref) // '.' // trim(suffix)
        end if
        open (luwire, file=trim(filename), status='old', err=22)
        close(luwire, status='delete')
   22   open (luwire, file=trim(filename), status='new',                 &  
     &                                     form='unformatted',           &  
     &                                     err=106)
                           write (nout, 62)   ciwiref, name
                           iwopen = 1
        write (luwire)  name
        write (luwire)  nx,nz,alx,alz,isym,ccon, nwire,mmax
        write (luwire)  (iwire(i),jwire(i),xwire(i),zwire(i), i=1,nwire)     
        write(luwire) (xpc1(m),zpc1(m),xpc2(m),zpc2(m),m=1,mmax)
                           dipdt = 0.0_R8
                           go to 80
!------------------------------------------------------------------------
        else
!
        write (nout, 32)   cijayph
   32   format (' Opening file : '  a8)
        if( numargs .lt. 1 ) then
           filename = trim(cijayph)
        else
           filename = trim(cijayph) // '.' // trim(suffix)
        endif
        open (lujayp, file=trim(filename), status='old',                 &  
     &                                     form='unformatted',           &  
     &                                     err=104)
        read (lujayp)   nx,nz,alx,alz,isym,ccon
        write (nout, 34)   nx,nz,alx,alz,isym,ccon
   34   format (' nx, nz, alx, alz, isym, ccon, nwire ='                 &  
     &            2i5, 2f8.3, i5, f8.3, 2i6)
!
        write (nout, 32)   ciplwall
        if( numargs .lt. 1 ) then
           filename = trim(ciplwall)
        else
           filename = trim(ciplwall) // '.' // trim(suffix)
        end if
        open (luplwall, file=trim(filename), status='old',               &  
     &                                       form='unformatted',         &  
     &                                       err=105)
        read (luplwall)  namdum
        read (luplwall)   nx,nz,alx,alz,isym,ccon, nwire, mplmax
        read (luplwall) (jplexv(m), kgrouplw(m), xplw1(m), zplw1(m),     &  
     &                   xplw2(m), zplw2(m), m=1,mplmax)
        write (nout, 34)  nx,nz,alx,alz,isym,ccon, nwire, mplmax
!
        kowir = iwayne / 20
        if (kowir.lt.1)   kowir = 1
        do 50  ii=1,iwayne-1
        read (lujayp, err=130, end=54) iwdum,timdum,                     &  
     &               indxdum,ind2dum,jmindum,jmaxdum, ajphimax, polmax
        if (iwdum.ne.ii)   go to 54
        iwrd = ind2dum - indxdum + 1
        jwrd = jmaxdum - jmindum + 1
        read (lujayp, err=130,end=54)  ((vecx(i,j),i=1,iwrd),j=1,jwrd)
        read (lujayp, err=130,end=54)  ((vecx(i,j),i=1,iwrd),j=1,jwrd)
        read (lujayp, err=130,end=54)  ((vecx(i,j),i=1,iwrd),j=1,jwrd)
!
        if (ii.ne.1.and.ii.lt.iwayne-2.and.mod(ii,kowir).ne.0) go to 50
        write (nout,  45)   iwdum, timdum, ajphimax, polmax
        write (nterm, 45)   iwdum, timdum, ajphimax, polmax
   45   format (' JAYPHI set:' i4,'    t, ajphimax, polmax ='1p3e10.3)
   50   continue
!
        go to 55
   54   backspace  lujayp
!----
   55    continue
        kowir = iwayne / 20
        if (kowir.lt.1)   kowir = 1
        do 56  ii=1,iwayne-1
        read (luplwall, err=135, end=58) iwdum,timdum, plcurdum,         &  
     &                                   plsegmax, kgrplmax, vpolsegmax
        read (luplwall, err=135, end=58) (polsegpl(m),  m=1,mplmax)
        read (luplwall, err=135, end=58) (vpolsegpl(m), m=1,mplmax)
        if (ii.ne.1.and.ii.lt.iwayne-2.and.mod(ii,kowir).ne.0) go to 56
        write (nout,  57)   iwdum, timdum, plcurdum,                     &  
     &                      plsegmax, vpolsegmax, kgrplmax
        write (nterm, 57)   iwdum, timdum, plcurdum,                     &  
     &                      plsegmax, vpolsegmax, kgrplmax
   56   continue
   57 format(' PLWALL set:' i4,'    t, Ip, Max I & V, grp ='1p4e10.3,i7)     
!
        go to 60
   58   backspace  luplwall
!----
   60   ijwname = ciwirer
        write (nterm, 32)   ciwirer
        write (nout, 32)   ciwirer
        if( numargs .lt. 1 ) then
          filename = trim(ciwirer)
        else
          filename = trim(ciwirer) // '.' // trim(suffix)
        end if
        open (luwirer, file=trim(filename), status='old',                &  
     &                                      form='unformatted',          &  
     &                                      err=106)
        iwopen = 1
        read (luwirer, err=130)  namdum
        write (nout, 61)   namdum
   61   format (1x, 10a8)
!
        if( numargs .lt. 1 ) then
           filename = trim(ciwiref)
        else
           filename = trim(ciwiref) // '.' // trim(suffix)
        endif
        open (luwire, file=trim(filename), status='old', err=961)
        close (luwire, status='delete')
  961   open (luwire, file=trim(filename), status='new',                 &  
     &                                     form='unformatted',           &  
     &                                     err=106)
        write (luwire, err=150)   name
        write (nout, 62)   ciwiref, name
   62   format (' Creating new wire file: ' a8/ 1x, 10a8)
!
        read (luwirer, err=130)  nx,nz,alx,alz,isym,ccon, nwire,mmax
        write(luwire,  err=150)  nx,nz,alx,alz,isym,ccon, nwire,mmax
        write (nout, 34)   nx,nz,alx,alz,isym,ccon, nwire,mmax
!
        read (luwirer,err=130)  (iwire(i),jwire(i),xwire(i),zwire(i),    &  
     &                          i=1,nwire)
        write(luwire, err=150)  (iwire(i),jwire(i),xwire(i),zwire(i),    &  
     &                          i=1,nwire)
        read(luwirer,err=130) (xpc1(m),zpc1(m),xpc2(m),zpc2(m),m=1,mmax)     
        write(luwire,err=150) (xpc1(m),zpc1(m),xpc2(m),zpc2(m),m=1,mmax)     
!
!
        write (nout, 63)
   63  format(/' Set      time    pl cur    xmag    zmag      Ro         &  
     & a                                                                 &  
     &    delt     eps    xsep    zsep    dIp/dt    diamag    cvvpol')
!......
        do 65  ii=1,iwayne-1
        read (luwirer,err=130, end=70)   iww, timeold, plascold,         &  
     &                xmagold,zmagold, r0old, aold, deltold,             &  
     &                epsold, xsepold, zsepold, dipdum, diadum, dumold
        write(luwire, err=150)           iww, timeold, plascold,         &  
     &                xmagold,zmagold, r0old, aold, deltold,             &  
     &                epsold, xsepold, zsepold, dipdum, diadum, dumold
!
        read (luwirer,err=130, end=70)   (cwires(i), bpolr(i),bpolz(i),  &  
     &                                  fwirer(i), fwirez(i), i=1,nwire)    
        write(luwire, err=150)           (cwires(i), bpolr(i),bpolz(i),  &  
     &                                  fwirer(i), fwirez(i), i=1,nwire)     
        read(luwirer,err=130, end=70) (polsegc(m),m=1,mmax)
        write(luwire, err=150)        (polsegc(m),m=1,mmax)
!
        if (ii.ne.1.and.ii.lt.iwayne-2.and.mod(ii,kowir).ne.0) go to 65
!
        write (nout, 64) iww, timeold, plascold,                         &  
     &                xmagold,zmagold, r0old, aold, deltold,             &  
     &                epsold, xsepold, zsepold, dipdum, diadum,          &  
     &                dumold(3)
   64   format (i4, f10.6, e10.3, 8f8.3, 3e10.2)
        write (nterm,100) iww, 1000._R8*timeold,1.E-6_R8*plascold,       &
     &  1.E-9_R8*dipdum
   65   continue
                           endif
   70   close (luwirer)
      
!------------------------------------------------------------------------
   80   if (iwayne.lt.1)   return
!
!.....calculate poloidal current array polseg
      do 210 i=1,nxp
      do 210 j=1,nzp
  210 itap1(i,j) = 0
      do 211 n=1,nwire
      iw = iwire(n)
      jw = jwire(n)
  211 itap1(iw,jw) = 1
!
      m = 0
      do 212 n=1,nwire
      iw = iwire(n)
      jw = jwire(n)
      if(itap1(iw+1,jw) .eq. 1) then
          m = m+1
          polsegc(m) =(-g(iw+1,jw)*xsqoj(iw+1)+g(iw+1,jw+1)*xsqoj(iw+1))  &  
!    &                                                                   &  
     &               *orient(m)
          endif
      if(itap1(iw,jw+1) .eq. 1) then
          m = m+1
          polsegc(m) = ( +g(iw+1,jw+1)*xsqoj(iw+1)-g(iw,jw+1)*xsqoj(iw))  &  
!    &                                                                   &  
     &                *orient(m)
          endif
!
  212 continue


      mmax = m
!           Calculate poloidal current between plasma & wall
      call  fcplwall (iwopen)
!
        call  ivvpol
        call  jpolo (ajphimax, polmax)
        write(lujayp) iwayne,times,iminn,imaxx,jminn,jmaxx,              &  
     &                             ajphimax, polmax
        write (lujayp)  ((ajphi(i,j)*udsi,i=iminn,imaxx),j=jminn,jmaxx)
        write (lujayp)        ((vecx(i,j),i=iminn,imaxx),j=jminn,jmaxx)
        write (lujayp)        ((vecz(i,j),i=iminn,imaxx),j=jminn,jmaxx)
!
        write (luplwall) iwayne,times, apl, plsegmax, kgrplmax,          &  
     &                   vpolsegmax
        write (luplwall) (polsegpl(m),  m=1,mplmax)
        write (luplwall) (vpolsegpl(m), m=1,mplmax)
!
!------------------------------------------------------------------------
        do 30  ii=1,nwire
        n = ncoil - nwire + ii
        cwires(ii) = ccoil(n) * udsi
        i = iwire(ii)
        j = jwire(ii)
        dpsir = 0.2_R8* (psi(i+1,j) - psi(i-1,j))                        &  
     &        + 0.4_R8* (psi(i+2,j) - psi(i-2,j))
!
        if (j.lt.3)   then
                      dpsiz = psi(i,j+1) - psi(i,j-1)
        else
        dpsiz = 0.2_R8* (psi(i,j+1) - psi(i,j-1))                        &  
     &        + 0.4_R8* (psi(i,j+2) - psi(i,j-2))
                      endif
        bpolr(ii) = + 0.5_R8* dpsiz / (deez * xwire(ii))
        bpolz(ii) = - 0.5_R8* dpsir / (deex * xwire(ii))
        fwirer(ii) = -0.5_R8* cwires(ii) * dpsir / deex
        fwirez(ii) = -0.5_R8* cwires(ii) * dpsiz / deez
   30   continue
!
!ccccccc        call  fcwset
!
        plascur = apl
!
        if (iwayne.gt.1)   dipdt = (plascur-plascold) / (times-timeold)
        timeold = times
        plascold = plascur
        write (luwire)   iwayne,times, plascur, xmag ,zmag,              &  
     &                   rmajor,rminor,shape3, shape5,shape6,shape7,     &  
     &                   dipdt, diamag, dumdum
        write (luwire)   (cwires(i), bpolr(i),bpolz(i),                  &  
     &            fwirer(i), fwirez(i), i=1,nwire)
        write (luwire)   (polsegc(m),m=1,mmax)
!
        if (mod(iwayne,5).eq.0)   write (nout, 63)
        write (nout, 64) iwayne,times, plascur, xmag ,zmag,              &  
     &                   rmajor,rminor,shape3, shape5,shape6,shape7,     &  
     &                   dipdt, diamag, cvvpol
        if (mod(iwayne,5).eq.0)   write (nterm, 99)
       write (nterm,100)   iwayne, 1000._R8*times, 1.E-6_R8*plascur,     &  
     &       1.E-9_R8*dipdt, 1.E-6_R8*plsegmax, kgrplmax, acoef(98)
   99   format ('WIREF: Set       Time        Ip    dIp/dt  plsegmax     &  
     &THALO')
  100   format ('WIREF:', i4, f11.3, 3f10.3, i10, f8.2)
        iwayne = iwayne + 1
        tjphi  = tjphi  + dtjphi
      acoef(30) = iwayne
      acoef(31) = tjphi
        if (ibump.gt.20)   return
        if (acoef(31).ge.0.0_R8)   return
        ibump = ibump + 1
        tjphi = tjphi - 0.5_R8*dtjphi
      acoef(31) = tjphi
        return
!------------------------------------------------------------------------
  104   ierropn = ierropn + 1
      if (ierropn.gt.20)   stop
      write (nout, 120)   cijayph
        write (nterm,120)   cijayph
        return
  105   ierropn = ierropn + 1
      if (ierropn.gt.20)   stop
      write (nout, 120)   ciplwall
        write (nterm,120)   ciplwall
        return
  106   ierropn = ierropn + 1
      if (ierropn.gt.20)   stop
      write (nout, 120)   ciwiref
        write (nterm,120)   ciwiref
        return
  110   ierropn = ierropn + 1
      if (ierropn.gt.20)   stop
      write (nout, 120)   ijwname
        write (nterm,120)   ijwname
  120   format (' Error opening file :' a8)
        return
  130   ierrin = ierrin + 1
      if (ierrin.gt.5)   stop
      write (nout, 140)   ijwname, ierrin
        write (nterm,140)   ijwname, ierrin
  140   format (' Error reading data from file :' a8, i10)
        return
  135   ierrin = ierrin + 1
      if (ierrin.gt.5)   stop
      write (nout, 140)   ciplwall, ierrin
        write (nterm,140)   ciplwall, ierrin
        return
  150   ierrout = ierrout + 1
      if (ierrout.gt.5)   stop
      write (nout, 160)   ciwiref, ierrout
        write (nterm,160)   ciwiref, ierrout
  160   format (' Error writing data to file :' a8, i10)
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
