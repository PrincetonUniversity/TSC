
!#include "f77_dcomplx.h"
      subroutine savprof
!**********************************************************************
!
!***********************************************************
!                                                          *
!....... saves data for profile plots in "plscr1"          *
!
!                                                          *
!.....scary(n,j) = midplane profile where :                *
!.....=====================================                *
!.....n =  1 ... x-coordinate                              *
!          2 ... toroidal field                            *
!          3 ... poloidal field                            *
!          4 ... pressure                               *
!          5 ... vertical field                            *
!          6 ... d(radial field)/dz                        *
!          7 ... field index                               *
!          8 ... toroidal current                          *
!          9 ... loop-volt                          *
!         10 ... poloidal flux                             *
!.......  11 ... safety factor                             *
!.......  12 ... electron temperature                      *
!.......  13 ... ion temp                   *
!         14 ... electron density
!         15 ... lamda  profile
!         16 ... resistivity
!         16+n ... divertor plate n
!
!                                                          *
!***********************************************************
!
      USE CLINAM
      USE SAPROP
      USE SCR1
      USE EQRUNS
      USE PROFCOM
      USE NEWPLOT

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,nrecor2,iadscr,ndivp,mwpl1
      INTEGER idisk,ios18,imin,imax,np2s,min,n,nhp,jmax,j,l,npmax
      INTEGER ii,irec,nptsx,jx,ihdb,js,jsm,ihd,lstop,lcolmn2,lcolmn3
      INTEGER len, imin1, imax1 
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 degx,degz,degy,frad,fwr,ss,ans,voltl,ps,pval
      REAL*8 ppval,eval,epval,rval,rpval,tival,rojmin,denval,bmod
      REAL*8 sinxa,cosxa,sinza,cosza,sinya,cosya,ymin,ymax
      REAL*8 xmin,xmax,yval,width,depth,height,xmapl,xmapr,xco,yco
      REAL*8 zco,xndivp,delx,dely,zcl,xcl,yloc,ycl,xcb,ycb,xc,yc
      REAL*8 xloc,xpltj,xpltjm,frac,xs,ys,xe,ye,yf,ylabset,xlabset
      REAL*8 plcprnt,alamda,scar5,scar6
      REAL*8 sum, AREAL
!============
!     common/eqruns/findex
!
!     integer pnplts,pnplt
!     parameter (pnplt=16,pnplts=pnplt+pnplat)
!
!     common/profcom/
!    .          xplt(pnx+ppsi+pnseg),yplt(pnx+ppsi+pnseg),
!    1          scary(pnplts+10,pnx+ppsi+pnseg),
!    .          psix(pnx),timesv(50),
!    .          ymaxsf(500),yminsf(500),xcs(4),
!    .          ycs(4),div(2),np2sav,npsitsv(50),ihds(4)
!
!     character*8 lab(pnplts+10)
      character*8 lab(100)
      character*1 num(40)
      character*8 cdiv(2)
!
      data num / "*","1","2","3","4","5","6","7","8","9",                &  
     &           "a","b","c","d","e","f","g","h","i","j",                &  
     &           "k","l","m","n","o","p","q","r","s","t",                &  
     &           "u","v","x","w","y","z","!","@","+","%"/
      data (lab(i),i=1,25)/'x','torfield','polfield','pres-kpa',         &  
     &             'bz      ','d(bx)/dz','index   ',                     &  
     &             'current ','loopvolt',                                &  
     &             'pol flx ','q-prof  ','te-prof ',                     &  
     &             'ti-prof ','ne-prof ','d-source',                     &  
     &             'resistiv','hflux-p1','hflux-p2',                     &
     &             'hflux-p3','hflux-p4','hflux-p5',                     &  
     &             'hflux-p6','hflux-p7','hflux-p8',                     &  
     &             'hflux-p9'/
!
      data nrecor2,iadscr /0,0/
      data degx / 0.0_R8/
      data degz / 15.0_R8/
      data degy / 90.0_R8/
      data frad / 0.8_R8/
      data fwr / 5._R8/
      data ndivp / 500 /
!
!...........................................................
!
!.....open scratch file
!
!...........................................................
      if(ifrst(6).ne.0) go to 50
      ss = 0
      mwpl1 = 15
      idisk = 1
!     ifiles3(1:6) = 'plscr1'
!     ifiles3(7:7) = isuffix(1:1)
      if( numargs .lt. 1 ) then
         ifiles3 = 'plscr1' // isuffix(1:1)
      else 
         ifiles3 = 'plscr1' // '.' // trim(suffix)
      end if
      len = ((ncycle-kcycle)/nskipr+2)*pnx*pnplts
      open(nsc3,file=trim(ifiles3),status='unknown',form='unformatted',  &  
     &     iostat=ios18)
      ifrst(6) = 1
   50 continue
      call xlimits(imin,imax)
      np2s = imax-imin+1
      np2sav = np2s
      imin1 = imin-1
      imax1 = imax+1
!...........................................................
!
!.....define external flux array on midplane
!
!...........................................................
      do 75 i=imin1,imax1
!
!.....upper half plane
!
      sum = 0
      do 78 n=1,nwire
      ans = 0._R8
      if(xary(i).eq.xary(iwire(n)).and.zary(jmag).eq.zary(jwire(n)))     &  
     &    go to 78
      call gf(ineg,nmult,xary(i),zary(jmag),xary(iwire(n))               &  
     &       ,zary(jwire(n)),ans)
   78 sum = sum + ans*ccoil(ncoil-nwire+n)/tpi
      if(ncoil.eq.nwire) go to 80
      do 79 n=1,ncoil-nwire
      call gf(ineg,nmult,xary(i),zary(jmag),xcoil(n),zcoil(n),ans)
   79 sum = sum + ans*ccoil(n)/tpi
   80 continue
      if(isym.eq.0) go to 180
!
!.....lower half plane
!
      do 178 n=1,nwire
      call gf(ineg,nmult,xary(i),zary(jmag),xary(iwire(n))               &  
     &       ,-zary(jwire(n)),ans)
      if(zary(jwire(n)).eq.0) ans = 0
  178 sum = sum + ans*ccoil(ncoil-nwire+n)/tpi
      if(ncoil.eq.nwire) go to 180
      do 179 n=1,ncoil-nwire
      call gf(ineg,nmult,xary(i),zary(jmag),xcoil(n),-zcoil(n),ans)
      if(zcoil(n).eq.0) ans = 0._R8
  179 sum = sum + ans*ccoil(n)/tpi
  180 continue
      psix(i) = sum
   75 continue
!...........................................................
!
!.....fill up the "scary" arrays
!
!...........................................................
      nhp=jmag+1
      jmax = np2s+1
      j = 0
      do 100 i=imin,imax
      j = j+1
      scary(1,j) = xary(i)
      scary(2,j) = 0.5_R8*(g(i,nhp)*xsqoj(i)+g(i+1,nhp)*xsqoj(i+1))      &  
     &            /xary(i)
      scary(3,j) =- (psi(i+1,jmag)-psi(i-1,jmag))/(xary(i)*(xary(i+1)    &  
     &            -xary(i-1)))
      scary(4,j) = .5_R8*(pr(i,nhp)                                      &  
     &                  +  pr(i+1,nhp) )*udsp
      scary(5,j) = -(psix(i+1)-psix(i-1))/(2._R8*deex*xary(i))
      scary(6,j) = -((psix(i+1)-2._R8*psix(i)+psix(i-1))/deex**2         &  
     &           + scary(5,j))/xary(i)
      scary(7,j) = -xary(i)*scary(6,j)/(scary(5,j)+1.0E-33_R8)
      scary(7,j) = max(-5.0_R8,min(scary(7,j),5.0_R8))
      scary(8,j) = ajphi(i,jmag)
      scary(9,j) = voltl(i,jmag)
      if(isurf.eq.1) then
      ps = psi(i,jmag)
      call peval(ps,2,pval,ppval,i,jmag)
      call eeval(ps,2,eval,epval,i,jmag)
      call reval(ps,idens,1,rval,rpval,i,jmag)
      call tieval(ps,2,tival,i,jmag,pval,eval,rval)
      scary(14,j) = rval*udsd
      scary(12,j) = (eval*udsh)/(rval*udsd)
      scary(13,j) = tival
      else
      scary(14,j) = .5_R8*(roj(i,nhp)                                    &  
     &                  +  roj(i+1,nhp) )*udsd
      rojmin = fracn0*r0
      if(scary(14,j) .lt. rojmin*udsd) scary(14,j) = rojmin*udsd
      scary(12,j) = .25_R8*(pr(i,nhp)                                    &  
     &                  +  pr(i+1,nhp) )*udsh/scary(14,j)
      scary(13,j) = .25_R8*(pr(i,nhp)                                    &  
     &                  +  pr(i+1,nhp) )*udsh/scary(14,j)
      endif
!
!.....October ,1996  SCJ
      if((idens.ne.0 .and. idens.ne.3)) go to 104
!.....special time integrated density source plot
      call deneval(ps,denval)
      scary(15,j) = denval
 104  continue
      scary(16,j) = etay(i,jmag)*udsr
  100 continue
!
!...more frequent profile output for DIII-D experiment matching
!
!     write(nout,1069) kcycle,times,dts
!     do 900 j=1,npsit
! 900 write(nout,2069) j,diffnema(j),diffary(j),ane(j), chiitima(j),     &
!    &  chiinca(j),chiicopi(j), ti(j), chietema(j), chienca(j),          &
!    &  chiecopi(j), te(j)  
!1069 format(1h1," cycle",i7,"    time",1pe12.4,"     dt",1pe12.4,//,4x, &
!    & "j   diffnema    D ~ 1/n    ne(mks)  chiitima   chiinca  ",       &
!    & "  chi ~ 1/n    ti(ev)   chietema    chienca   chi ~ 1/n",        &
!    & "     te" )
!2069 format(i5,1p11e11.3)
!
!...end of more frequent profile output for DIII-D experiment matching
!
!
!
!.....special output needed for ET
      write(nout,1117) kcycle, times
 1117 format(" special R. Taylor output: cycle=",i7,                     &  
     &       "   time = ",1pe12.4,"(s)",/,                               &  
     &" I     R        P(Kpa)      Te(keV)    Ti(kev)     ne(mks)",      &  
     &"  J(A/m**2)  Pol flx(W)    BMOD")
      do j=1,jmax
      bmod = sqrt(scary(2,j)**2 + scary(3,j)**2)
      write(nout,1118)j,scary(1,j),scary(4,j)*1.E-3_R8,scary(12,j)*      &
     &          1.E-3_R8,                                                &  
     &          scary(13,j)*1.E-3_R8, scary(14,j),scary(8,j),            &  
     &          psi(imin+j-1,jmag)*tpi, bmod
 1118 format(i3,1p8e11.3)
      enddo
!.....end of special output needed for ET
!
      if(isurf.ne.0) go to 110
      if(lrswtch.eq.0) call qcalc
      if(lrswtch .gt. 0) go to 120
      j = 0
      do 105 i=1,npts
      j = j+1
      scary(10,j) = .5_R8*(polflx(j)+polflx(j+1))-polflx(1)
      scary(11,j) = qprof(j)
  105 continue
      go to 119
  110 continue
      do 106 j=1,npsi
      scary(10,j) = xsv2(j) - xsv2(1)
      scary(11,j) = qprof2(j)
  106 continue
  119 continue
      if(iplate.eq.0 .or. nplate.eq.0) go to 120
      do 118 n=1,nplate
      do 117 l=1,nseg(n)
      scary(pnplt+n,l) = hplate(n,l)
  117 continue
  118 continue
!
  120 continue
!...........................................................
!
!.....write onto scratch disk plscr1
!
!...........................................................
      if(nrecor2 .gt. 40) return
      nrecor2 = nrecor2 + 1
      timesv(nrecor2) = time*udst
      npsitsv(nrecor2) = npsit
      call bufout(nsc3,scary(1,1),scary(pnplts,pnx+ppsi+pnseg))
!
      return
!
      entry plotit1
!
!***********************************************************
!                                                          *
!....... plots data for profile plots                      *
!                                                          *
!***********************************************************
!
      if(ifrst(6).eq.0) return
!...........................................................
!
!.....start loop on making profile plots
!
!...........................................................
      sinxa = sin(degx*pi/180._R8)
      cosxa = cos(degx*pi/180._R8)
      sinza = sin(degz*pi/180._R8)
      cosza = cos(degz*pi/180._R8)
      sinya = sin(degy*pi/180._R8)
      cosya = cos(degy*pi/180._R8)
!
      npmax = pnplts
      if(lrswtch .gt. 0) npmax = 9
      do 400 ii=2,npmax
      if(ii.eq.2.and.noplot(80).gt.0._R8) goto 400
      if(ii.eq.3.and.noplot(81).gt.0._R8) goto 400
      if(ii.eq.5.and.noplot(82).gt.0._R8) goto 400
      if((ii.eq.6._R8.or.ii.eq.7).and.noplot(83).gt.0._R8) goto 400
      if(ii.eq.9.and.noplot(84).gt.0._R8) goto 400
      if(ii.eq.12.and.noplot(85).gt.0._R8) goto 400
      if(ii.eq.13.and.noplot(86).gt.0._R8) goto 400
      if(ii.eq.14.and.noplot(87).gt.0._R8) goto 400
      if(ii.eq.10) go to 400
!     write(nterm,8881) ii,npmax
!8881 format("ii,npmax =",2i5)
      if(ii.eq.15 .and. (idens.ne.0 .and. idens.ne.3)) go to 400
      if(ii.gt.pnplt+nplate) go to 400
      if(ii.eq.11.and.noplot(90).gt.0._R8) goto 400
!
!.....first determine ymin and ymax
!
      ymin = 1.E40_R8
      ymax =-1.E40_R8
      xmin = 1.E40_R8
      xmax =-1.E40_R8
      rewind nsc3
!
      if(nrecor2.gt.40) then
      write(nsc1,1090) nrecor2
      nrecor2 = 40
 1090 format(' more than 40 profiles requested, nrecor2=',i4)
      call frscj(6)
                        endif
      do 450 irec=1,nrecor2
      np2s = np2sav
      if(ii.eq.11 .and. isurf.ne.0) np2s = npsitsv(irec)
      if(ii.eq.11 .and. isurf.eq.0) np2s=npts
      if(ii.gt.pnplt) np2s = nseg(ii-pnplt)
!
      call bufin(nsc3,scary(1,1),scary(pnplts,pnx+ppsi+pnseg))
!
      do 460 j=1,np2s
      yval = scary(ii,j)
! pressure units in kN/m^2
      if(ii.eq.4) yval=scary(ii,j)*1.E-3_R8
! toroidal current density units: kA/m^2
      if(ii.eq.8) yval=scary(ii,j)*udsi*1.E-3_R8
      if(yval.gt.ymax) ymax = yval
      if(yval.lt.ymin) ymin = yval
  460 continue
      if(ii.gt.pnplt) go to 454
      if(ii.eq.11 ) go to 455
      if(scary(1,np2s).gt.xmax) xmax = scary(1,np2s)
      if(scary(1,1).lt.xmin) xmin = scary(1,1)
      go to 450
  455 continue
      xmin = 0._R8
      if(scary(10,np2s) .gt. xmax) xmax = scary(10,np2s)
      go to 450
  454 continue
      if(noplot(9) .gt. 0) go to 400
      xmin = dplate(ii-pnplt,1)
      xmax = dplate(ii-pnplt,nseg(ii-pnplt))
  450 continue
      if(ymax .le. ymin + 1.E-8_R8*abs(ymin)) go to 400
!
!.....now plot the profiles
!
      degz = 15.0_R8
!     if(ii.eq.11) degz = 165.
      cosza = cos(degz*pi/180._R8)
      degx = 0.0_R8
      if(ii.eq.11) degx = 180._R8
      cosxa = cos(degx*pi/180._R8)
      width = fwr*(xmax-xmin)
      depth = 0.8_R8*(width-(xmax-xmin))/abs(cosza )
      height = (sinza*depth + (xmax-xmin))/0.8_R8
      xmapl = 0.0_R8
      if(ii.eq.11) xmapl = -(xmax-xmin)
      xmapr = xmapl + width
      xco = 0.1_R8*width + xmapl
      yco = 0.1_R8*height + (xmax-xmin)
      zco = 0.0_R8
      xndivp = ndivp
      delx = width/xndivp
      dely = ymax - ymin
      do 3015 j = 1,ndivp
      ymaxsf(j) = -1.E30_R8
      yminsf(j) = 1.E30_R8
 3015 continue
!
      if(acoef(901).eq.0) then
        select case(ii)
          case(7)
          call mapg(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.290_R8,.6_R8)
!
          case(16)
          call mapgsl(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.685_R8,1._R8)
!
          case default
          call mapg(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.685_R8,1._R8)
        end select
      else
        select case(ii)
          case(4)
          call mapg(xmin,xmax,ymin,ymax,.1_R8,.45_R8,.55_R8,.9_R8)
!
          case(7)
          call mapg(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.290_R8,.6_R8)
!
          case(8)
          call mapg(xmin,xmax,ymin,ymax,.6_R8,.95_R8,.55_R8,.9_R8)
!
          case(16)
          call mapgsl(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.685_R8,1._R8)
          case default
            call mapg(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.685_R8,1._R8)
        end select
      endif
!
      rewind nsc3
      do 470 irec=1,nrecor2
      np2s = np2sav
      if(ii.eq.11 .and. isurf.ne.0) np2s = npsitsv(irec)
      if(ii.eq.11 .and. isurf.eq.0) np2s=npts
      if(ii.gt.pnplt) np2s = nseg(ii-pnplt)
!
      call bufin(nsc3,scary(1,1),scary(pnplts,pnx+ppsi+pnseg))
!
      do 480 j=1,np2s
      yplt(j) = scary(ii,j)
! pressure units in kN/m^2
      if(ii.eq.4) yplt(j)=scary(ii,j)*1.E-3_R8
! toroidal current density units: kA/m^2
      if(ii.eq.8) yplt(j)=scary(ii,j)*udsi*1.E-3_R8
      xplt(j) = scary(1,j)
      if(ii.eq.11 ) xplt(j) = scary(10,j)
      if(ii.gt.pnplt) xplt(j) = dplate(ii-pnplt,j)
  480 continue
      if(acoef(901).eq.0) then
        select case(ii)
          case(7)
          call map(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.290_R8,.6_R8)
!
          case(16)
          call mapsl(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.685_R8,1._R8)
!
          case default
          call map(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.685_R8,1._R8)
        end select
      else
        select case(ii)
          case(4)
          call map(xmin,xmax,ymin,ymax,.1_R8,.45_R8,.55_R8,.9_R8)
!
          case(7)
          call map(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.290_R8,.6_R8)
!
          case(8)
          call map(xmin,xmax,ymin,ymax,.6_R8,.95_R8,.55_R8,.9_R8)
!
          case(16)
          call mapsl(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.685_R8,1._R8)
          case default
            call map(xmin,xmax,ymin,ymax,.155_R8,.550_R8,.685_R8,1._R8)
        end select
      endif

      if(acoef(901).eq.0._R8)                                            &  
     &   call tracec(num(irec),xplt,yplt,np2s,-1,-1,0._R8,0._R8)
      if(acoef(901).gt.0._R8)                                            &  
     &   call trace(xplt,yplt,np2s,-1,-1,0._R8,0._R8)
!
      if(ii.eq.6 .or. ii.eq.7) go to 470
      if(nrecor2 .le. 1 .or. timesv(nrecor2).eq.timesv(1)) go to 500
      call map (xmapl,xmapr,0.0_R8,height,.05_R8,.95_R8,.25_R8,.75_R8)
!
      zcl = zco + depth*(timesv(irec)-timesv(1))                         &  
     &                 /(timesv(nrecor2)-timesv(1))
!
      do 520 n = 1,3
      nn = 2
      nptsx = (xplt(np2s)-xplt(1))/delx + 1
      if(n.gt.2) go to 525
      nn = 1
      if(n.eq.2) nn = 3
      xcs(nn) = xcs(nn+1)
      ycs(nn) = ycs(nn+1)
      ihds(nn) = ihds(nn+1)
      j =  1
      if(n.eq.2) j = np2s
      xcl = (xplt(j)-xmin) + xco
      yloc = (yplt(j)-ymin)/dely
      ycl = yco - (xmax-xmin)*min(1._R8,1._R8-frad*yloc)
      xcs(nn+1) = xcl*cosxa + zcl*cosza + ycl*cosya
      ycs(nn+1) = xcl*sinxa + zcl*sinza + ycl*sinya
      nptsx = (xcs(nn+1)-xcs(nn))/delx + 1
      ihds(nn+1) = 1
      jx = min(xndivp,xndivp*(xcs(nn+1)-xmapl)/(xmapr-xmapl))
      if(ycs(nn+1).lt.ymaxsf(jx)                                         &  
     &           .and. ycs(nn+1).gt.yminsf(jx)) go to  525
      if(ycs(nn+1).ge.ymaxsf(jx)) go to 526
      ihds(nn+1) = -1
      yminsf(jx) = ycs(nn+1)
      if(ycs(nn+1).ge.ymaxsf(jx)) ymaxsf(jx) = ycs(nn+1)
      go to 525
  526 ihds(nn+1) = 0
      ymaxsf(jx) = ycs(nn+1)
      if(ycs(nn+1).le.yminsf(jx)) yminsf(jx) = ycs(nn+1)
  525 ihdb = ihds(nn)
      xcb = xcs(nn)
      ycb = ycs(nn)
      js = 2
      if(irec.eq.1 .and. n.lt.3) go to 520
      do 530 i = 2,nptsx
      go to(527,527,528),n
  527 xc = (xcs(nn+1)-xcs(nn))*AREAL(i-1)/AREAL(nptsx-1) + xcs(nn)
      yc = (ycs(nn+1)-ycs(nn))*AREAL(i-1)/AREAL(nptsx-1) + ycs(nn)
      go to 545
  528 xloc = (xplt(np2s)-xplt(1))*AREAL(i-1)/AREAL(nptsx-1) + xplt(1)
      xcl = (xloc-xmin) + xco
      do 535 j = js,np2s
      xpltj = (xplt(j)-xmin) + xco
      if(xcl.gt.xpltj) go to  535
      js = j
      go to 540
  535 continue
      js = np2s
  540 jsm = js-1
      xpltjm = (xplt(jsm)-xmin) + xco
      frac = max(0._R8,(xcl-xpltjm)/(xpltj-xpltjm))
      frac = min(1._R8,frac)
      yloc = (frac*(yplt(js)-yplt(jsm))+yplt(jsm)-ymin)/dely
      ycl = yco - (xmax-xmin)*min(1._R8,1._R8-frad*yloc)
      xc = xcl*cosxa + zcl*cosza + ycl*cosya
      yc = xcl*sinxa + zcl*sinza + ycl*sinya
  545 ihd = 1
      jx = min(xndivp,xndivp*(xc-xmapl)/(xmapr-xmapl))
      if(yc.lt.ymaxsf(jx) .and. yc.gt.yminsf(jx) ) go to 560
      if(yc.ge.ymaxsf(jx)) go to 555
      ihd = -1
      yminsf(jx) = yc
      if(yc.ge.ymaxsf(jx)) ymaxsf(jx) = yc
      if(ihdb.gt.0) go to 560
      call points(xcb,ycb,2,0,0,xc-xcb,yc-ycb)
      go to 560
  555 ihd = 0
      ymaxsf(jx) = yc
      if(yc.le.yminsf(jx)) yminsf(jx) = yc
      if(ihdb.gt.0) go to 560
      call setcrt(xcb,ycb)
      call vector(xc ,yc )
  560 ihdb = ihd
      xcb = xc
      ycb = yc
  530 continue
  520 continue
  470 continue
      if(ii.eq.6 .or. ii.eq.7) go to 500
!
! * * draw axes
!
      zcl = zco
      ycl = yco-(xmax-xmin)
      xcl = xco
      if(ii.eq.11) xcl = xcl + (xmax-xmin)
      xs = xcl*cosxa + zcl*cosza + ycl*cosya
      ys = xcl*sinxa + zcl*sinza + ycl*sinya
      xe = xs
      ycl = yco - (1._R8-frad)*(xmax-xmin)
      ye = xcl*sinxa + zcl*sinza + ycl*sinya
      div(1) = ymin
      div(2) = ymax
      call gaxisf(xs, ys, xe, ye, 1, 1, 1,"e10.2", 2, div)
      call setcrt(xs-0.05_R8*(xmax-xmin),ys)
      call vector(xs,ys)
      call vector(xe,ye)
      call vector(xe-0.05_R8*(xmax-xmin),ye)
      yf = ye
      xcl = xco + (xmax-xmin)
      if(ii.eq.11) xcl = xco
      ycl = yco - (xmax-xmin)
      xe = xcl*cosxa + zcl*cosza + ycl*cosya
      ye = ys
      div(1) = xmin
      div(2) = xmax
      if(ii.ne.11 ) then
      call gaxisf(xs, ys, xe, ye, 0, 1, 0,"f6.3", 2, div)
      cdiv(1) = 'major ra'
      cdiv(2) = 'dius (m)'
                   else
      call gaxisf(xe, ye, xs, ys, 0, 1, 1,"f6.3", 2, div)
      cdiv(1) = 'poloidal'
      cdiv(2) =   'flux '
                   endif
      if(ii.gt.16) then
      cdiv(1) = 'plate di'
      cdiv(2) = 'stance-m'
                   endif
      xc = xs - 0.20_R8*(xe-xs)
      yc = ys + 0.20_R8*(yf-ys)
      call setold(xc,yc,1,0,1,1)
      write(s100,2000) lab(ii)
      call gtext(s100,80,0)
 2000 format(2a8)
      xc = xs + 0.2_R8*(xe-xs)
      yc = ys - 0.35_R8*(yf-ys)
      call setold(xc,yc,1,0,1,0)
      write(s100,2000) cdiv(1),cdiv(2)
      call gtext(s100,80,0)
      xs = xe
      ys = ye
      zcl = zco + depth
      xe = xcl*cosxa + zcl*cosza + ycl*cosya
      ye = xcl*sinxa + zcl*sinza + ycl*sinya
      div(1) = timesv(1)
      div(2) = timesv(nrecor2)
      call gaxisf(xs, ys, xe, ye, 0, 1, 0,"f8.4", 0, div)
      xc = 0.5_R8*(xs+xe) + 0.20_R8*(xe-xs)
      yc = 0.5_R8*(ys+ye) - 0.15_R8*(ye-ys)
      call setold(xc,yc,1,0,1,0)
      cdiv(1) = 'time-sec'
      write(s100,2000) cdiv(1)
      call gtext(s100,80,0)
!
!.....label graphs
!
  500 if(ii.eq.7) go to 505
      if(acoef(901).eq.0._R8) then
         call map (xmin,xmax,ymin,ymax,.150_R8,.550_R8,.685_R8,1._R8)
         call setold(xmin,ymin-0.2_R8*(ymax-ymin),1,0,1,0)
         if(ii.ne.11 ) then
            if(ii.le.16) then
               write(s100,1301)
               call gtext(s100,80,0)
            endif
            if(ii.gt.16) then
               write(s100,1303)
               call gtext(s100,80,0)
            endif
         else
            write(s100,1302)
            call gtext(s100,80,0)
         endif
      endif
 
         if(acoef(901).gt.0._R8) then
            if(ii.eq.4) then
               call map(xmin,xmax,ymin,ymax,.1_R8,.45_R8,.55_R8,.9_R8)
               ylabset=ymax+(ymax-ymin)*.05_R8
               xlabset=xmin+(xmax-xmin)*.1_R8
               call setold(xlabset,ylabset,1,0,1,0)
               write(s100,13001)
13001          format("    p[kN/m^2] vs R[m]")
               call gtext(s100,80,0)
            endif
            if(ii.eq.8) then
               call map(xmin,xmax,ymin,ymax,.6_R8,.95_R8,.55_R8,.9_R8)
               ylabset=ymax+(ymax-ymin)*.05_R8
               xlabset=xmin+(xmax-xmin)*.1_R8
               call setold(xlabset,ylabset,1,0,1,0)
               write(s100,13002)
13002          format("   Jphi[kA/m^2] vs R[m]")
               call gtext(s100,80,0)
            endif
      endif
!
 1301 format(" major radius ")
 1302 format(" poloidal flux ")
 1303 format(" plate distance ")
!
      if(acoef(901).eq.0._R8) then
         call setold(xmax,ymax,1,0,1,0)
         write(s100,1484)
         call gtextm(s100,80,0,1,2)
 1484 format(/,6x,"time")
!
      lstop = 12
      if(nrecor2 .lt. 12) lstop=nrecor2
      lcolmn2 = nrecor2 - 11
      lcolmn3 = nrecor2 - 23
      do 485 l=1,lstop
      lcolmn2 = lcolmn2 - 1
      lcolmn3 = lcolmn3 - 1
      if(nrecor2.gt.12 .and. lcolmn2.gt.0 .and. lcolmn3.gt.0) then
      write(s100,1485) num(l),timesv(l),num(l+12),timesv(l+12),          &  
     &                 num(l+24),timesv(l+24)
      call gtext(s100,80,0)
                                                              endif
      if(nrecor2.gt.12 .and. lcolmn2.gt.0 .and. lcolmn3.le.0) then
      write(s100,1485) num(l),timesv(l),num(l+12),timesv(l+12)
      call gtext(s100,80,0)
                                                              endif
      if(nrecor2.le.12 .or. lcolmn2.le.0) then
      write(s100,1485) num(l),timesv(l)
      call gtext(s100,80,0)
                                          endif
  485 continue
 1485 format(1x,a2,1pe10.4,2x,a1,1pe10.4,2x,a1,1pe10.4)
      endif
!
!
      if(acoef(901).eq.0._R8) then
      if(ii.lt.10 .and. ii.ne.6) write(nsc1,1088) lab(ii)
      if(ii.ge.12 .and. ii.le.16) write(nsc1,1088) lab(ii)
      if(ii.eq.6) write(nsc1,1089) lab(ii+1),lab(ii)
      if(ii.eq.11 ) write(nsc1,1188) lab(ii)
      if(ii.gt.pnplt) write(nsc1,1288) lab(ii)
 1088 format(4x,a8,"  vs major radius")
 1089 format(1x,a8,"&",a8,"  vs major radius")
 1188 format(4x,a8,"  vs poloidal flux")
 1288 format(4x,a8,"  vs plate dis")
      call setld(3._R8,50._R8,0,0,2,1)
      call crtbcd(lab(ii))
      if(ii.eq.6) go to 400
      go to 510
      endif
  505 call map (xmin,xmax,ymin,ymax,.150_R8,.550_R8,.290_R8,.6_R8)
      call setld(2._R8,22._R8,0,0,2,1)
      call crtbcd(lab(ii))
!
      call setold(xmin+0.2_R8*(xmax-xmin),ymin-.25_R8*(ymax-ymin),1,0,1,  &  
     & 0)
      write(s100,1301)
      call gtext(s100,80,0)
      if(acoef(901).gt.0._R8) then
      if(ii.eq.8) write(nsc1,10088)
10088 format(" p and Jphi vs R")
      endif
!
  510 continue
!
      if(acoef(901).gt.0._R8) then
      call setld(1._R8,11._R8,1,0,1,0)
      plcprnt = tcurdtp*tpi*udsi
      write(s100,2005) gzero/xplas,plcprnt, fluxlink
 2005 format(" Bt[T]=",1pe10.3,"  Ip[A]  =",1pe10.3,                     &  
     &   " Flux Linkage[W]=",1pe10.3)
      call gtext(s100,80,0)
!
      alamda = ali2+betapol
      write(s100,2167) 2._R8*ali2,betapol,beta*100._R8
 2167 format(" Li(3)=",1pe10.3,"  Betapol=",1pe10.3,                     &  
     &   " Beta[%]        =",1pe10.3)
      call gtext(s100,80,0)
!
      write(s100,2168) rmajor,rminor
 2168 format(" R0[m]=   ",f7.3,"  a[m]   =   ",f7.3)
      call gtext(s100,80,0)
!
      write(s100,3168) shape5,shape3
 3168 format(" k    =   ",f7.3,"  d      =   ",f7.3)
      call gtext(s100,80,0)
!
      write(s100,21680) el95,del95
21680 format(" k95  =   ",f7.3,"  d95    =   ",f7.3)
      call gtext(s100,80,0)
!
      write(s100,2169) shape6,shape7
 2169 format(" Rx[m]=   ",f7.3,"  Zx[m]  =   ",f7.3)
      call gtext(s100,80,0)
!
      write(s100,2170) qprof2(1),global(15)
 2170 format(" q0   =   ",f7.3,"  q100   =   ",f7.3)
      call gtext(s100,80,0)
      write(s100,3170) global(61),global(62)
 3170 format(" q95  =   ",f7.3,"  qstar  =   ",f7.3)
      call gtext(s100,80,0)
!
      scar5 = -(psix(imag+1)-psix(imag-1))/(2._R8*deex*xary(imag))
      scar6 = -((psix(imag+1)-2._R8*psix(imag)+psix(imag-1))/deex**2     &  
     &           + scar5)/xary(imag)
      findex = -xary(imag)*scar6/(scar5+1.0E-33_R8)
      write(s100,3171) findex
      write(nterm,3171) findex
 3171 format(" nindx=   ",f7.3)
      call gtext(s100,80,0)
!      write(nsc1,1019)
! 1019 format(" poloidal flux")
!
      endif
 
      if(acoef(901).eq.0._R8) call frscj(6)
      if(acoef(901).gt.0._R8) then
      if(ii.eq.8) call frscj(6)
      endif
  400 continue
!...........................................................
!
!.....close disk files
!
!...........................................................
!     write(nterm,9331)
!9331 format("surfplot called")
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
