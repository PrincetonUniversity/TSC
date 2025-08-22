      subroutine plotobs(tmin,tmaxx,big1,big2,big3,big4,big5)
!     ------------------
!.....plot observation point time historys
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,iflopx,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tmaxx,big1,big2,big3,big4,big5,tmin,xmin,xmax,ymax
      REAL*8 ymin,x1,x2,y1,y2,amval,abval,amval10,abval10,dely,ybar
      REAL*8 gam1,gam100,gam2,gam200,tau1,tau100,tau2,tau200,szrms
      REAL*8 zzrms,tauz100,ymx,ymn
!============
      dimension big1(1),big2(1),big3(1),big4(1),big5(1)
!============      
!
      xmin = tmin
      xmax = tmaxx
!.....begin big loop:
      do 8051 n=1,nobs-1,2
      if(npltobs(n).eq.0) go to 8051
!       do plots twice for flux and voltage if(iwayne.gt.0)
       iflopx = 1
 9999  continue
!
!.....plot unprocessed fluxes:
      ymax = -1.E40_R8
      ymin =  1.E40_R8
      do 6001 i=1,nrecord
!      write unprocessed fluxes in array big1, big3:
!      time is stored in big2:
       big1(i) = fluxu(n,i)
       if (iflopx .eq. 2 .and. i .lt. 15) big1(i) = 0._R8
       if (iflopx .eq. 2 .and. i .ge. 15)                                &  
     &big1(i)=6.2832_R8*(big1(i)-fluxu(n,i-5))/(big2(i)-big2(i-5))
       big3(i) = fluxl(n,i)
       if (iflopx .eq. 2 .and. i .lt. 15) big3(i) = 0._R8
       if (iflopx .eq. 2 .and. i .ge. 15)                                &  
     &big3(i)=6.2832_R8*(big3(i)-fluxl(n,i-5))/(big2(i)-big2(i-5))
       ymax = max(ymax,big1(i))
       ymin = min(ymin,big1(i))
       ymax = max(ymax,big3(i))
       ymin = min(ymin,big3(i))
 6001 continue
!
      x1 = .1_R8
      x2 = .3_R8
      y1 = .7_R8
      y2 = 1._R8
      if(iflopx.eq.2) then
      x1 = .2_R8
      x2 = .8_R8
      y1 = .4_R8
      y2 = 1._R8
                      endif
      if( abs(ymax-ymin) .le. 1.E-8_R8*max(abs(ymax),abs(ymin)))         &  
     &go to 8000
      call mapg(xmin,xmax,ymin,ymax,x1,x2,y1,y2)
      call setold(xmin-(xmax-xmin)*.40_R8,ymin,1,0,1,1)
      if (iflopx .eq. 1) then
      write(s100,9001)
      call gtext(s100,80,0)
                         endif
 9001 format("unprocessed fluxes")
       if (iflopx .eq. 2) then
      write(s100,9099)
      call gtext(s100,80,0)
                          endif
 9099  format("unprocessed volts ")
      call trace(big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
      call tracep(big2(1),big3(1),nrecord,-1,-1,-1)
 8000 continue
       if (iflopx .eq. 2) go to 9998
!
!.....compute least squares fit to growth rate
!
!     write flux difference in big1:
      do 9006 i=1,nrecord
       big1(i) = log(abs(fluxu(n,i)-fluxl(n,i)) + 1.E-12_R8)
 9006 continue
!     determine growth rates from flux differences:
!     big4: fit over all time points, big5: last 10
      call fitls(big2,big1,nrecord,big4,big5,                            &  
     &           amval,abval,amval10,abval10,10)
!
!.....plot and write growth rates from flux differences:
!
       ymax = max(big4(nrecord),big5(nrecord))
       ymin = min(big4(nrecord),big5(nrecord))
       dely = ymax - ymin
       ybar = (ymax + ymin)*0.5_R8
       ymax = ybar + 5._R8*dely
       ymin = ybar - 5._R8*dely
       do 9008 i=1,nrecord
        big4(i) = min(ymax,big4(i))
        big4(i) = max(ymin,big4(i))
        big5(i) = min(ymax,big5(i))
        big5(i) = max(ymin,big5(i))
 9008  continue
!
      gam1 = amval
      gam100 = amval10
!
      if( abs(ymax-ymin) .le. 1.E-8_R8*max(abs(ymax),abs(ymin)))         &  
     &go to 8001
      call mapg(xmin,xmax,ymin,ymax,0.5_R8,0.7_R8,.7_R8,1._R8)
      call setold(xmin-(xmax-xmin)*.70_R8,ymin,1,0,1,1)
      write(s100,7777)
      call gtext(s100,80,0)
 7777 format('flux differences:')
      call setold(xmin-(xmax-xmin)*.50_R8,ymin,1,0,1,1)
      write(s100,8009) gam1
      call gtext(s100,80,0)
 8009 format("gam1   =",1pe12.4," sec-1")
      call setold(xmin-(xmax-xmin)*.40_R8,ymin,1,0,1,1)
      write(s100,8010) gam100
      call gtext(s100,80,0)
 8010 format("gam100 =",1pe12.4," sec-1")
!
      call trace(big2(1),big4(1),nrecord,-1,-1,0._R8,0._R8)
      call tracep(big2(1),big5(1),nrecord,-1,-1,-1)
 8001 continue
!
      write(nout,8011) n,gam1,gam100
 8011 format(' n =',i3,3x,'gam1 =',1pe12.4,3x,'gam100 =',1pe12.4)
!
!.....fits to flux differences:
!     write fit to flux differences in big4, big5:
!     flux difference is in big1:
      ymax = -1.E40_R8
      ymin =  1.E40_R8
      do 8007 i=1,nrecord
       big1(i) = abs(fluxu(n,i)-fluxl(n,i))
       big4(i) = exp(amval  *big2(i) + abval)
       big5(i) = exp(amval10*big2(i) + abval10)
       ymax = max(ymax,big1(i))
       ymin = min(ymin,big1(i))
 8007 continue
       if(ymin.lt.1.E-3_R8*ymax) ymin = 1.E-3_R8*ymax
      do 8015 i=1,nrecord
       big1(i) = max(ymin,big1(i))
 8015 continue
!
!.....plot fits to flux differences:
      if( abs(ymax-ymin) .le. 1.E-8_R8*max(abs(ymax),abs(ymin)))         &  
     &go to 8002
      if(ymax.le.0 .or. ymin.le.0) go to 8002
      call mapgsl(xmin,xmax,ymin,ymax,0.8_R8,1.0_R8,.7_R8,1._R8)
      call trace(big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
!     call trace(big2(1),big4(1),nrecord,-1,-1,0.,0.)
      call tracep(big2(1),big5(1),nrecord,-1,-1,-1)
 8002 continue
!
!.....compute growth rates from fits to the time derivatives
!     of flux differences:
!              big2: time
!              big4: contains fitted curve (all points)
!              big5: contains fitted curve (last 10 points)
!
      do 300 i=1,nrecord
       big1(i)= fluxu(n,i)-fluxl(n,i)
 300  continue
      call dif1(big2,big1,big3,nrecord)
!     store derivative of flux difference in big1:
      do 310 i=1,nrecord
       big1(i) = big3(i)/udst
       big3(i) = abs(big3(i))/udst
       big3(i) = log(big3(i) + 1.E-12_R8)
 310  continue
!     determine growth rates from flux time derivative:
!     (fit over all and last 10 points)
      call fitls(big2,big3,nrecord,big4,big5,                            &  
     &           amval,abval,amval10,abval10,10)
!
!.....plot and write growth rates determined from flux time derivatives:
      ymax = -1.E40_R8
      ymin =  1.E40_R8
!
       ymax = max(big4(nrecord),big5(nrecord))
       ymin = min(big4(nrecord),big5(nrecord))
       dely = ymax - ymin
       ybar = (ymax + ymin)*0.5_R8
       ymax = ybar + 10._R8*dely
       ymin = ybar - 10._R8*dely
!
       do 9009 i=1,nrecord
        big4(i) = min(ymax,big4(i))
        big4(i) = max(ymin,big4(i))
        big5(i) = min(ymax,big5(i))
        big5(i) = max(ymin,big5(i))
 9009  continue
!
      gam2   = amval
      gam200 = amval10
!
      if( abs(ymax-ymin) .le. 1.E-8_R8*max(abs(ymax),abs(ymin)))         &  
     &go to 8003
!     if(ymin.le.0 .or. ymax.le.0) go to 8003
      call mapg(xmin,xmax,ymin,ymax,0.5_R8,0.7_R8,.30_R8,.60_R8)
      call setold(xmin-(xmax-xmin)*.70_R8,ymin,1,0,1,1)
      write(s100,7778)
      call gtext(s100,80,0)
 7778 format('flux derivative:')
      call setold(xmin-(xmax-xmin)*.50_R8,ymin,1,0,1,1)
      write(s100,8039) gam2
      call gtext(s100,80,0)
 8039 format("gam2   =",1pe12.4," sec-1")
      call setold(xmin-(xmax-xmin)*.40_R8,ymin,1,0,1,1)
      write(s100,8040) gam200
      call gtext(s100,80,0)
 8040 format("gam200 =",1pe12.4," sec-1")
!
      call trace(big2(1),big4(1),nrecord,-1,-1,0._R8,0._R8)
      call tracep(big2(1),big5(1),nrecord,-1,-1,-1)
 8003 continue
!
      write(nout,8041) n,gam2  ,gam200
 8041 format(' n =',i3,3x,'gam2   =',1pe12.4,3x,'gam200 =',1pe12.4)
!
!.....fits to the time derivatives of flux differences:
!     write fit to flux differences in big4, big5,
!                            d(psi)/dt is in big1:
      ymax = -1.E40_R8
      ymin =  1.E40_R8
      do 8042 i=1,nrecord
       big1(i) = abs(big1(i))
       big4(i) = exp(amval  *big2(i) + abval)
       big5(i) = exp(amval10*big2(i) + abval10)
       ymax = max(ymax,big1(i))
       ymin = min(ymin,big1(i))
 8042 continue
      if(ymin.lt.1.E-3_R8*ymax) ymin = 1.E-3_R8*ymax
      do 8043 i=1,nrecord
       big1(i) = max(ymin,big1(i))
 8043 continue
!
!.....plot fit to the time derivative of flux differences:
      if( abs(ymax-ymin) .le. 1.E-8_R8*max(abs(ymax),abs(ymin)))         &  
     &go to 8004
      if(ymax.le.0 .or. ymin.le.0) go to 8004
      call mapgsl(xmin,xmax,ymin,ymax,0.8_R8,1.0_R8,.30_R8,.60_R8)
      call trace(big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
!     call trace(big2(1),big4(1),nrecord,-1,-1,0.,0.)
      call tracep(big2(1),big5(1),nrecord,-1,-1,-1)
 8004 continue
!
!.....write additional informations:
!
      tau1 = 0._R8
      tau100 = 0._R8
      tau2 = 0._R8
      tau200 = 0._R8
      if(gam1.ne.0._R8) tau1 = 1._R8/gam1
      if(gam100.ne.0._R8) tau100 = 1._R8/gam100
      if(gam2  .ne.0._R8) tau2 = 1._R8/gam2
      if(gam200.ne.0._R8) tau200 = 1._R8/gam200
!
      call setld(47._R8,27._R8,1,0,1,0)
      write(s100,9050)
      call gtext(s100,80,0)
 9050 format('growth rates')
      call setld(74._R8,27._R8,1,0,1,0)
      write(s100,9051)
      call gtext(s100,80,0)
 9051 format('fits')
      call setld(1._R8,8._R8,1,0,1,0)
      write(s100,9052)
      call gtext(s100,80,0)
 9052 format('observation pair:')
      write(s100,9004) xobs(n),zobs(n),xobs(n+1),zobs(n+1)
      call gtext(s100,80,0)
 9004 format('x,y=(',f5.2,',',f5.2,'),(',f5.2,',',f5.2,')' )
      write(s100,9053)
      call gtext(s100,80,0)
 9053 format('      ')
!     write(s100,9054) tau1
!     call gtext(s100,80,0)
 9054 format('tau1 = ',1pe12.4,' sec')
      write(s100,9055) tau100
      call gtext(s100,80,0)
 9055 format('tau100 = ',1pe12.4,' sec (flux fit over 10 points)')
!     write(s100,9056) tau2
!     call gtext(s100,80,0)
 9056 format('tau2 = ',1pe12.4,' sec')
      write(s100,9057) tau200
      call gtext(s100,80,0)
 9057 format('tau200 = ',1pe12.4,' sec (flux der over 10 points)')
      write(nterm,9055) tau100
      write(nterm,9057) tau200
!
!
!..calculate growth rate from zmag, and plot
      szrms=0.0_R8
      do 7000 i=2,nrecord
      big1(i)=log(abs(pltsav((i-1)*lenscr+7))+1.E-12_R8)
      szrms=szrms+(pltsav((i-1)*lenscr+7))**2                            &  
     &*(big2(i)-big2(i-1))
7000  continue
      big1(1)=big1(2)
      zzrms=sqrt(abs(szrms/(big2(nrecord)-big2(1))))
      call fitls(big2,big1,nrecord,big4,big5,                            &  
     &   amval,abval,amval10,abval10,10)
      if(amval10.ne.0._R8) tauz100=1._R8/amval10
      call setld(1._R8,10._R8,1,0,1,0)
      write(s100,9058) tauz100,zzrms
      call gtext(s100,80,0)
9058  format('tauz100 = ',1pe12.4,' sec',2x,'Z(RMS)= ',1pe12.4)
      ymx=-1.E40_R8
      ymn=1.E40_R8
      do 7001 i=1,nrecord
      big1(i)=abs(pltsav((i-1)*lenscr+7))
      big4(i)=exp(amval*big2(i)+abval)
      big5(i)=exp(amval10*big2(i)+abval10)
      ymx=max(ymx,big1(i))
      ymn=min(ymn,big1(i))
7001  continue
      if(ymn.lt.1.E-3_R8*ymx) ymn=1.E-3_R8*ymx
      do 7002 i=1,nrecord
      big1(i)=max(ymn,big1(i))
7002  continue
      if( abs(ymx-ymn) .le. 1.E-8_R8*max(abs(ymx),abs(ymn)))             &  
     &go to 8005
      if(ymn.le.0.0_R8.or. ymx.le.0.0_R8) go to 8005
      call mapgsl(xmin,xmax,ymn,ymx,0.1_R8,0.3_R8,0.3_R8,0.6_R8)
      call setold(xmin-(xmax-xmin)*.40_R8,ymn,1,0,1,1)
      write(s100,9037)
      call gtext(s100,80,0)
9037  format("zmag")
      call trace(big2(1),big1(1),nrecord,-1,-1,0._R8,0._R8)
!     call trace(big2(1),big4(1),nrecord,-1,-1,0.,0.)
      call tracep(big2(1),big5(1),nrecord,-1,-1,-1)
! kdm for moviedata.f90
      obsdiag(1,(n+1)/2)=big1(1)
      obsdiag(2,(n+1)/2)=big1(2)
 9998 if(iflopx.eq.2) then
      write(s100,9991) xobs(n),zobs(n),xobs(n+1),zobs(n+1)
      call gtextm(s100,80,0,1,3)
                      endif
 9991  format ("xy/",f5.2,",",f5.2,"/",f5.2,",",f5.2)
 8005 continue
      nn = (n+1)/2
      write(nsc1,1019) nn
 1019 format(" observation pair ",i3)
      call frscj(6)
      if(iwayne.eq.0) go to 9996
       if (iflopx .eq. 1) go to 9995
       if (iflopx .eq. 2) go to 9996
 9995  iflopx=2
       go to 9999
 9996  continue
!.....end of big loop
8051  continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
