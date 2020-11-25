!#include "f77_dcomplx.h"
          subroutine diagbol(radpow,nd)
!***************************************************************************
!
!...subroutine to calculate the radiation collected at bolometer
!...detectors by summing along pre-defined chords through the plasma
!...ckessel  1997
!
!         rw          -          radial position of window
!         zw          -          vertical position of window
!         dpolw       -          window dimension in poloidal direction
!         dtorw       -          window dimension in toroidal direction
!         ppolw       -          projected window dimension in poloidal
!                                direction
!         ptorw       -          projected window dimension in toroidal
!                                direction
!         thetw       -          theta of window normal
!         trans       -          transmission factor
!         nd          -          number of detectors
!         radd        -          radial distance from window center to
!                                detector
!         aread       -          detector area for intercepting radiation
!         thetd       -          poloidal theta of detector (measured from
!                                outboard midplane)
!         thetd2      -          toroidal theta of detector from center
!                                of window
!         thets       -          spread angle of view from detector
!                                through window
!         ns          -          number of segments used along a detector
!                                chord
!         ds          -          segment length along detector chord
!         rs          -          radial point along detector chord
!         zs          -          vertical point along detector chord
!         psis        -          flux at point (rs,zs) along detector chord
!         radpow      -          total radiated power intercepted by detector
!
!***NOTE: this subroutine will not give correct answers for isym=1***
!
!
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nd,i,ns,j,l,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 radpow,thetd,thetd2,ppolw,tanth,tanph,coswd,ptorw,rs
      REAL*8 zs,psis,fv,rw,zw,aread,dpolw,dtorw,thetw,trans,radd,ds
      REAL*8 rso,zso,srad,thets,dsfac,sumfv,pinterp,vlfac,sourcp
      REAL*8 sourcm,sourc,areax,sa,tant,tanp,ddx,ddy,radfac
      REAL*8 AREAL
!============
      dimension thetd(50),thetd2(50)
      dimension ppolw(50),tanth(50),tanph(50),coswd(50),ptorw(50)
      dimension rs(10),zs(10),psis(10),fv(10)
      dimension radpow(nd)
!
!...window data
      data rw/3.104_R8/
      data zw/1.477_R8/
      data aread/0.0001_R8/
      data dpolw/0.01_R8/
      data dtorw/0.04_R8/
      data thetw/225.0_R8/
      data trans/0.3_R8/
!
!...bolometer detector data
      data radd/0.20_R8/
!     data thetd/-34.54,-38.63,-42.87,-46.96,-51.20,-55.30,-59.54,
!    .-63.63,-67.87,-71.96,-76.20,-80.30,-84.54,-88.63,-92.87,
!    .-96.96,-101.20,-105.29,-109.54,-90.75,-94.92/
      data (thetd(i),i=1,21)                                             &  
     &/214.54_R8,218.63_R8,222.87_R8,226.96_R8,231.20_R8,235.30_R8,      &  
     & 239.54_R8,                                                        &  
     & 243.63_R8,246.87_R8,251.96_R8,256.20_R8,260.30_R8,264.54_R8,      &  
     & 268.63_R8,272.87_R8,                                              &  
     & 276.96_R8,281.20_R8,285.29_R8,289.54_R8,270.75_R8,174.92_R8/
      data thetd2/0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,  &  
     & 0.000_R8,0.000_R8,0.000_R8,                                       &  
     & 0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,   &  
     & 0.000_R8,0.000_R8,0.000_R8,0.000_R8,                              &  
     & 0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,   &  
     & 0.000_R8,0.000_R8,0.000_R8,0.000_R8,                              &  
     & 0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,   &  
     & 0.000_R8,0.000_R8,0.000_R8,0.000_R8,                              &  
     & 0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,0.000_R8,   &  
     & 0.000_R8/
!
!...number of segments along detector chord used in summing up total
!...radiation to detector
      data ns/100/
!============      
!
!
!...calculate projected area of window, tangent of theta (poloidal),
!...and tangent of phi (toroidal)
      do 100 i=1,nd
      ppolw(i)=dpolw*cos((thetd(i)-thetw)/(2._R8*pi))
      ptorw(i)=dtorw*cos(thetd2(i)/(2._R8*pi))
      tanth(i)=ppolw(i)/radd/2._R8
      tanph(i)=ptorw(i)/radd/2._R8
 100  continue
!
!...calculate length of segment used in summing up total radiation
!...to detector
      ds=sqrt((xlim2-xlim)**2+(zlim+zlim)**2)/AREAL(ns)
!
!...loop over detectors
      do 300 i=1,nd
      rso=rw
      zso=zw
      srad=0.0_R8
      thets=(360._R8/(2._R8*pi))*atan(tanth(i))
      dsfac=sqrt(1._R8+tanth(i)**2)
!...loop over segments along a detector chord
      do 301 j=1,ns
!...(rs(1),zs(1)) is the central point for the segment volume
      rs(1)=rso+ds*cos((thetd(i))/(2._R8*pi))
      zs(1)=zso+ds*sin((thetd(i))/(2._R8*pi))
      rs(2)=rso+(ds/2._R8)*cos((thetd(i))/(2._R8*pi))
      zs(2)=zso+(ds/2._R8)*sin((thetd(i))/(2._R8*pi))
      rs(3)=rso+(3._R8*ds/2._R8)*cos((thetd(i))/(2._R8*pi))
      zs(3)=zso+(3._R8*ds/2._R8)*sin((thetd(i))/(2._R8*pi))
      rs(4)=rso+ds*dsfac*cos((thetd(i)-thets)/(2._R8*pi))
      zs(4)=zso+ds*dsfac*sin((thetd(i)-thets)/(2._R8*pi))
      rs(5)=rso+(ds/2._R8)*dsfac*cos((thetd(i)-thets)/(2._R8*pi))
      zs(5)=zso+(ds/2._R8)*dsfac*sin((thetd(i)-thets)/(2._R8*pi))
      rs(6)=rso+(3._R8*ds/2._R8)*dsfac*cos((thetd(i)-thets)/(2._R8*pi))
      zs(6)=zso+(3._R8*ds/2._R8)*dsfac*sin((thetd(i)-thets)/(2._R8*pi))
      rs(7)=rso+ds*dsfac*cos((thetd(i)+thets)/(2._R8*pi))
      zs(7)=zso+ds*dsfac*sin((thetd(i)+thets)/(2._R8*pi))
      rs(8)=rso+(ds/2._R8)*dsfac*cos((thetd(i)+thets)/(2._R8*pi))
      zs(8)=zso+(ds/2._R8)*dsfac*sin((thetd(i)+thets)/(2._R8*pi))
      rs(9)=rso+(3._R8*ds/2._R8)*dsfac*cos((thetd(i)+thets)/(2._R8*pi))
      zs(9)=zso+(3._R8*ds/2._R8)*dsfac*sin((thetd(i)+thets)/(2._R8*pi))
      rso=rs(1)
      zso=zs(1)
!...determine number of points out of 9 that exist inside the plasma
!...and weight volume accordingly
      sumfv=0._R8
      do 303 l=1,9
      fv(l)=1._R8
      psis(l)=pinterp(rs(l),zs(l),0,0)
      if(psis(l) .gt. psilim) fv(l)=0._R8
      sumfv=sumfv+fv(l)
 303  continue
      if(sumfv .eq. 0._R8) go to 301
      vlfac=sumfv/9._R8
!...sum over all radiation sources in a given segment volume
      do 302 k=1,npsit-1
      if(psis(1) .ge. xsv2(k) .and. psis(1) .lt. xsv2(k+1)) then
      sourcp=sradion(k+1)+savebre(k+1)+savecyc(k+1)
      sourcm=sradion(k)+savebre(k)+savecyc(k)
      sourc=sourcm+((sourcp-sourcm)/(xsv(k+1)-xsv(k))*                   &  
     &(psis(1)-xsv(k)))
      endif
 302  continue
      areax=aread
      sa=sqrt((rs(1)-rw)**2+(zs(1)-zw)**2)
      tant=ppolw(i)/(2._R8*sa)
      tanp=ptorw(i)/(2._R8*sa)
      ddx=2._R8*(sa+radd)*tant
      ddy=2._R8*(sa+radd)*tanp
      if(ddx .lt. sqrt(aread)) areax=sqrt(aread)*ddx
      if(ddy .lt. sqrt(aread)) areax=sqrt(aread)*ddy
      if(ddx .lt. sqrt(aread) .and. ddy .lt. sqrt(aread)) areax=ddx*ddy
      radfac=4._R8*tanth(i)*tanph(i)*areax*ds/(4._R8*pi)
      srad=srad+radfac*sourc*vlfac
 301  continue
      radpow(i)=trans*srad
 300  continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
