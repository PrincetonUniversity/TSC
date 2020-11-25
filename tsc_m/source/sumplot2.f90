!#include "f77_dcomplx.h"
      subroutine sumplot2
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}   ROS   }}}}}}}}}}}}}}}}}}}}}}
!
!.....produce a summary plot
!
      USE CLINAM
      USE SCRATCH
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ios,nxbig,nzbig,ii,nplrecm,irec,npl,icount,minc,mloop
      INTEGER imin2,imax2,jmin2,jmax2,i,j,imagt,jmagt,jj,isw,kk,k1
      INTEGER k2,n,l
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cval,div
      REAL*8 pds,pfx,pfy,pfxy,psik
      REAL*8 xmin,xmax,zmin,zmax,pmin,xmg,zmg,ptemp
      REAL*8 pinterp,ptempp,xzmg,xx,dum1,dum2,dum3,dum4,psval,dum5
      REAL*8 dum6,dum7,sym
      REAL*8 AREAL
!============
      dimension cval(10),div(10)
!     dimension tary(pnsave),xmary(pnsave),zmary(pnsave),
!    1     asfary(pnsave),abary(pnsave),aipary(pnsave),aefary(pnsave)
!     dimension xbig(4*pnx),zbig(4*pnz),pbig(4*pnx,4*pnz)
      dimension pds(6),pfx(6),pfy(6),pfxy(6,6),psik(10)
!     dimension xpsiko(10,4*pnz),xpsiki(10,4*pnz),zpsiko(10,4*pnz),
!    +zpsiki(10,4*pnz)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xmary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zmary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: asfary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: abary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aipary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aefary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xbig
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zbig
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pbig
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xpsiko
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xpsiki
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: zpsiko
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: zpsiki
!============      
      IF(.not.ALLOCATED(tary)) ALLOCATE( tary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(xmary)) ALLOCATE( xmary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(zmary)) ALLOCATE( zmary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(asfary)) ALLOCATE( asfary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(abary)) ALLOCATE( abary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(aipary)) ALLOCATE( aipary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(aefary)) ALLOCATE( aefary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(xbig)) ALLOCATE( xbig(4*pnx), STAT=istat)
      IF(.not.ALLOCATED(zbig)) ALLOCATE( zbig(4*pnz), STAT=istat)
      IF(.not.ALLOCATED(pbig)) ALLOCATE( pbig(4*pnx,4*pnz), STAT=istat)
      IF(.not.ALLOCATED(xpsiko)) ALLOCATE( xpsiko(10,4*pnz), STAT=istat)
      IF(.not.ALLOCATED(xpsiki)) ALLOCATE( xpsiki(10,4*pnz), STAT=istat)
      IF(.not.ALLOCATED(zpsiko)) ALLOCATE( zpsiko(10,4*pnz), STAT=istat)
      IF(.not.ALLOCATED(zpsiki)) ALLOCATE( zpsiki(10,4*pnz), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : sumplot2  ' 
!============      
!
      if( numargs .lt. 1 ) then
         filename = 'solout'
      else
         filename = 'solout' // '.' // trim(suffix)
      end if
      open(90,file=trim(filename),status='unknown',iostat=ios)
!
      if(ifrst(4) .eq. 1) return
!
      xmin = acoef(751)
      xmax = acoef(752)
      zmin = acoef(753)
      zmax = acoef(754)
      if(xmin.lt.ccon .or. xmax .gt. alx .or.                            &  
     &   zmin.lt.-alz .or. zmax .gt. alz .or.                            &  
     &   xmin.ge.xmax .or. zmin .ge. zmax) return
      call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.285_R8,1._R8)
!
      if(acoef(755) .gt. 4._R8) return
!
      if(acoef(755) .lt. 1.0_R8) go to 997
      nxbig = int(acoef(755))*nx
      nzbig = int(acoef(755))*nz
      do 999 ii=1,nxbig
 999  xbig(ii) = xlim + (ii-1)*(xlim2-xlim)/(nxbig-1)
      do 998 ii=1,nzbig
 998  zbig(ii) = zlim*(isym-1) + (ii-1)*((2-isym)*zlim)/(nzbig-1)
 997  continue
!
      nplrecm = nplrec
      rewind nsc2
      irec = pnx*pnz
      do 100 npl=1,nplrecm
      call bufin(nsc2,vecx(1,1),vecx(pnx,pnz))
      icount = icount+irec
!
!.....start coding to find mag axis in vecx array
      pmin = 1.E20_R8
      minc = 2
      mloop = 0
!
      imin2 = nx/2 - minc
      imax2 = nx/2 + minc
      jmin2 = nz/2 - minc
      if(isym.eq.1) jmin2 = 2
      jmax2 = nz/2 + minc
  121 continue
      if(imin2 .lt. iminn+2) imin2 = iminn+2
      if(imax2 .gt. imaxx-2) imax2 = imaxx-2
      if(jmin2 .lt. jminn  ) jmin2 = jminn
      if(jmax2 .gt. jmaxx  ) jmax2 = jmaxx
!
!
      do 120 i=imin2,imax2
      do 130 j=jmin2,jmax2
      if(iexv(i,j).eq.1) go to 130
      if(vecx(i,j).gt.pmin) go to 130
      imagt = i
      jmagt = j
      pmin = vecx(i,j)
  130 continue
  120 continue
      if(imagt.ne.imin2) go to 131
      imin2 = imin2 -minc
      mloop = mloop + 1
      if(mloop.lt.nx) go to 121
      stop
  131 continue
      if(imagt.ne.imax2) go to 132
      imax2 = imax2 +minc
      mloop = mloop + 1
      if(mloop.lt.nx) go to 121
      stop
  132 continue
      if(isym.eq.1) go to 133
      if(jmagt.ne.jmin2) go to 133
      jmin2 = jmin2 -minc
      mloop = mloop + 1
      if(mloop.lt.nx) go to 121
      stop
  133 continue
      if(jmagt.ne.jmax2) go to 134
      jmax2 = jmax2 +minc
      mloop = mloop + 1
      if(mloop.lt.nx) go to 121
      stop
  134 continue
      call axm2d2(vecx,penx,imagt,jmagt,xary(imagt),zary(jmagt),         &  
     &xmg,zmg,xmg,zmg,pds,pfx,pfy,pfxy,6,deex,deez,0)
!
!.....end of coding to find magnetic axis in vecx array, xmg,zmg
!
      do 993 ii=nx/2,nx
      ptemp = pinterp(xary(ii),zmg,0,0)
      ptempp = pinterp(xary(ii+1),zmg,0,0)
      if(ptemp .le. vecx(1,1) .and. ptempp .gt. vecx(1,1)) then
      xzmg = xary(ii) + (xary(ii+1)-xary(ii))/(ptempp-ptemp)             &  
     &*(vecx(1,1)-ptemp)
      endif
 993  continue
      psik(1)=vecx(1,1)
      do 992 ii=2,7
      xx = xzmg + AREAL(ii-1)*0.01_R8
      psik(ii) = pinterp(xx,zmg,0,0)
 992  continue
!
      if(acoef(755) .lt. 1.0_R8) go to 994
      do 996 ii=1,nxbig
      do 995 jj=1,nzbig
      isw = 0
      if(ii.eq.1 .and. jj.eq.1) isw = 1
      call grap2(1,zbig(jj),xbig(ii),dum1,dum2,dum3,dum4,psval,          &  
     &dum5,dum6,dum7,isw)
      pbig(ii,jj)=psval
 995  continue
 996  continue
 994  continue
!
      if(acoef(755) .lt. 1.0_R8) then
      nxbig=nx
      nzbig=nz
      do 988 ii=2,nx-1
 988  xbig(ii)=xary(ii)
      do 987 ii=2,nz-1
 987  zbig(ii)=zary(ii)
      do 986 ii=2,nx-1
      do 985 jj=2,nz-1
      pbig(ii,jj)=vecx(ii,jj)
 985  continue
 986  continue
      endif
      do 991 ii=2,nzbig-1
      do 990 jj=2,nxbig-1
      do 989 kk=1,7
      if(pbig(jj,ii) .le. psik(kk) .and. pbig(jj+1,ii) .gt.              &  
     &psik(kk)) then
      xpsiko(kk,ii)=xbig(jj)+(xbig(jj+1)-xbig(ii))/(pbig(jj+1,ii)-       &  
     &pbig(jj,ii))*(psik(kk)-pbig(jj,ii))
      zpsiko(kk,ii)=zbig(ii)
      endif
      if(pbig(jj,ii) .gt. psik(kk) .and. pbig(jj+1,ii) .le.              &  
     &psik(kk)) then
      xpsiki(kk,ii)=xbig(jj+1)+(xbig(jj)-xbig(jj+1))/(pbig(jj,ii)        &  
     &-pbig(jj+1,ii))*(psik(kk)-pbig(jj+1,ii))
      zpsiki(kk,ii)=zbig(ii)
      endif
 989  continue
 990  continue
 991  continue
!
      do 983 kk=1,7
      write(90,984) (xpsiki(kk,ii),zpsiki(kk,ii),xpsiki(kk,ii+1),        &  
     &zpsiki(kk,ii+1),kk,ii=2,nzbig-1)
      write(90,984) (xpsiko(kk,ii),zpsiko(kk,ii),xpsiko(kk,ii+1),        &  
     &zpsiko(kk,ii+1),kk,ii=2,nzbig-1)
 984  format(4e15.5,i5)
 983  continue
!
      sym = 1._R8
      do 50 j=2,nzp
   50 vzy(j) = sym*zary(j)
!
!.....define k1,k2,cval to plot psilim contour only
      k1 = 1
      k2 = 1
      cval(1) = vecx(1,1)
      if(acoef(755) .lt. 1.0_R8) then
      call rcontr(k1,cval,k2,vecx,penx,xary,iminy(npl),imaxy(npl),1,     &  
     &                         vzy,jminy(npl),jmaxy(npl),1)
      else
      call rcontr(k1,cval,k2,pbig,4*pnx,xbig,1,nxbig,1,                  &  
     &                        zbig,1,nzbig,1)
      endif
  100 continue
!     do 81 ii=1,nwire
!  81 call coildr(iwire(ii),jwire(ii),ilwire(ii))
      if(iplate.eq.0 .or. nplate.eq.0) go to 8100
      do 8110 n=1,nplate
      call setcrt(xsega(n,1),zsega(n,1))
      do 8111 l=2,nseg(n)+1
 8111 call vector(xsega(n,l),zsega(n,l))
 8110 continue
 8100 continue
!
      call frscj(9)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
