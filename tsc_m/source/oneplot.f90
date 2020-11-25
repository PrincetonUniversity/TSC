      subroutine oneplot(dumx,plotmin,plotmax,iminp,imaxp,               &  
     &                                         jminp,jmaxp)
!
!.....produce a single plot
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iminp,imaxp,jminp,jmaxp,itmax,it,j,k1,k2
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 plotmin,plotmax,dumx,cval,xmin,zmin,xmax,zmax,xave
      REAL*8 xmaxpl,xminpl,zmaxpl,zminpl,sym
!============
      dimension dumx(penx,penz)
      dimension cval(10)
!============      
!
      xmin = ccon
      zmin = -alz
      xmax = alx
      zmax = alz
!
      xave = .5_R8*(xmax+xmin)
      if((zmax-zmin).gt.xmax-xmin) xmax = xave+.5_R8*(zmax-zmin)
      if((zmax-zmin).gt.xmax-xmin) xmin = xave-.5_R8*(zmax-zmin)
      if((xmax-xmin).gt.(zmax-zmin)) zmax = .5_R8*(xmax-xmin)
      if((xmax-xmin).gt.(zmax-zmin)) zmin = -.5_R8*(xmax-xmin)
      xmaxpl=xmax+.05_R8*(xmax-xmin)
      xminpl=xmin-.05_R8*(xmax-xmin)
      zmaxpl=zmax+.05_R8*(zmax-zmin)
      zminpl=zmin-.05_R8*(zmax-zmin)
      call maps(xminpl,xmaxpl,zminpl,zmaxpl,.200_R8,.800_R8,.400_R8,     &  
     & 1._R8)
!
      itmax = 1+isym
      do 100 it=1,itmax
      sym = 1._R8
      if(it.eq.2) sym = -1._R8
      do 50 j=2,nzp
   50 vzy(j) = sym*zary(j)
!
!.....define k1,k2,cval to plot psilim contour only
      k1 = -10
      k2 = 0
      cval(1) = plotmin
      cval(2) = plotmax
      call rcontr(k1,cval,k2,dumx,penx,xary,iminp,imaxp,1,               &  
     &                         vzy,jminp,jmaxp,1)
  100 continue
!
      call frscj(17)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
