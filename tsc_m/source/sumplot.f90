      subroutine sumplot
!
!.....produce a summary plot
!
 
      USE CLINAM
      USE SCRATCH
      USE EQRUNS
      USE ITER1

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,nplrecm,irec,npl,icount,itmax,it,j,k1,k2,ii,l
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cval,div
      REAL*8 xmin,zmin,xmax,zmax,xave,xmaxpl,xminpl
      REAL*8 zmaxpl,zminpl,sym,dxc,dzc,x1,z1,x2,z2,x3,z3,x4,z4
      REAL*8 z1neg,zlneg,xp1,zp1,xp2,zp2,xp3,zp3,xp4,zp4,xp5,zp5
      REAL*8 xn1,zn1,xn2,zn2,xn3,zn3,xn4,zn4,t1,ccdens,cnorm
      REAL*8 sum
!============
!     common/iter1/ rz(pncoil,20),zz(pncoil,20),bpolav(pncoil),
!    +     bpolr(5,pncoil,pncoil),bpolz(5,pncoil,pncoil),
!    1     sumr(5),sumz(5),ansr(5),ansz(5),
!    1     rfz(pncoil,5), zfz(pncoil,5),bpolmax(pncoil)
!     common/eqruns/findex
      dimension cval(10),div(10)
!     dimension tary(pnsave),xmary(pnsave),zmary(pnsave),
!    1     asfary(pnsave),abary(pnsave),aipary(pnsave),aefary(pnsave)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xmary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zmary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: asfary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: abary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aipary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aefary
!============      
      IF(.not.ALLOCATED(tary)) ALLOCATE( tary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(xmary)) ALLOCATE( xmary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(zmary)) ALLOCATE( zmary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(asfary)) ALLOCATE( asfary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(abary)) ALLOCATE( abary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(aipary)) ALLOCATE( aipary(pnsave), STAT=istat)
      IF(.not.ALLOCATED(aefary)) ALLOCATE( aefary(pnsave), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : sumplot  ' 
!============      
!
!
!     write(nterm,9331)
!9331 format("sumplot called")
!
!.....special added April 2000 to write new type 10 cards
      call newinput
!
      if(ifrst(4) .eq. 1) return
!
      xmin = ccon
      zmin = -alz
      xmax = alx
      zmax = alz
!.....external coils
!      if(itemp.eq.0) go to 761
      if(dxcoil(1).eq.0._R8) goto 761
      do 760 n=1,ncoil-nwire
      xmax = max(xmax,xcoil(n)+0.5_R8*dxcoil(n))
      xmin = min(xmin,xcoil(n)-0.5_R8*dxcoil(n))
      zmax = max(zmax,zcoil(n)+0.5_R8*dzcoil(n))
      zmin = min(zmin,-zmax)
  760 continue
  761 continue
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
!
      nplrecm = nplrec
      rewind nsc2
      irec = pnx*pnz
      do 100 npl=1,nplrecm
      call bufin(nsc2,vecx(1,1),vecx(pnx,pnz))
      icount = icount+irec
      itmax = 1+isym
      do 100 it=1,itmax
      sym = 1._R8
      if(it.eq.2) sym = -1._R8
      do 50 j=2,nzp
   50 vzy(j) = sym*zary(j)
!
!.....define k1,k2,cval to plot psilim contour only
      k1 = 1
      k2 = 1
      cval(1) = vecx(1,1)
      call rcontr(k1,cval,k2,vecx,penx,xary,iminy(npl),imaxy(npl),1,     &  
     &                         vzy,jminy(npl),jmaxy(npl),1)
  100 continue
!      if(itemp.eq.0) go to 790
      if(dxcoil(1).eq.0._R8) goto 790
!
!.....draw and label external coils
      do 780 n=1,ncoil-nwire
      if(xcoil(n).ge.100) go to 780
      if(acoef(901).eq.0._R8) then
      dxc=dxcoil(n)
      dzc=dzcoil(n)
      endif
!
      if(acoef(901).gt.0._R8) then
         if(fcu(n).gt.0._R8) then
! assume fcu(ico) kA/cm^2  current density sizing
            if(fss(n).eq.0._R8) then
      print *,'***warning***fss on type39 should be coil aspect ratio'
! assume fss(ic0) coil aspect ratio for sizing
            endif
            dxc=sqrt(abs(ccoil(n))*udsi*1.E-3_R8/(fcu(n)*fss(n)))/       &  
     & 100._R8
            dzc=dxc*fss(n)
         endif
         if(fcu(n).le.0._R8) then
            dxc=dxcoil(n)
            dzc=dzcoil(n)
         endif
      endif
      x1 = xcoil(n)-.5_R8*dxc
      z1 = zcoil(n)-.5_R8*dzc
      x2 = x1+dxc
      z2 = z1
      x3 = x2
      z3 = z2+dzc
      x4 = x1
      z4 = z3
      call setcrt(x1,z1)
      call vector(x2,z2)
      call vector(x3,z3)
      call vector(x4,z4)
      call vector(x1,z1)
      if(isym.eq.0) go to 780
      call setcrt(x1,-z1)
      call vector(x2,-z2)
      call vector(x3,-z3)
      call vector(x4,-z4)
      call vector(x1,-z1)
      call setlch(x2,z1,0,1,0,-1)
      write(s100,898) igroupc(n)
 898  format(i2)
      call gtext(s100,80,0)
  780 continue
!
  790 continue
      do 81 ii=1,nwire
   81 call coildr(iwire(ii),jwire(ii),ilwire(ii))
      print *,'iplate,nplate',iplate,nplate
      if(iplate.eq.0 .or. nplate.eq.0) go to 8100
      do 8110 n=1,nplate
      call setcrt(xsega(n,1),zsega(n,1))
      do 8111 l=2,nseg(n)+1
 8111 call vector(xsega(n,l),zsega(n,l))
      if(isym.eq.0) go to 8110
      z1neg = -zsega(n,1)
      call setcrt(xsega(n,1),z1neg)
      do 8112 l=2,nseg(n)+1
      zlneg = -zsega(n,l)
 8112 call vector(xsega(n,l),zlneg)
 8110 continue
 8100 continue
!
      if(acoef(901).gt.0._R8) then
      do n=1,ncnt
      xp1=xcon0(istart,n)-deex/2._R8
      zp1=zcon0(istart,n)-deez/2._R8
      call setcrt(xp1,zp1)
      xp2=xp1+deex
      zp2=zp1
      call vector(xp2,zp2)
      xp3=xp2
      zp3=zp2+deez
      call vector(xp3,zp3)
      xp4=xp3-deex
      zp4=zp3
      call vector(xp4,zp4)
      xp5=xp4
      zp5=zp4-deez
      call vector(xp5,zp5)
      enddo
      if(acoef(901).eq.3._R8.or.acoef(901).eq.4._R8) then
      xn1=acoef(903)-deex/2._R8
      zn1=acoef(904)-deez/2._R8
      call setcrt(xn1,zn1)
      xn2=xn1+deex
      zn2=zn1+deez
      call vector(xn2,zn2)
      xn3=xn2-deex
      zn3=zn2
      call setcrt(xn3,zn3)
      xn4=xn3+deex
      zn4=zn3-deez
      call vector(xn4,zn4)
      endif
      endif
!
      call frscj(9)
      if(acoef(901).gt.0._R8) then
! calculate poloidal magnetic field at coils due to all other coils
      call iterbp(1)
!
      call map(0.0_R8,1._R8,0.0_R8,1._R8,0._R8,1._R8,0._R8,1._R8)
! write coil current values
      call setlch(0._R8,1._R8,0,1,0,-1)
 
      write(s100,3601)
 3601 format(" grp",2x,"xc[m]",2x,"zc[m]",2x,"dx[m]",2x,"dz[m]",3x,      &  
     &"Wt",5x,"I[kA]",2x,                                                &  
     &     "J[kA/cm^2]",1x,"Bmax[T]",2x,"Bav[T]")
      write(nterm,3601)
      call gtext(s100,80,0)
      do 86 n=nc0,ncoil-nwire
 
      t1 = ccoil(n)*udsi*1.E-3_R8
      ccdens=abs(t1)*1.E-4_R8/(dxcoil(n)*dzcoil(n))
      write(s100,3602) n,xcoil(n),zcoil(n),dxcoil(n),dzcoil(n),          &  
     &gcur(2,n)*1.E-3_R8,t1,                                             &  
     &   ccdens,bpolmax(n),bpolav(n)
 3602 format(i3,5f7.2,f10.2,1x,f7.2,2x,f7.2,2x,f7.2)
      write(nterm,3602) n,xcoil(n),zcoil(n),dxcoil(n),dzcoil(n),         &  
     &gcur(2,n)*1.E-3_R8,t1,                                             &  
     &   ccdens,bpolmax(n),bpolav(n)
      call gtext(s100,80,0)
 86   continue
! calculate Cnorm
      sum=0._R8
      do n=nc0,ncoil-nwire
      t1 = ccoil(n)*udsi*1.E-3_R8
      sum=sum+t1*t1
      enddo
      cnorm=sqrt(sum)
      write(s100,3603) cnorm
 3603 format(" norm of coil current vector [kA]=",f12.2)
      call gtext(s100,80,0)
      write(s100,3171) findex
 3171 format(" nindx=   ",f7.3)
      call gtext(s100,80,0)
!
      call frscj(17)
      endif
!
      if (idiv.ne.0)    call  sumplotros
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
