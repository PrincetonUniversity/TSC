      subroutine sumplotros
!           Mod 30 Dec 94 for color.   ROS
!           Mod 22 Dec 96 for call from scj sumplot.   ROS
!.....produce a summary plot in a separate file 'sumgmeta'
!
      USE CLINAM
      USE SCRATCH
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ipass,n,nplrecm,irec,npl,icount,itmax,it,j,k1,k2,ii,l
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cval,div
      REAL*8 xmin,zmin,xmax,zmax,xave,xleft,xright,sym,x1,z1,x2,z2
      REAL*8 x3,z3,x4,z4,z1neg,zlneg
!============
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
      if (istat .ne. 0) stop 'Allocation Error : sumplotros  ' 
!============      
!
!ccccccc    write (nterm,*)  'entering sumplotros'
         call plote
         if ( numargs .lt. 1 ) then
            filename = 'sumgmeta'
         else
            filename = 'sumgmeta' // '.' // trim(suffix)
         end if
         call ncarcgm (1, trim(filename))
!
      do 9000  ipass=2,2
!
      xmin = ccon
      zmin = -alz
      xmax = alx
      zmax = alz
!.....external coils
      if(itemp.eq.0) go to 761
      do 760 n=1,ncoil-nwire
      xmax = max(xmax,xcoil(n)+0.5_R8*dxcoil(n))
      xmin = min(xmin,xcoil(n)-0.5_R8*dxcoil(n))
      zmax = max(zmax,zcoil(n)+0.5_R8*dzcoil(n))
      zmin = min(zmin,-zmax)
  760 continue
  761 continue
      xave = .5_R8*(xmax+xmin)
      if (ipass.eq.1)   then
         if((zmax-zmin).gt.xmax-xmin) xmax = xave+.5_R8*(zmax-zmin)
         if((zmax-zmin).gt.xmax-xmin) xmin = xave-.5_R8*(zmax-zmin)
         if((xmax-xmin).gt.(zmax-zmin)) zmax = .5_R8*(xmax-xmin)
         if((xmax-xmin).gt.(zmax-zmin)) zmin = -.5_R8*(xmax-xmin)
         call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.285_R8,1._R8)
      else
!ccccccc    xmax = xave+.5*(zmax-zmin)
!ccccccc    xmin = xave-.5*(zmax-zmin)
      xleft  = 0.5_R8- 0.4_R8* (xmax-xmin)/(zmax-zmin)
      xright = 0.5_R8+ 0.4_R8* (xmax-xmin)/(zmax-zmin)
      if (xleft .lt.0.1_R8)   xleft  = 0.1_R8
      if (xright.gt.0.9_R8)   xright = 0.9_R8
      call maps(xmin,xmax,zmin,zmax,xleft,xright,.15_R8,.95_R8)
            endif
!
      nplrecm = nplrec
      rewind nsc2
      irec = pnx*pnz
      do 100 npl=1,nplrecm
      if (npl.eq.1)   call colora("blue")
      if (npl.eq.2)   call colora("cyan")
      if (npl.eq.3)   call colora("green")
      if (npl.eq.4)   call colora("yellow")
      if (npl.eq.5)   call colora("red")
      if (npl.eq.6)   call colora("white")
      if (npl.eq.7)   call colora("blue")
      if (npl.eq.8)   call colora("cyan")
      if (npl.eq.9)   call colora("green")
      if (npl.eq.10)   call colora("yellow")
      if (npl.eq.11)   call colora("red")
      call bufin(nsc2,vecx(1,1),vecx(pnx,pnz))
!ccccccc      buffer in(nsc2,1) (vecx(1,1),vecx(pnx,pnz))
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
      call colora("white")
      if(itemp.eq.0) go to 790
!
!.....draw external coils
      do 780 n=1,ncoil-nwire
      if(xcoil(n).ge.100) go to 780
      x1 = xcoil(n)-.5_R8*dxcoil(n)
      z1 = zcoil(n)-.5_R8*dzcoil(n)
      x2 = x1+dxcoil(n)
      z2 = z1
      x3 = x2
      z3 = z2+dzcoil(n)
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
  780 continue
  790 continue
      do 81 ii=1,nwire
   81 call coildr(iwire(ii),jwire(ii),ilwire(ii))
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
      if (ipass.eq.1)   then
            call frscj(9)
!
      else
!ccccccc    write (nterm,*)  'start of id write'
!.....write identifier
      call setld(77._R8,5._R8,1,0,1,0)
!ccccc      write(s100,8485) nframe
!ccccc      call gtext(s100,80,0)
 8485 format(1x,i8)
      call setld(2._R8,48._R8,1,0,1,0)
      write(s100,8484) name
      call gtext(s100,80,0)
      call setld(65._R8,1._R8,1,0,1,0)
      write(s100,8484) idate
      call gtext(s100,80,0)
      call setld(74._R8,1._R8,1,0,1,0)
      write(s100,8484) itime
      call gtext(s100,80,0)
 8484 format(8a10 )
!ccccccc 8484 format(1x,10a8 )
!
      call frame(0)
!ccccccc    write (nterm,*)  'frame(0) called'
            endif
!                    Deactivate workstation 1
!ccccccc    call  gdawk (1)
!                    Close workstation 1
!ccccccc    call  gclwk (1)
 9000   continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
