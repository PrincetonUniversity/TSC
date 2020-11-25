!#include "f77_dcomplx.h"
      subroutine facr(iop,ibnd,xmin,xmax,zmin,zmax,nx,nz,psi,isym,       &  
     &    itype,inegp,dert,derb,derl,derr,nout)
!.....4.11 facr
      USE PARAM
      USE SCR8  
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!*****************************************************************
!                                                                *
!                                                                *
!                        elliptic solver                         *
!                        ===============                         *
!                                                                *
!.....this subroutine solves various elliptic equations on an    *
!.....equally spaced rectangular grid using the method of        *
!.....fourier analysis and tridiagonal matrix inversion.                     *
!...(see hockney...methods in computational physics,vol.9)       *
!                                                                *
!                                                                *
!...input:                                                       *
!...=====                                                        *
!                                                                *
!...iop  = 1...delstar   operator in cyclindrical coordinates    *
!....... = 2...laplacian operator in cyclindrical coordinates    *
!....... = 3...laplacian operator in cartesian    coordinates    *
!                                                                *
!...ibnd = 1...dirichlet boundary condition                      *
!...ibnd = 2...neuman    boundary condition  *
!                                                                *
!...xmin = x-position of left   grid boundary                    *
!...xmax = x-position of right  grid boundary                    *
!...zmin = z-position of bottom grid boundary                    *
!...zmax = z-position of top    grid boundary                    *
!                                                                *
!...nx   = number of x grid points                  *
!...nz   = number of z grid points                  *
!                                                                *
!...psi  = right side plus boundary cond's (ibnd=1)              *
!....... = right side                      (ibnd=2)              *
!     => important note:  rhs must be pre-multiplied by dx*dz
!                         bc must be pre-multiplied by dx/dz
!
!                                                                *
!...isym is 1 for symmetry about midplane
!           0 for no symmetry
!
!...itype = 1 for point centered with even symmetry
!           2 for point centered with odd symmetry
!           3 for cell centered with even symmetry
!
!...inegp is returned non-zero if error is detected
!
!...dert,derb,derl,derr are finite differences at boundary for ibnd=2
!      => note:  must be pre-multiplied by dx/dz
!
!                                                                *
!...output:                                                      *
!...======                                                       *
!                                                                *
!...psi  = solution times (dx/dz)                                              *
!                                                                *
!                                                                *
!...notes:                                                       *
!...=====                                                        *
!                                                                *
!...1) in the calling routine psi must be dimensioned            *
!...........psi(0:penx-1,0:penz-1).....where penx-1 and penz-1           *
!......are given in the parameter statement in this routine      *
!                                                                *
!                                                                *
!*****************************************************************
!
!
!.......................................................................
!...form grid dependent quantities
!.......................................................................
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ibnd,nx,nz,isym,itype,inegp,nout,iop,istrt,ni,nj,i,j
      INTEGER k,l,imax,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xmin,xmax,zmin,zmax,pi,dx,dz
      REAL*8 rfac,fj,factor,fac
      REAL*8 ta,tb,alam,bfcc
      REAL*8 cfcc,phitni,diff,del
      REAL*8 AREAL
!============
      dimension istrt(6)
      data istrt / 6*0 /
!============
      REAL*8, DIMENSION(0:penx) :: dert, derb
      REAL*8, DIMENSION(0:penz) :: derl, derr
      REAL*8, DIMENSION(0:penx-1,0:penz-1) :: psi
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: av, fv, gama, g
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: phit
!============      
      INTEGER :: istat = 0 
!============
!
      if (.NOT.ALLOCATED(av)) ALLOCATE (av(penx,penz), STAT=istat)
      if (.NOT.ALLOCATED(fv)) ALLOCATE (fv(penx,penz), STAT=istat)
      if (.NOT.ALLOCATED(gama)) ALLOCATE (gama(penx,penz), STAT=istat)
      if (.NOT.ALLOCATED(g)) ALLOCATE (g(penx,penz), STAT=istat)
      if (.NOT.ALLOCATED(phit)) ALLOCATE (phit(0:penx-1,0:penz-1),       &  
     &                                    STAT=istat)
      if (istat .ne. 0) stop "Allocation Error : facr"

      pi = 3.1415926535897_R8
!
      ni = nx - 1
      nj = nz - 1
      dx = (xmax-xmin)/ni
      dz = (zmax-zmin)/nj
      rfac = (dz/dx)**2
      fj = AREAL(nj)
!
      do 2 i = 0,ni
    2 x(i) = xmin + i*dx
      do 3 j = 0,nj
    3 z(j) = zmin + j*dz
!
!...trigonometric quantities
!
      if(istrt(itype) .gt. 0 ) go to 15
      istrt(itype) = 1
      if(isym.eq.1) go to 6
      if(itype.eq.3) go to 30
      do 5 k = 0,nj
      do 5 j = 0,nj
      sn(itype,k,j) = sin(pi*k*j/fj)
      cn(itype,k,j) = cos(pi*k*j/fj)
    5 continue
      go to 15
    6 continue
      go to(10,20,30),itype
   10 continue
      do 4 l=0,nj
      do 4 j=0,nj
      sn(itype,l,j) = sin(pi*(l+.5_R8)*(fj-j)/fj)
      cn(itype,l,j) = cos(pi*(l+.5_R8)*(fj-j)/fj)
    4 continue
      go to 15
   20 continue
      do 7 l=0,nj
      do 7 j=0,nj
      sn(itype,l,j) = sin(pi*l*(fj-j)/fj)
      cn(itype,l,j) = cos(pi*l*(fj-j)/fj)
    7 continue
      go to 15
   30 continue
      do 8 l=0,nj
      do 8 j=0,nj
      sn(itype,l,j) = sin(l*pi*(j+.5_R8)/(fj+1.0_R8))
    8 cn(itype,l,j) = cos(l*pi*(j+.5_R8)/(fj+1.0_R8))
   15 continue
!
!
!
!...operator dependent quantities
!
      if(istrt(3+iop) .gt. 0) go to 100
      istrt(3+iop) = 1
      do 85 i=0,ni
      go to(86,87,88),iop
   86 afac(iop,i) = -(1._R8/(1._R8+.5_R8*dx/x(i))+1._R8/(1._R8-.5_R8*dx/  &  
     & x(i)))
      bfac(iop,i) = 1._R8/(1._R8-.5_R8*dx/x(i))
      cfac(iop,i) = 1._R8/(1._R8+.5_R8*dx/x(i))
      go to 85
   87 afac(iop,i) = -2._R8
      bfac(iop,i) = 1._R8-.5_R8*dx/x(i)
      cfac(iop,i) = 1._R8+.5_R8*dx/x(i)
      go to 85
   88 afac(iop,i) = -2._R8
      bfac(iop,i) = 1._R8
      cfac(iop,i) = 1._R8
   85 continue
  100 factor = 0._R8
      if(iop.eq.1) factor = -.5_R8*dx
      if(iop.eq.2) factor =  .5_R8*dx
      if(ibnd.eq.2) go to 200
!.......................................................................
!...dirichlet boundary conditions
!.......................................................................
!
!...form abnd and bbnd vectors
!
      if(isym.ne.0 .and. itype.eq.1) go to 106
      do 105 i = 0,ni
      abnd(i) =  psi(i,0)
      bbnd(i) = (psi(i,nj) - psi(i,0))/fj
  105 continue
      go to 107
  106 continue
      do 108 i=0,ni
      abnd(i) = psi(i,nj)
      bbnd(i) = 0._R8
  108 continue
  107 continue
!
!...modify right side to give homogenious top and bottom
!...boundary conditions
!
      do 110 i = 1,ni-1
      fac = factor/x(i)
      ta  = bfac(iop,i)*abnd(i-1)+afac(iop,i)*abnd(i)                    &  
     &     +cfac(iop,i)*abnd(i+1)
      tb  = bfac(iop,i)*bbnd(i-1)+afac(iop,i)*bbnd(i)                    &  
     &     +cfac(iop,i)*bbnd(i+1)
      do 110 j=1-isym,nj-1
      psi(i,j) = psi(i,j) - rfac*(ta + j*tb)
  110 continue
!
!...modify left and right boundary conditions
!
      do 115 j = 0,nj
      psi(0 ,j) = psi(0 ,j) - abnd(0 ) - j*bbnd(0 )
      psi(ni,j) = psi(ni,j) - abnd(ni) - j*bbnd(ni)
  115 continue
!
!...transform right side and boundary conditions at i = 0,ni
!
      do 125 k = 0,nj-1
      if(isym.eq.1 .and. itype.eq.1) go to 116
      do 119 i=0,ni
      phit(i,k) = 0._R8
  119 continue
      go to 117
  116 continue
      do 118 i=0,ni
  118 phit(i,k) = (-1._R8)**k*psi(i,0)/2._R8
  117 continue
      do 120 j = 1,nj-1
      do 120 i = 0,ni
      phit(i,k) = phit(i,k) + psi(i,j)*sn(itype,k,j)
  120 continue
      do 125 i = 0,ni
      phit(i,k) = 2._R8*phit(i,k)/fj
  125 continue
!
!...form off diagonal elements
!
      do 130 i = 1,ni-1
      fac = factor/x(i)
      bv(i) = rfac*bfac(iop,i)
      cv(i) = rfac*cfac(iop,i)
  130 continue
!
!...invert transformed flux
!
      do 140 k = 0,nj-1
      i = 1
      fv(i,k+1) = phit(i,k) - bv(i)*phit(i-1,k)
      do 135 i = 2,ni-2
      fv(i,k+1) = phit(i,k)
  135 continue
      i = ni-1
      fv(i,k+1) = phit(i,k) - cv(i)*phit(i+1,k)
      do 140 i = 1,ni-1
      av(i,k+1) = 2._R8*(cn(itype,k,1+isym*(nj-2))-1._R8)+rfac*afac(iop,  &  
     & i)
  140 continue
      imax = ni-1
!
!
!
!...invert the block tridiagnal system
!
!
!...form gama and g vectors
!
      do 499 k = 1,nj
  499 g(1,k)  = fv(1,k)/av(1,k)
      do 500 i = 2,imax
!cj dir$ ivdep
      do 500 k = 1,nj
      gama(i-1,k) = cv(i-1)/av(i-1,k)
      av(i,k)     = av(i,k) - bv(i)*gama(i-1,k)
      g(i,k)      = (fv(i,k) - bv(i)*g(i-1,k))/av(i,k)
  500 continue
!
!...back substitute to obtain x
!
!cj dir$ ivdep
      do 145 k = 1,nj
  145 phit(imax,k-1) = g(imax,k)
      do 150 ii = 1,imax-1
      i = imax - ii
!cj dir$ ivdep
      do 150 k = 1,nj
      phit(i,k-1) = g(i,k) - gama(i,k)*phit(i+1,k-1)
  150 continue
!
!...real space solution
!
      do 160 j = 1-isym,nj-1
      do 158 i = 1,ni-1
  158 psi(i,j) =  abnd(i) + j*bbnd(i)
      do 155 k = 0,nj-1
      do 155 i = 1,ni-1
      psi(i,j) = psi(i,j) + phit(i,k)*sn(itype,k,j)
  155 continue
  160 continue
!
      go to 400
!.......................................................................
!...neuman boundary conditions
!.......................................................................
  200 continue
!.
!.....form abnd and bbnd vectors
      do 205 i=0,ni
      if(isym.ne.0) derb(i) = 0
      bbnd(i) = (dert(i)-derb(i))/(fj+1)
      abnd(i) = derb(i) + 0.5_R8*bbnd(i)
  205 continue
!
!.....modify right side to give homogeneous boundary conditions
!     on top and bottom
!
      abnd(-1) = 0
      bbnd(-1) = 0
      abnd(ni+1) = 0
      bbnd(ni+1) = 0
      do 210 i=0,ni
      ta=bfac(iop,i)*abnd(i-1)+afac(iop,i)*abnd(i)+cfac(iop,i)*abnd(i+1)     
       
      tb=bfac(iop,i)*bbnd(i-1)+afac(iop,i)*bbnd(i)+cfac(iop,i)*bbnd(i+1)     
       
      do 210 j=0,nj
  210 psi(i,j) = psi(i,j) - rfac*(ta*j+0.5_R8*tb*j**2) - bbnd(i)
!
!.....modify left and right boundary conditions
!
      do 215 j=0,nj
      derl(j) = derl(j) - j*abnd(0) - 0.5_R8*j**2*bbnd(0)
      derr(j) = derr(j) + j*abnd(ni)+ 0.5_R8*j**2*bbnd(ni)
  215 continue
!
!...transform right side and boundary conditions
!
      do 225 k = 0,nj
      derlt(k) = 0._R8
      derrt(k) = 0._R8
      do 219 i = 0,ni
      phit(i,k) = 0._R8
  219 continue
      do 220 j = 0,nj
      derlt(k) = derlt(k) + derl(j)*cn(itype,k,j)
      derrt(k) = derrt(k) + derr(j)*cn(itype,k,j)
      do 220 i = 0,ni
      phit(i,k) = phit(i,k) + psi(i,j)*cn(itype,k,j)
  220 continue
      derlt(k) = 2._R8*derlt(k)/(fj+1)
      derrt(k) = 2._R8*derrt(k)/(fj+1)
      do 225 i = 0,ni
      phit(i,k) = 2._R8*phit(i,k)/(fj+1.0_R8)
  225 continue
!
!...calculate the transformed flux
!
!........special treatment for k=0
      k = 0
!.....temporary storage
      do 230 i = 0,ni
      psi(i,k) = phit(i,k)
  230 continue
!
      alam = -2._R8*rfac
      phit(0,k) = 0._R8
      bfcc = rfac*(1._R8-factor/x(0))
      cfcc = rfac*(1._R8+factor/x(0))
      phit(1,k) = phit(0,k) + (psi(0,0)+bfcc*derlt(0))/cfcc
      do 235 i = 2,ni
      bfcc = rfac*(1._R8- factor/x(i-1))
      cfcc = rfac*(1._R8+ factor/x(i-1))
      phit(i,k) = (psi(i-1,k)-alam*phit(i-1,k)-bfcc*phit(i-2,k))/cfcc
  235 continue
      bfcc = rfac*(1._R8-factor/x(ni))
      cfcc = rfac*(1._R8+factor/x(ni))
      phitni = phit(ni-1,k)+(-psi(ni,k)+cfcc*derrt(0))/bfcc
      diff = abs(phitni-phit(ni,k))
      if(diff.le.1.E-6_R8*abs(phitni).or.diff.le.1.E-12_R8) go to 238
      write(nout,1000) phitni,phit(ni,k)
 1000 format(' incompatible bc',1p2e12.4)
      inegp = 19
  238 continue
!........tridiagonal for k = 1,2,...,nj
      cv(1) = rfac*(1._R8+factor/x(0))
      do 245 i = 2,ni
      del   = factor/x(i-1)
      bv(i) = rfac*(1._R8-del)
      cv(i) = rfac*(1._R8+del)
  245 continue
      bv(ni+1) = rfac*(1._R8-factor/x(ni))
!
      do 241 k = 1,nj
      alam = 2._R8*(cos(pi*k/(fj+1))    - 1._R8- rfac)
      do 240 i = 1,ni+1
      av(i,k) = alam
      fv(i,k) = phit(i-1,k)
  240 continue
      fv(1,k) = fv(1,k) + rfac*(1._R8-factor/x(0))*derlt(k)
      fv(ni+1,k) = fv(ni+1,k) - rfac*(1._R8+factor/x(ni))*derrt(k)
      av(1,k) = -cv(1)+2._R8*(cos(pi*k/(fj+1)) - 1._R8)
  241 av(ni+1,k) = -bv(ni+1) + 2._R8*(cos(pi*k/(fj+1._R8))-1._R8)
      imax = ni + 1
!
!...invert the block tridiagnal system
!
!
!...form gama and g vectors
!
      do 509 k = 1,nj
  509 g(1,k)  = fv(1,k)/av(1,k)
      do 510 i = 2,imax
!cj dir$ ivdep
      do 510 k = 1,nj
      gama(i-1,k) = cv(i-1)/av(i-1,k)
      av(i,k)     = av(i,k) - bv(i)*gama(i-1,k)
      g(i,k)      = (fv(i,k) - bv(i)*g(i-1,k))/av(i,k)
  510 continue
!
!...back substitute to obtain x
!
      do 249 k = 1,nj
  249 phit(imax-1,k) = g(imax,k)
      do 250 ii = 1,imax-1
      i = imax - ii
!cj dir$ ivdep
      do 250 k = 1,nj
      phit(i-1,k) = g(i,k) - gama(i,k)*phit(i,k)
  250 continue
!
!...real space solution
!
      do 260 j = 0,nj
!
!.....restore boundary conditions
      derl(j) = derl(j) + j*abnd(0) + 0.5_R8*j**2*bbnd(0)
      derr(j) = derr(j) - j*abnd(ni)- 0.5_R8*j**2*bbnd(ni)
      do 258 i = 0,ni
  258 psi(i,j) = .5_R8* phit(i,0) + j*abnd(i) + 0.5_R8*j**2*bbnd(i)
      do 255 k = 1,nj
      do 255 i = 0,ni
      psi(i,j) = psi(i,j) + phit(i,k)*cn(itype,k,j)
  255 continue
  260 continue
!
!
  400 continue

!     DEALLOCATE (phit, av, fv, gama,g)

      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
