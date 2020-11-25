      subroutine stepon
!......2.10    stepon
!
!.....advance solution 1 step in time
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE CLINAM
      USE SCR1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ifirstso,i,j,nzp1,ip
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 oi,oj
!============
      data ifirstso/0/
!
      if(lrswtch .eq. 0) go to 4
      do 3 i=2,nxp
      do 3 j=2,nzp
      utsave(i,j)=0._R8
      gtsave(i,j)=0._R8
      ptsave(i,j)=0._R8
      u(i,j) = 0._R8
      uo(i,j)=0._R8
      b(i,j) = 0._R8
      bo(i,j) = 0._R8
      w(i,j) = 0._R8
      wo(i,j) = 0._R8
    3 continue
    4 continue
      if(ekino .gt. ekin.and. acoef(896).gt.0) then
      write(nout,9111) kcycle,times,ekin,ekino
      write(nterm,9111) kcycle,times,ekin,ekino
 9111 format(" velocity zeroed:  kcycle=",i7, "time=",1pe16.8,           &  
     & "ekin,ekino=",     1p2e12.4)
      ekin = 0._R8
      ekino=0._R8
      if(kcycle.eq.0) ekino = 0._R8
      do 801 i=2,nxp
      do 801 j=2,nzp
      u(i,j) = 0._R8
      uo(i,j)=0._R8
      b(i,j) = 0._R8
      bo(i,j) = 0._R8
      w(i,j) = 0._R8
      wo(i,j) = 0._R8
  801 continue
      endif
!
      nzp1 = nzp + 1
!
      do 5 j=1,nzp1
      piz(j) = psi(1,j)
      wiz(j) = w(1,j)
      giz(j) = g(1,j)
      uiz(j) = u(1,j)
      aiz(j) = abig(1,j)
      biz(j) = b(1,j)
      riz(j) = r(1,j)
      qiz(j) = q(1,j)
!
      pip(j) = psi(2,j)
      wip(j) = w(2,j)
      gip(j) = g(2,j)
      uip(j) = u(2,j)
      abip(j) = abig(2,j)
      bip(j) = b(2,j)
      rip(j) = r(2,j)
      qip(j) = q(2,j)
    5 continue
!
!.....set left boundary fluxes for magnetic axis
!
      do 6 j=2,nzp1
      wf1pa(j) = 0._R8
      gf1pa(j) = 0._R8
      uf1pa(j) = 0._R8
      bf1pa(j) = 0._R8
      rf1pa(j) = 0._R8
      qf1pa(j) = 0._R8
    6 continue
!     define velocity arrays vecx and vecz
      do 1006 i=2,nxp
      do 1006 j=2,nzp
      oi = (omeg(i+1,j+1)+omeg(i+1,j)-omeg(i,j+1)-omeg(i,j))
      oj = (omeg(i+1,j+1)+omeg(i,j+1)-omeg(i,j)-omeg(i+1,j))
      vecx(i,j) = (abig(i,j+1)-abig(i,j-1))/xary(i)*.5_R8/deez+oi*.5_R8/  &  
     & deex
      vecz(i,j) =-(abig(i+1,j)-abig(i-1,j))/xary(i)*.5_R8/deex+oj*.5_R8/  &  
     & deez
 1006 continue
!
!...define and save bo2
      do 1007 i=2,nxp
      do 1007 j=1,nzp
      bo2(i,j) = bo(i,j)
 1007 continue
!
!.....start loop on constant i lines
!
      do 200 i=2,nxp
      icol  = i
      ip = i+1
!
      do 10 j=1,nzp1
      pim(j) = piz(j)
      wim(j) = wiz(j)
      gim(j) = giz(j)
      aim(j) = aiz(j)
      uim(j) = uiz(j)
      bim(j) = biz(j)
      rim(j) = riz(j)
      qim(j) = qiz(j)
   10 continue
!
      do 11 j=1,nzp1
      piz(j) = pip(j)
      wiz(j) = wip(j)
      giz(j) = gip(j)
      aiz(j) = abip(j)
      uiz(j) = uip(j)
      biz(j) = bip(j)
      riz(j) = rip(j)
      qiz(j) = qip(j)
   11 continue
!
      do 12 j=1,nzp1
      pip(j) = psi(ip,j)
      wip(j) = w(ip,j)
      gip(j) = g(ip,j)
      uip(j) = u(ip,j)
      abip(j) = abig(ip,j)
      bip(j) = b(ip,j)
      rip(j) = r(ip,j)
      qip(j) = q(ip,j)
   12 continue
!
!
!.....zone centered variables
!
      if(lrswtch.ge.1) go to 201
      call advanc
!
!.....point centered variables
!
  201 call advanp
  200 continue
!
!.....hyper resistivity contribution to psi and g
      if(times.gt.acoef(57) .and. times.lt.acoef(58)                     &  
     &   .and. isaw.lt.2) hypermult = acoef(64)*udsi
      if(  (hypermult.gt.0 .and. lrswtch.eq.0 .and.                      &  
     &     ifirstso.ge.1 .and. isaw   .ge. 2       )                     &  
     &       .or.    (  irfp.eq.1                  )                     &  
     &       .or.    (  hypermult.gt.0 .and.                             &  
     & times.gt.acoef(57) .and. times.lt.acoef(58) ) ) then
!
      call hyper
!
!.....hyper called
      write(nterm,6656) kcycle,nloop, delpmx, delgmx,hypermult
 6656 format("Hyper called..kcycle,nloop,delpmx,delgmx=",i7,i4,1p3e10.2)    
!
      if(iplt.gt.nskipl.and.(imovie.eq.0 .or. imovie.ge.10)              &  
     &     .and.irfp.eq.1) then
                                         call lingin1
                                         call lingin2
                                         call lingin3
                                         endif
 
       endif
      ifirstso = 1
!
!
!.....substepping for toroidal field , vel divergence, and poloidal flux
      call advan23
!
!
!.....define auxialliary arrays
      call auxdef
!
!.....extrapolate values at outer boundary
!
      call boundc
!
      if(acoef(1).gt. 0 .and. acoef(1) .le. 4._R8) go to 400
!
!.....invert del-star operator
      if(lrswtch.ge.1) return
      call elliptic(2,abig)
!
!.....invert del-sqrared operator
      call elliptic(3,omeg)
      return
  400 continue
      if(lrswtch.ge.1) return
      call iterate(2,abig)
      call iterate(3,omeg)
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
