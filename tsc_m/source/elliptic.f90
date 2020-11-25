      subroutine elliptic(itype,aa1)
!.......................................................................
!
!.....sets up arrays and calls facr to invert elliptic operators
!
!....... itype=1 for poloidal flux solution
!....... itype=2 for stream function solution
!....... itype=3 for velocity potential solution
!....... itype=4 for Von-Hagenow boundary conditions setup
!
!.......................................................................
!
      USE CLINAM
      USE SCR7, except_this => aa1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,imin,jmin,imax,jmax,itypesw,i,j,ii,jj,iop,ibnd
      INTEGER nxx,nzz
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dxz,dzx,xmin,xmax,zmin,zmax
      REAL*8 bcsum,usum,csum
      REAL*8, DIMENSION(penx,penz) :: aa1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pass
      INTEGER :: istat
!============
      imin=3
      jmin=3
      imax = nx
      jmax = nz
!
      dxz = deex/deez
      dzx = deez/deex
!
!......define boundary conditions
!
!
      itypesw = itype
      if(itype.eq.4) itypesw = 1
      if(itype .eq. 4) go to 11
      if(itype .eq. 3) go to 20
      call bounda(itype)
      if(lrswtch.ge.1 .and. itype.eq.2) return
!.......................................................................
!....... * * * delstar * * *
!.......................................................................
      go to 14
   11 continue
!
!.....special boundary conditions to implement Von-Hagenow method
      do 12 i=imin,imax
      psibb(i)=0._R8
      psibt(i)=0._R8
   12 continue
      do 13 j=jmin-1,jmax+1
      psibl(j)=0._R8
      psibr(j)=0._R8
   13 continue
   14 continue

      if (.NOT.ALLOCATED(pass)) ALLOCATE(pass(0:penx-1,0:penz-1),        &  
     &                                  STAT=istat)
      if (istat .ne. 0) stop "Allocation Error : elliptic" 

      if(itype.ne.2) go to 50
      do 39 i=3,nx
      do 39 j=3-isym,nz
      vd(i,j) = b(i,j)*xary(i)
   39 continue
   50 continue
      do 5 i = imin,imax
      ii = i - 2
      if(isym.eq.0) pass(ii,0)      = dxz*psibb(i)
      if(isym.ne.0) pass(ii,0) = vd(i,2)*(2._R8-itypesw)
      pass(ii,jmax-1) = dxz*psibt(i)
      aa1(i,nzp) = psibt(i)
      if(isym.eq.0) aa1(i,2) = psibb(i)
      do 5 j = jmin,jmax
      jj = j - 2
      pass(ii,jj) =  vd(i,j)
    5 continue
!
!cj dir$ ivdep
      do 10 j = jmin-1,jmax+1
      jj = j - 2
      pass(0     ,jj) = dxz*psibl(j)
      pass(imax-1,jj) = dxz*psibr(j)
      aa1(2,j) = psibl(j)
      aa1(nxp,j) = psibr(j)
   10 continue
!
      iop  = 1
      ibnd = 1
      xmin = xary(2)
      xmax = xary(imax+1)
      zmin = zary(2)
      zmax = zary(jmax+1)
!
!cj      write(91,1089) iop,ibnd,nx,nz,isym,itypesw,ineg,nout
!cj      write(91,1090) xmin,xmax,zmin,zmax
!cj      write(91,1090) dert
!cj      write(91,1090) derb
!cj      write(91,1090) derl
!cj      write(91,1090) derr
!cj      write(91,1091)
!cj      write(91,1092) pass
!
!cj      write(*,*) " -- penx = ", penx, pnx
!cj      write(*,*) " -- penz = ", penz, pnz
!
!cj      write(91) iop,ibnd,nx,nz,isym,itypesw,ineg,nout
!cj      write(91) xmin,xmax,zmin,zmax
!cj      write(91) dert
!cj      write(91) derb
!cj      write(91) derl
!cj      write(91) derr
!cj      write(91)
!cj      write(91) pass
!  
      call facr(iop,ibnd,xmin,xmax,zmin,zmax,nx,nz,pass,isym,            &  
     &    itypesw,ineg,dert,derb,derl,derr,nout)
!
!cj      write(92,1089) iop,ibnd,nx,nz,isym,itypesw,ineg,nout
!cj      write(92,1090) xmin,xmax,zmin,zmax
!cj      write(92,1090) dert
!cj      write(92,1090) derb
!cj      write(92,1090) derl
!cj      write(92,1090) derr
!cj      write(92,1091)
!cj      write(92,1092) pass 
!cj      stop
!cj 1089 format(10i5)
!cj 1090 format(1p10e12.4)
!cj 1091 format(" pass")
!cj 1092 format(1p10e12.5)
!
      do 16 i = imin,imax
      ii = i - 2
      if(isym.ne.0) aa1(i,2) = dzx*pass(ii,0)
      do 15 j = jmin,jmax
      jj = j - 2
      aa1(i,j) = dzx*pass(ii,jj)
   15 continue
      if(isym.ne.0) aa1(i,1) = aa1(i,3)*(3-2*itypesw)
   16 continue
      do 17 i=2,nxp
      aa1(i,nz+2) = 2._R8*aa1(i,nz+1)-aa1(i,nz)
      if(isym.eq.0) aa1(i,1) = 2._R8*aa1(i,2) - aa1(i,3)
   17 continue
      do 18 j=2,nzp
      aa1(nx+2,j) = 2._R8*aa1(nx+1,j) - aa1(nx,j)
   18 aa1(   1,j) = 2._R8*aa1(   2,j) - aa1( 3,j)
!
      go to 40
!
!
   20 continue
!.......................................................................
!....... * * * laplacian * * *
!.......................................................................
!.....force u array to be compatible with neuman boundary conditions
!
      if (.NOT.ALLOCATED(pass)) ALLOCATE(pass(0:penx-1,0:penz-1),        &  
     &                                  STAT=istat)
      if (istat .ne. 0) stop "Allocation Error in elliptic" 


      bcsum = 0._R8
      do 331 j=3,nzp
  331 bcsum = bcsum + dzx*(xary(nxp)*omderr(j)-xary(2)*omderl(j))
      do 332 i=3,nxp
  332 bcsum = bcsum + dxz*xarh(i)*(omdert(i)-omderb(i))
!
      usum = 0._R8
      do 335 i=3,nxp
      do 335 j=3,nzp
      usum = usum + u(i,j)
  335 continue
!
      csum = 0
      if(acoef(20).gt.0) go to 350
      do 336 i=3,nxp
  336 csum = csum + (2-isym)*xarh(i)
      do 337 j=4,nz
  337 csum = csum + xarh(3) + xarh(nxp)
      ucor = (usum-bcsum)/csum
      do 340 i=3,nxp
  340 u(i,nzp) = u(i,nzp) - xarh(i)*ucor
      if(isym.eq.1) go to 343
      do 342 i=3,nxp
  342 u(i,3) = u(i,3) - xarh(i)*ucor
      go to 344
  343 continue
      do 345 i=3,nxp
  345 u(i,2) = u(i,3)
  344 continue
      do 341 j=4,nz
      u(3,j) = u(3,j) - xarh(3)*ucor
  341 u(nxp,j) = u(nxp,j) - xarh(nxp)*ucor
      go to 360
!
!.....spread correction over domain for acoef(20).gt.0
  350 continue
      do 351 i=3,nxp
      do 351 j=3,nzp
  351 csum = csum+xarh(i)
      ucor = (usum-bcsum)/csum
      do 353 i=3,nxp
      do 352 j=3,nzp
  352 u(i,j) = u(i,j) - xarh(i)*ucor
      if(isym.eq.1) u(i,2) = u(i,3)
  353 continue
!
  360 continue
!
      do 25 i = imin,imax+1
      ii = i - 3
      do 25 j = jmin,jmax+1
      jj = j - 3
      pass(ii,jj) =  u(i,j)/xarh(i)
   25 continue
!
      do 26 i=imin,imax+1
      ii = i - 3
      dert(ii) = omdert(i)*dxz
   26 derb(ii) = omderb(i)*dxz
      do 27 j=jmin,jmax+1
      jj = j-3
      derl(jj) = omderl(j)*dxz
   27 derr(jj) = omderr(j)*dxz
!
      iop  = 2
      ibnd = 2
      xmin = xary(imin) - .5_R8*deex
      xmax = xary(imax) + .5_R8*deex
      zmin = zary(jmin) - .5_R8*deez
      zmax = zary(jmax) + .5_R8*deez
      nxx  = nx - 1
      nzz  = nz - 1
!
      call facr(iop,ibnd,xmin,xmax,zmin,zmax,nxx,nzz,pass,isym,          &  
     &    itypesw,ineg,dert,derb,derl,derr,nout)
!
      do 30 i = imin,imax+1
      ii = i - 3
      do 30 j = jmin,jmax+1
      jj = j - 3
      aa1(i,j) = dzx*pass(ii,jj)
   30 continue
!
!
!
!*********************************************************************
!
!.....define boundary values for itype=3
!
!*********************************************************************
!
!...top-left corner
!
      aa1(2,nz+2)=aa1(3,nzp)
!
!...top boundary
!
      j = nz+2
!cj dir$ ivdep
      do 721 i=3,nxp
      aa1(i,j) = aa1(i,j-1) + omdert(i)
  721 continue
!
!...top-right corner
!
      aa1(nx+2,nz+2) = aa1(nxp,nzp)
!
!...right boundary
!
      i = nx+2
!cj dir$ ivdep
      do 722 j=3,nzp
      aa1(i,j) = aa1(i-1,j) + omderr(j)
  722 continue
!
!...bottom-right corner
!
      aa1(nx+2,2)=aa1(nxp,3)
!
!...bottom boundary
!
      j=2
!cj dir$ ivdep
      do 723 i=3,nxp
      aa1(i,j)=aa1(i,j+1) - omderb(i)
  723 continue
!
!...bottom-left corner
!
      aa1(2,2)=aa1(3,3)
!
!...left boundary
!
      i=2
!cj dir$ ivdep
      do 724 j=3,nzp
      aa1(i,j)=aa1(i+1,j) - omderl(j)
  724 continue
   40 continue
!
!
!     DEALLOCATE(pass)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
