!#include "f77_dcomplx.h"
      subroutine residcal(residd,imax,jmax,residmax)
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER imax,jmax,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 residmax,residd,amag,dpdz
      REAL*8 dpdx,ajmid,gmid,dfdz,xmidsq,dfdx,denom,termd
      REAL*8 AREAL
!============
!     dimension residz(pnx,pnz), residx(pnx,pnz),extra1(pnx,pnz),
!    1          extra2(pnx,pnz)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: residz
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: residx
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: extra1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: extra2
!============      
      IF(.not.ALLOCATED(residz)) ALLOCATE( residz(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(residx)) ALLOCATE( residx(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(extra1)) ALLOCATE( extra1(pnx,pnz), STAT=istat)
      IF(.not.ALLOCATED(extra2)) ALLOCATE( extra2(pnx,pnz), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : residcal  ' 
!============      
!
      residd = 0._R8
      amag = 0._R8
      residmax = 0._R8
      do 333 i=iminn,imaxx
      do 333 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 333
      dpdz = .25_R8*(psi(i,j+1)+psi(i-1,j+1)-psi(i,j-1)-psi(i-1,j-1))
      dpdx = .25_R8*(psi(i+1,j)+psi(i+1,j-1)-psi(i-1,j)-psi(i-1,j-1))
      ajmid = .5_R8*(ajphi(i,j)+ajphi(i-1,j))*xarh(i)
      gmid = .5_R8*(g(i,j+1)+g(i,j))*xsqoj(i)
      dfdz = (ajmid*dpdz                                                 &  
     &         + (gmid)*(g(i,j+1)-g(i,j))*xsqoj(i)                       &  
     &         + xarh(i)**2*(pr(i,j+1)-pr(i,j)))/deez
      ajmid = .5_R8*(ajphi(i,j)+ajphi(i,j-1))*xary(i)
      gmid = .5_R8*(g(i+1,j)*xsqoj(i+1)+g(i,j)*xsqoj(i))
      xmidsq = xary(i)**2
      dfdx = (ajmid*dpdx                                                 &  
     &         + (gmid)*(g(i+1,j)*xsqoj(i+1)-g(i,j)*xsqoj(i))            &  
     &         + xmidsq*(pr(i+1,j)-pr(i,j)))/deex
      residz(i,j) = dfdz
      residx(i,j) = dfdx
      denom = (dpdz**2 + dpdx**2) + 1.E-12_R8
      termd =  (dpdz*dfdz+dpdx*dfdx)**2/denom
      residd = residd + (dpdz*dfdz+dpdx*dfdx)**2/denom
      amag = amag + (ajmid*(dpdz/deez + dpdx/deex))**2
      if(termd.gt.residmax) then
      residmax = termd
      imax = i
      jmax = j
      endif
  333 continue
      if(amag.gt.0) residd = sqrt(residd/amag)
      if(igone.eq.0) return
      do i=2,nxp
      do j=2,nzp
      extra1(i,j) = gs(i,j) - gzero
      extra2(i,j) = AREAL(iexv(i,j))
      enddo
      enddo
      call diag(residx,' residx   ')
      call diag(residz,' residz   ')
      call diag(b     ,' vorticity')
      call diag(ajphi ,' current  ')
      call diag(extra1,'small g   ')
      call diag(extra2,'iexv      ')
      call diag(abig  ,'abig      ')
      do i=2,nxp
      do j=2,nzp
      extra2(i,j) = g(i,j) - gzero/xsqoj(i)
      extra1(i,j) = AREAL(iexvc(i,j))
      enddo
      enddo
!
      call diag(u     ,' div-V    ')
      call diag(w     ,' tor V    ')
      call diag(extra2,' tor flux ')
      call diag(omeg  ,' omeg     ')
      call diag(extra1,' iexvc    ')
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
