!#include "library_names.h"
      subroutine multipol
!
!.....Ref:  M.F.Reusch, G.H.Neilson, J. of Comput.Phys, 64 p416 (1986)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,icount,lln,l,ll,ig,lmax,i,ifail,j
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 gsumg,grrsumg,pmom,wk,deltij,grho,grho2,grho4,pminv
!     REAL*8 aident,summ,pmoms,grsumg,gzsumg,dsq,xminsq,xsq,xpos
      REAL*8 dsq,xminsq,xsq,xpos
      REAL*8 zpos,ratio1,ratio2,ratio3,ratio4,det,summin,summax
      REAL*8 sum
!============
!     dimension gsumg(pngroup,9),grrsumg(pngroup,9),
!    1    pmom(pngroup,6),wk(6),deltij(6,6),grho(pngroup,9),
!    3    grho2(pngroup,9), ngroupy(pngroup),grho4(pngroup,9),
!    4    pminv(pngroup,6),aident(6,6),summ(6),pmoms(pngroup,6),
!    5    grsumg(pngroup,9),gzsumg(pngroup,9)
      character*8 namul(6)
      data namul/'constant','dipole  ','quadpole',                       &  
     &    'hexapole','octapole','decapole'/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: gsumg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grrsumg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pmom
      REAL*8, ALLOCATABLE, DIMENSION(:) :: wk
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: deltij
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grho
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grho2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ngroupy
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grho4
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pminv
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: aident
      REAL*8, ALLOCATABLE, DIMENSION(:) :: summ
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pmoms
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grsumg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: gzsumg
!============      
      INTERFACE

      SUBROUTINE f03aae(a1,n1,n2,a2,a3,n3)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n1,n2,n3
      REAL*8 a2,d
      REAL*8 a3(*),a1(n1,*)
      INTEGER ipiv(n1)
      INTEGER i
      END SUBROUTINE f03aae

      SUBROUTINE f04aae(a1,n1,a2,n2,n3,n4,a3,n5,a4,n6)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n1,n2,n3,n4,n5,n6
      REAL*8 a1(n1,*),a2(n2,*),a3(n5,*),a4(*)
      INTEGER ipiv(n1)
      INTEGER i, j
      REAL*8 d
      END SUBROUTINE f04aae

      END INTERFACE

      IF(.not.ALLOCATED(gsumg)) ALLOCATE( gsumg(pngroup,9), STAT=istat)
      IF(.not.ALLOCATED(grrsumg)) ALLOCATE( grrsumg(pngroup,9),          &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(pmom)) ALLOCATE( pmom(pngroup,6), STAT=istat)
      IF(.not.ALLOCATED(wk)) ALLOCATE( wk(6), STAT=istat)
      IF(.not.ALLOCATED(deltij)) ALLOCATE( deltij(6,6), STAT=istat)
      IF(.not.ALLOCATED(grho)) ALLOCATE( grho(pngroup,9), STAT=istat)
      IF(.not.ALLOCATED(grho2)) ALLOCATE( grho2(pngroup,9), STAT=istat)
      IF(.not.ALLOCATED(ngroupy)) ALLOCATE( ngroupy(pngroup),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(grho4)) ALLOCATE( grho4(pngroup,9), STAT=istat)
      IF(.not.ALLOCATED(pminv)) ALLOCATE( pminv(pngroup,6), STAT=istat)
      IF(.not.ALLOCATED(aident)) ALLOCATE( aident(6,6), STAT=istat)
      IF(.not.ALLOCATED(summ)) ALLOCATE( summ(6), STAT=istat)
      IF(.not.ALLOCATED(pmoms)) ALLOCATE( pmoms(pngroup,6), STAT=istat)
      IF(.not.ALLOCATED(grsumg)) ALLOCATE( grsumg(pngroup,9),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(gzsumg)) ALLOCATE( gzsumg(pngroup,9),            &  
     &                                 STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : multipol  ' 
!============      
!
      dsq = .5_R8*xplas*deex
      xminsq = xplas**2 - 4._R8*dsq
      do 100 n=1,9
      xsq = xminsq + (n-1)*dsq
      xpos = sqrt(xsq)
      zpos = zplas
  100 call fieldg(xpos,zpos,gsumg(1,n),grsumg(1,n),gzsumg(1,n),          &  
     &     grrsumg(1,n),grho(1,n),                                       &  
     &           grho2(1,n))
      icount = 0
      do 200 lln=1,ngroupt
      l = nogroupt(lln)
      sum=0._R8
      do 210 n=1,9
  210 sum = sum + abs(gsumg(l,n))
      if(sum.eq.0) go to 200
!..rxw/23/04/87
!     check for linear dependence
      if(icount.eq.0) go to 211
      do 212 ll=1,icount
      ig=ngroupy(ll)
      if(l.eq.ig) go to 200
      ratio1=0._R8
      ratio2=0._R8
      ratio3=0._R8
      ratio4=0._R8
      if(gsumg(l ,6).ne.0._R8) ratio1=gsumg(l ,4)/gsumg(l ,6)
      if(gsumg(ig,6).ne.0._R8) ratio2=gsumg(ig,4)/gsumg(ig,6)
      if(gsumg(l ,6).ne.0._R8) ratio3=gsumg(l ,5)/gsumg(l ,6)
      if(gsumg(ig,6).ne.0._R8) ratio4=gsumg(ig,5)/gsumg(ig,6)
      if(ratio1.eq.ratio2.and. ratio3.eq.ratio4) go to 200
 212  continue
 211  continue
!..rxw/end
      icount = icount + 1
      ngroupy(icount) = l
!
      do 260 n=4,6
      grho4(l,n) = (grho2(l,n+1)-2._R8*grho2(l,n)+grho2(l,n-1))/dsq**2
  260 continue
      pmom(icount,1) = gsumg(l,5)/xplas**2
      pmom(icount,2) = grho(l,5)*2._R8
      pmom(icount,3) = grho2(l,5)*4._R8*xplas**2
      pmom(icount,4) = (grho2(l,6)-grho2(l,4))/(2._R8*dsq)*4._R8*xplas**  &  
     & 4
      pmom(icount,5) = (grho2(l,6)-2._R8*grho2(l,5)+grho2(l,4))/(dsq**2)  &  
!    &                                                                   &  
     &            *8._R8*xplas**6/3._R8
      pmom(icount,6) = (grho4(l,6)-grho4(l,4))/(2._R8*dsq)*4._R8*xplas**  &  
     & 8/3._R8
  200 continue
!
      write(nout,1001)
 1001 format(//,"group       null      dipole       quad        hex   ",  &  
     &                                                                   &  
     &  "      oct          dec")
      lmax = icount
      do 280 l=1,lmax
      write(nout,1002) ngroupy(l),(pmom(l,i),i=1,6)
 1002 format(i3,6e12.4)
      do 280 i=1,6
      pmoms(l,i) = pmom(l,i)
  280 continue
!
      ifail = 1
      n = 6
      if(lmax.lt.6) n=lmax
      do 385 i=1,n
      do 385 j=1,n
  385 deltij(i,j) = pmom(i,j)
      call f03aae(deltij,6,n,det,wk,ifail)
      if(det.eq.0 .or. ifail.ne.0) go to 411
!
      do 400 i=1,n
      do 390 j=1,n
      deltij(i,j) = 0._R8
  390 continue
      deltij(i,i) = 1._R8
  400 continue
      ifail=1
      call f04aae(pmom,pngroup,deltij,6,n,n,pminv,pngroup,wk,ifail)
      if(ifail.ne.0) ineg=15
!
      write(nout,2001)
 2001 format(//," inverse matrix   ... multipolar decomposition of"      &  
     & ," coil groups",/)
      do 290 l=1,n
  290 write(nout,2002) namul(l),(pminv(l,i),i=1,6),ngroupy(l)
 2002 format(a8,1p6e10.2,"  gr(",i2,")" )
 2003 format(a8,1p5e10.2,"  gr(",i2,")" )
!
!.....repeat with constant term omitted
      do 380 l=1,lmax
      do 380 i=1,5
      pmom(l,i) = pmoms(l,i+1)
  380 continue
      ifail = 1
      n = 5
      if(lmax.lt.5) n=lmax
      do 383 l=1,n
      summin = 0._R8
      summax = 0._R8
      sum = 0._R8
      do 382 i=1,n
  382 sum = sum + abs(pmom(l,i))
      summin = min(sum,summin)
      summax = max(sum,summax)
  383 continue
      if(n.le.1 .or. summax .gt. 1.E8_R8*summin) go to 411
!
      do 384 i=1,n
      do 384 j=1,n
  384 deltij(i,j) = pmom(i,j)
      call f03aae(deltij,6,n,det,wk,ifail)
      if(det.eq.0 .or. ifail.ne.0) go to 411
!
      do 500 i=1,n
      do 490 j=1,n
      deltij(i,j) = 0._R8
  490 continue
      deltij(i,i) = 1._R8
  500 continue
      ifail = 1
      call f04aae(pmom,pngroup,deltij,6,n,n,pminv,pngroup,wk,ifail)
      if(ifail.ne.0) ineg=15
      write(nout,2101)
 2101 format(//," multipolar decomposition without constant term")
      do 410 l=1,n
  410 write(nout,2003) namul(l+1),(pminv(l,i),i=1,5),ngroupy(l)
  411 continue
!
!
!
!.....print preprogrammed coil currents
      if(ntpts.le.0) return
      write(nout,3001)
 3001 format(//," preprogrammed field moments (KA equivalent) ",/,       &  
     &" time pts    time         nul         dip        quad         ",  &  
     &"hex         oct         dec")
      do 3030 l=1,ntpts
      do 3005 i=1,6
 3005 summ(i) = 0._R8
      do 3020 n=1,lmax
      ig = ngroupy(n)
      do 3010 i=1,6
 3010 summ(i) = summ(i) + pmoms(n,i)*gcur(l,ig)*.001_R8
 3020 continue
      write(nout,3031) l,tpros(l),(summ(i),i=1,6)
 3031 format(i5,1p7e12.4)
 3030 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
