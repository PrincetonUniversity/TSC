      subroutine fieldg(xpos,zpos,gsumg,grsumg,gzsumg,                   &  
     &     grrsumg,grho,grho2)
!
!.....field from coil groups
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER m,n,ii,i2,j2,npass,iig,ig,iabs
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zpos,gsumg,grsumg,gzsumg,grrsumg,grho,grho2,xpos
      REAL*8 gsum,grsum,gzsum,grzsum
      REAL*8 gzzsum,grrsum
!============
!     dimension r1(2*pncoil),z1(2*pncoil),r2(2*pncoil),z2(2*pncoil),
!    1          gfun(2*pncoil),gr(2*pncoil),gz(2*pncoil),grz(2*pncoil),
!    2          gzz(2*pncoil),grr(2*pncoil)
!
!     dimension gsumg(pngroup),grrsumg(pngroup),grho(pngroup),           &  
!    &  grsumg(pngroup),gzsumg(pngroup)                                  &  
!    &          ,grho2(pngroup)
      dimension gsumg(*),grrsumg(*),grho(*),                             &  
     &  grsumg(*),gzsumg(*)                                              &  
     &          ,grho2(*)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: r1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: z1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: r2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: z2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gfun
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gz
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grz
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gzz
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grr
!============      
      IF(.not.ALLOCATED(r1)) ALLOCATE( r1(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(z1)) ALLOCATE( z1(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(r2)) ALLOCATE( r2(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(z2)) ALLOCATE( z2(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(gfun)) ALLOCATE( gfun(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(gr)) ALLOCATE( gr(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(gz)) ALLOCATE( gz(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(grz)) ALLOCATE( grz(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(gzz)) ALLOCATE( gzz(2*pncoil), STAT=istat)
      IF(.not.ALLOCATED(grr)) ALLOCATE( grr(2*pncoil), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : fieldg  ' 
!============      
      gsumg(1:pngroup)=0.0_R8
      grsumg(1:pngroup)=0.0_R8
      gzsumg(1:pngroup)=0.0_R8
      grrsumg(1:pngroup)=0.0_R8
      grho(1:pngroup)=0.0_R8
      grho2(1:pngroup)=0.0_R8
      
      m = 0
      do 10 n=1,ncoil
      m = m+1
      r1(m) = xpos
      z1(m) = zpos
      r2(m) = xcoil(n)
      z2(m) = zcoil(n)
      if(r1(m).eq.r2(m) .and. z1(m) .eq. z2(m)) return
      ii = n-ncoil+nwire
      if(ii.le.0) go to 9
      i2 = iwire(ii)
      j2 = jwire(ii)
      r2(m) = xary(i2)
      z2(m) = zary(j2)
      if(r1(m).eq.r2(m) .and. z1(m) .eq. z2(m)) return
    9 continue
      if(isym.eq.0) go to 10
      if(z2(m).eq.0) go to 10
      m = m+1
      r1(m) = xpos
      z1(m) = zpos
      r2(m) = r2(m-1)
      z2(m) =-z2(m-1)
!
   10 continue
      npass = m
      call gvect(r1,z1,r2,z2,npass,gfun,gr,gz,grz,gzz,grr,nmult,ineg)
      do 20 iig=1,ngroupt
      ig=nogroupt(iig)
      gsum = 0._R8
      grsum = 0._R8
      gzsum = 0._R8
      grzsum = 0._R8
      gzzsum = 0._R8
      grrsum = 0._R8
      ii = 0
      if(ncoil.eq.nwire) go to 31
      do 30 n=1,ncoil-nwire
      ii = ii + 1
      if(iabs(igroupc(n)) .ne. ig) go to 29
      gsum = gsum + gfun(ii)*aturnsc(n)/tpi
      grsum = grsum + gr(ii)*aturnsc(n)/tpi
      gzsum = gzsum + gz(ii)*aturnsc(n)/tpi
      grzsum = grzsum + grz(ii)*aturnsc(n)/tpi
      gzzsum = gzzsum + gzz(ii)*aturnsc(n)/tpi
      grrsum = grrsum + grr(ii)*aturnsc(n)/tpi
   29 if(isym.eq.0) go to 30
      if(z2(ii).eq.0) go to 30
      ii = ii + 1
      if(iabs(igroupc(n)) .ne. ig) go to 30
      gsum = gsum + gfun(ii)*aturnsc(n)/tpi
      grsum = grsum + gr(ii)*aturnsc(n)/tpi
      gzsum = gzsum + gz(ii)*aturnsc(n)/tpi
      grzsum = grzsum + grz(ii)*aturnsc(n)/tpi
      gzzsum = gzzsum + gzz(ii)*aturnsc(n)/tpi
      grrsum = grrsum + grr(ii)*aturnsc(n)/tpi
   30 continue
   31 continue
      if(nwire.le.0) go to 41
      do 40 n=1,nwire
      ii = ii + 1
      if(iabs(igroupw(n)) .ne. ig) go to 39
      gsum = gsum + gfun(ii)*aturnsw(n)/tpi
      grsum = grsum + gr(ii)*aturnsw(n)/tpi
      gzsum = gzsum + gz(ii)*aturnsw(n)/tpi
      grzsum = grzsum + grz(ii)*aturnsw(n)/tpi
      gzzsum = gzzsum + gzz(ii)*aturnsw(n)/tpi
      grrsum = grrsum + grr(ii)*aturnsw(n)/tpi
   39 if(isym.eq.0) go to 40
      if(z2(ii).eq.0) go to 40
      ii = ii + 1
      if(iabs(igroupw(n)) .ne. ig) go to 40
      gsum = gsum + gfun(ii)*aturnsw(n)/tpi
      grsum = grsum + gr(ii)*aturnsw(n)/tpi
      gzsum = gzsum + gz(ii)*aturnsw(n)/tpi
      grzsum = grzsum + grz(ii)*aturnsw(n)/tpi
      gzzsum = gzzsum + gzz(ii)*aturnsw(n)/tpi
      grrsum = grrsum + grr(ii)*aturnsw(n)/tpi
   40 continue
   41 continue
      if(xpos .eq. 0) go to 20
      gsumg(ig) = gsum
      grsumg(ig) = grsum
      gzsumg(ig) = gzsum
      grrsumg(ig) = grrsum
      grho(ig) = .5_R8*grsum/xpos
      grho2(ig) = -.25_R8*grsum/xpos**3 + .25_R8*grrsum/xpos**2
   20 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
