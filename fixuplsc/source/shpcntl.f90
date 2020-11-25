!     subroutine shpcntl(ig,nl,alphac,term)
!     include 'clinam.i'
!     ineg = 50
!     write(*,*) "subroutine shpcntl called but not supplied"
!     return
!     end
      subroutine minv(a,n,ndim,wrk1,det,epsminv,imode,jmode)
      USE PARAM
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,ndim,imode,jmode,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 wrk1,det,epsminv,a,epsmin,d
!============
      dimension a(ndim,n),wrk1(ndim)
!    &         ,wk(pngroup,pngroup),
!    1          a1(pngroup,pngroup),a2(pngroup,pngroup)
!     integer   indx(pngroup)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: wk, a1, a2
      INTEGER, ALLOCATABLE, DIMENSION(:)  :: indx
      INTEGER :: istat =0

      if(.not.ALLOCATED(wk)) ALLOCATE(wk(pngroup,pngroup),STAT=istat)
      if(.not.ALLOCATED(a1)) ALLOCATE(a1(pngroup,pngroup),STAT=istat)
      if(.not.ALLOCATED(a2)) ALLOCATE(a2(pngroup,pngroup),STAT=istat)
      if(.not.ALLOCATED(indx)) ALLOCATE(indx(pngroup),                   &  
     &                                STAT=istat)
      if( istat .ne. 0 ) stop 'Allocation Error : minv'
      if(n.eq.0) return
      do 20 i=1,n
      do 10 j=1,n
      a1(i,j) = a(i,j)
   10 continue
   20 continue
      epsmin = 1.E-13_R8
!     call f01abf(a1,pngroup,n,a2,pngroup,wrk1,ifail)
!     do 100 i=1,n
!     do 90  j=1,i
!     a(i,j) = a1(i+1,j)
!  90 continue
! 100 continue
!kuma st
      call ludcmp(a1,n,pngroup,indx,d)
      do 102 i=1,n
        do 101 j=1,n
          a2(i,j) = 0.0_R8
101     continue
        a2(i,i) = 1.0_R8
102   continue
      do 103 j=1,n
          call lubksb(a1,n,pngroup,indx,a2(1,j))
103   continue
      do 100 i=1,n
      do 90  j=1,i
      a(i,j) = a2(i,j)
   90 continue
  100 continue
!kuma en
      do 110 i=1,n-1
      do 105 j=i+1,n
      a(i,j) = a(j,i)
  105 continue
  110 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
