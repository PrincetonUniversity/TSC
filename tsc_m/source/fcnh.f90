      subroutine fcnh(x,y,f)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE CINTG

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,j
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 :: x
!     dimension y(201),f(201),veca(201),vecb(201),vecf(201)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: y
      REAL*8, ALLOCATABLE, DIMENSION(:) :: f
      REAL*8, ALLOCATABLE, DIMENSION(:) :: veca
      REAL*8, ALLOCATABLE, DIMENSION(:) :: vecb
      REAL*8, ALLOCATABLE, DIMENSION(:) :: vecf
!============      
      IF(.not.ALLOCATED(y)) ALLOCATE( y(201), STAT=istat)
      IF(.not.ALLOCATED(f)) ALLOCATE( f(201), STAT=istat)
      IF(.not.ALLOCATED(veca)) ALLOCATE( veca(201), STAT=istat)
      IF(.not.ALLOCATED(vecb)) ALLOCATE( vecb(201), STAT=istat)
      IF(.not.ALLOCATED(vecf)) ALLOCATE( vecf(201), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : fcnh  ' 
!============      
!     common/cintg/ahat(201,201),bhat(201,7),chat(7,201),fhat(201,7),
!    +gapdw(6),dgapdw(6)
!
      do 100 i=1,201
      veca(i)=0.0_R8
      vecb(i)=0.0_R8
      vecf(i)=0.0_R8
      do 101 j=1,201
      veca(i)=veca(i)+ahat(i,j)*y(j)
 101  continue
      do 102 j=1,6
      vecb(i)=vecb(i)+bhat(i,j)*gapdw(j)
      vecf(i)=vecf(i)+fhat(i,j)*dgapdw(j)
 102  continue
      f(i)=veca(i)+vecb(i)+vecf(i)
 100  continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
