      subroutine firevv
!
!...routine to generate a NO plasma zone on grid for use in
!...disruption simulations of FIRE design           3/7/2000
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j,jbot,jtop,ii,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rstart,zbot,ztop
!============
!     dimension ileft(200),iright(200)
      data rstart,zbot,ztop/1.750_R8,-1.450_R8,1.450_R8/
!============      
      INTEGER :: istat = 0 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ileft
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iright
!============      
      IF(.not.ALLOCATED(ileft)) ALLOCATE( ileft(200), STAT=istat)
      IF(.not.ALLOCATED(iright)) ALLOCATE( iright(200), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : firevv  ' 
!============      
!
      do 100 i=2,nxp-1
      if(xary(i) .le. rstart .and. xary(i+1) .gt. rstart) then
      istart=i
      go to 90
      endif
 100  continue
!
  90  continue
      do 101 j=2,nzp-1
      if(zary(j) .le. zbot .and. zary(j+1) .gt. zbot) then
      jbot=j
      go to 91
      endif
 101  continue
  91  continue
      do 102 j=2,nzp-1
      if(zary(j) .le. ztop .and. zary(j+1) .gt. ztop) then
      jtop=j
      go to 92
      endif
 102  continue
!
  92  continue
      do 103 j=jbot,jtop
      do 104 i=1,nxp
      ii=istart-(i-1)
      do 105 k=1,nwire
      if(iwire(k) .eq. ii .and. jwire(k) .eq. j) then
      ileft(j)=ii
      go to 106
      endif
 105  continue
 104  continue
 106  continue
!
      do 107 i=1,nxp
      ii=istart+i
      do 108 k=1,nwire
      if(iwire(k) .eq. ii .and. jwire(k) .eq. j) then
      iright(j)=ii
      go to 109
      endif
 108  continue
 107  continue
 109  continue
 103  continue
!
      do 110 j=2,nzp
      do 110 i=2,nxp
      iexv(i,j)=0
      if(j .lt. jbot .or. j .gt. jtop) then
      iexv(i,j)=1
      go to 110
      endif
      if(i .le. ileft(j)) iexv(i,j)=1
      if(i .ge. iright(j)) iexv(i,j)=1
 110  continue
!
      do 120 k=1,nwire
 120  iexv(iwire(k),jwire(k))=4
      call maskout
!
      do 121 k=1,nwire
 121  iexv(iwire(k),jwire(k))=1
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
