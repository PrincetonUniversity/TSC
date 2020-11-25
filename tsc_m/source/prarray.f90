      subroutine prarray
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER jj,j,i
!============
      character*1 char1(2),char2(3)
      data char1/".","*"/
!     dimension jflag(pnx,pnz)
      data char2/"0","X","V"/
!============      
      INTEGER :: istat = 0 
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jflag
!============      
      IF(.not.ALLOCATED(jflag)) ALLOCATE( jflag(pnx,pnz), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : prarray  ' 
!============      
!
      write(nout,1010)
      do 200 jj=2,nzp
      j = 2+nzp-jj
      do 190 i=2,nxp
      jflag(i,j) = 0
      if(icoil(i,j) .ne. 0) jflag(i,j) = 1
  190 continue
      write(nout,1011) (char1(iexv(i,j)+1),i=2,nxp)
      if(j.eq.2) go to 200
      write(nout,1012) (char2(iexvc(i,j)+1),i=3,nxp)
  200 continue
!
 1010 format(1h1,"masking arrays")
 1011 format(1x,65(1a1,1x))
 1012 format(2x,64(1a1,1x))
      write(nout,1020)
 1020 format(1h1," coils")
!
      do 500 jj=2,nzp
      j = 2+nzp-jj
      write(nout,1011) (char1(jflag(i,j)+1),i=2,nxp)
      write(nout,1012) (char2(iexvc(i,j)+1),i=3,nxp)
  500 continue
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
