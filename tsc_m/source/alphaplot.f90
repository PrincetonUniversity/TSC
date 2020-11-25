      subroutine alphaplot
!
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ps,pval,ppval,alprmax,parymax
!============
!     dimension rave(ppsi), raveh(ppsi) ,yplot(ppsi), pary(ppsi)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rave
      REAL*8, ALLOCATABLE, DIMENSION(:) :: raveh
      REAL*8, ALLOCATABLE, DIMENSION(:) :: yplot
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pary
!============      
      IF(.not.ALLOCATED(rave)) ALLOCATE( rave(ppsi), STAT=istat)
      IF(.not.ALLOCATED(raveh)) ALLOCATE( raveh(ppsi), STAT=istat)
      IF(.not.ALLOCATED(yplot)) ALLOCATE( yplot(ppsi), STAT=istat)
      IF(.not.ALLOCATED(pary)) ALLOCATE( pary(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : alphaplot  ' 
!============      
!
      do 100 j=1,npsit
      rave(j) = sqrt((j-0.99999_R8)/(npsit-1))
      if(j.lt.2) go to 100
      raveh(j)= sqrt((j-1.5_R8)/(npsit-1))
      ps = xsv(j)
      call peval(ps,2,pval,ppval,imag,jmag)
      pary(j) = pval
  100 continue
      pary(1) = pary(2)
!
      call maps(0._R8,1.1_R8,0._R8,1.1_R8,.200_R8,.800_R8,.200_R8,       &  
     & .700_R8)
      call scalea(alphapr,yplot,npsit-1,1.0_R8,alprmax,1._R8)
      call colora("magenta")
      call tracec(1hA,raveh(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(pary(2),yplot(2),npsit-1,1.0_R8,parymax,1._R8)
      call colora("green")
      call tracec(1hT,raveh(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
!
      call setch(10._R8,31._R8,1,2,0,-1)
!
      call colora("magenta")
      write(s100,9001)alprmax
      call gtextm(s100,80,0,1,1)
!
      call colora("green")
      write(s100,9002) parymax
      call gtextm(s100,80,0,1,1)
!
      call setch(10._R8,4.0_R8,1,2,0,-1)
      call colora("red")
      write(s100,9007) kcycle,times
      call gtextm(s100,80,0,1,1)
 9007 format(" cycle=",i7,"  times=",1pe12.4)
      call colora("yellow")
!
      call frscj(14)
      return
 9001 format("A...Alpha pressure scaled to ",1pe12.4)
 9002 format("T...total pressure scaled to ",1pe12.4)
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
