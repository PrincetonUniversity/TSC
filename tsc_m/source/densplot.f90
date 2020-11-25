      subroutine densplot
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
      REAL*8 dperpmax,anemax,scalef,sravebmax
      REAL*8 sravemax,simpemax,sredgmax,srjetmax
!============
!     dimension rave(ppsi), raveh(ppsi) ,yplot(ppsi)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rave
      REAL*8, ALLOCATABLE, DIMENSION(:) :: raveh
      REAL*8, ALLOCATABLE, DIMENSION(:) :: yplot
!============      
      IF(.not.ALLOCATED(rave)) ALLOCATE( rave(ppsi), STAT=istat)
      IF(.not.ALLOCATED(raveh)) ALLOCATE( raveh(ppsi), STAT=istat)
      IF(.not.ALLOCATED(yplot)) ALLOCATE( yplot(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : densplot  ' 
!============      
!
      do 100 j=1,npsit
      rave(j) = sqrt((j-0.99999_R8)/(npsit-1))
      if(j.le.2) go to 100
      raveh(j)= sqrt((j-1.5_R8)/(npsit-1))
  100 continue
!
      call maps(0._R8,1.1_R8,0._R8,1.1_R8,.200_R8,.800_R8,.200_R8,       &  
     & .700_R8)
      call scalea(dperpa,yplot,npsit-1,1.0_R8,dperpmax,1._R8)
      call colora("magenta")
      call tracec(1hD,rave,yplot,npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(ane(2),yplot(2),npsit-1,1.0_R8,anemax,1._R8)
      call colora("green")
      call tracec(1hn,raveh(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
!.....scalefactor to convert from dimensionless units to MKS
      scalef = udsd/udst
      call scalea(sraveb(2),yplot(2),npsit-1,0.5_R8,sravebmax,scalef)
      call colora("red")
      call tracec(1hB,raveh(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(srave (2),yplot(2),npsit-1,0.5_R8,sravemax,scalef )
      call colora("blue")
      call tracec(1hP,raveh(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(simpe (2),yplot(2),npsit-1,0.5_R8,simpemax,scalef )
      call colora("cyan")
      call tracec(1hI,raveh(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(sraveedg (2),yplot(2),npsit-1,0.5_R8,sredgmax,scalef )   
      call colora("yellow")
      call tracec(1hE,raveh(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call scalea(sravejet (2),yplot(2),npsit-1,0.5_R8,srjetmax,scalef )    
      call colora("white")
      call tracec(1hJ,raveh(2),yplot(2),npsit-1,-1,-1,0._R8,0._R8)
!
      call setch(10._R8,31._R8,1,2,0,-1)
!
      call colora("magenta")
      write(s100,9001)dperpmax
      call gtextm(s100,80,0,1,1)
!
      call colora("green")
      write(s100,9002) anemax
      call gtextm(s100,80,0,1,1)
!
      call colora("red")
      write(s100,9003) sravebmax
      call gtextm(s100,80,0,1,1)
!
      call colora("blue")
      write(s100,9004) sravemax
      call gtextm(s100,80,0,1,1)
!
      call colora("cyan")
      write(s100,9005) simpemax
      call gtextm(s100,80,0,1,1)
!
      call colora("yellow")
      write(s100,9006) sredgmax
      call gtextm(s100,80,0,1,1)
!
      call colora("white")
      write(s100,9008) srjetmax
      call gtextm(s100,80,0,1,1)
!
      call setch(10._R8,4.0_R8,1,2,0,-1)
      call colora("red")
      write(s100,9007) kcycle,times
      call gtextm(s100,80,0,1,1)
 9007 format(" cycle=",i7,"  times=",1pe12.4)
      call colora("yellow")
!
      call frscj(13)
      return
 9001 format("D...diffusion coeff scaled to ",1pe12.4)
 9002 format("n...density profile scaled to ",1pe12.4)
 9003 format("B...beam source scaled to ",1pe12.4)
 9004 format("P...pellet source scaled to ",1pe12.4)
 9005 format("I...impurity source scaled to ",1pe12.4)
 9006 format("E...edge source scaled to ",1pe12.4)
 9008 format("J...jet source scaled to ",1pe12.4)
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
