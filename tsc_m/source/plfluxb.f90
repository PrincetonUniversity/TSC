      subroutine plfluxb(xb,zb,ans)
!......4.61 plfluxb
!
!***********************************************************************
!                                                                      *
!...calculate plasma contribution to poloidal flux at boundary         *
!...point (xb,zb) using full volume integral method.                   *
!                                                                      *
!***********************************************************************
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER m,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zb,ans,xb,curf
      REAL*8 sum
!============
!     dimension xtpass(2*pnx*pnz),ztpass(2*pnx*pnz),xspass(2*pnx*pnz),
!    1          zspass(2*pnx*pnz),anspass(2*pnx*pnz)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xtpass
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ztpass
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xspass
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zspass
      REAL*8, ALLOCATABLE, DIMENSION(:) :: anspass
!============      
      IF(.not.ALLOCATED(xtpass)) ALLOCATE( xtpass(2*pnx*pnz),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(ztpass)) ALLOCATE( ztpass(2*pnx*pnz),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(xspass)) ALLOCATE( xspass(2*pnx*pnz),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(zspass)) ALLOCATE( zspass(2*pnx*pnz),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(anspass)) ALLOCATE( anspass(2*pnx*pnz),          &  
     &                                 STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : plfluxb  ' 
!============      
!
!
      if(lrswtch.gt.0) return
      m = 0
      do 99 i=iminn,imaxx
      do 99 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 99
      if(psi(i,j).ge.psilim) go to 99
      m = m+1
      xtpass(m) = xb
      ztpass(m) = zb
      xspass(m) = xary(i)
      zspass(m) = zary(j)
      if(isym.eq.0) go to 99
      if(zspass(m).eq.0) go to 99
      m = m+1
      xtpass(m) = xb
      ztpass(m) = zb
      xspass(m) = xary(i)
      zspass(m) =-zary(j)
   99 continue
!
      call grnfncv(ineg,nmult,xtpass,ztpass,xspass,zspass,m,anspass)
!
      m = 0
      sum = 0
      do 98 i=iminn,imaxx
      do 98 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 98
      if(psi(i,j).ge.psilim) go to 98
      m = m+1
      curf = ajphi(i,j)*dxdz/(2.0_R8*pi)
      sum = sum + curf*anspass(m)
      if(isym.eq.0) go to 98
      if(zary(j).eq.0) go to 98
      m = m + 1
      sum = sum + curf*anspass(m)
   98 continue
      ans = sum
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
