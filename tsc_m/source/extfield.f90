      subroutine extfield(val)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ixplas,i,n
      INTEGER imin1, imax1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 val,ans
      REAL*8 sum
!============
!     dimension psix(pnx)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psix
!============      
      IF(.not.ALLOCATED(psix)) ALLOCATE( psix(pnx), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : extfield  ' 
!============      
!
!...........................................................
!
!.....define external flux array on midplane
!
!...........................................................
      ixplas = (xplas - ccon)/deex + 2
      imin1 = ixplas - 1
      imax1 = ixplas + 1
!
      do 75 i=imin1,imax1
!
!.....upper half plane
!
      sum = 0
      do 78 n=1,nwire
      ans = 0._R8
      if(xary(i).eq.xary(iwire(n)).and.zary(jmag).eq.zary(jwire(n)))     &  
     &    go to 78
      call gf(ineg,nmult,xary(i),zary(jmag),xary(iwire(n))               &  
     &       ,zary(jwire(n)),ans)
   78 sum = sum + ans*ccoil(ncoil-nwire+n)/tpi
      if(ncoil.eq.nwire) go to 80
      do 79 n=1,ncoil-nwire
      call gf(ineg,nmult,xary(i),zary(jmag),xcoil(n),zcoil(n),ans)
   79 sum = sum + ans*ccoil(n)/tpi
   80 continue
      if(isym.eq.0) go to 180
!
!.....lower half plane
!
      do 178 n=1,nwire
      call gf(ineg,nmult,xary(i),zary(jmag),xary(iwire(n))               &  
     &       ,-zary(jwire(n)),ans)
      if(zary(jwire(n)).eq.0) ans = 0
  178 sum = sum + ans*ccoil(ncoil-nwire+n)/tpi
      if(ncoil.eq.nwire) go to 180
      do 179 n=1,ncoil-nwire
      call gf(ineg,nmult,xary(i),zary(jmag),xcoil(n),-zcoil(n),ans)
      if(zcoil(n).eq.0) ans = 0._R8
  179 sum = sum + ans*ccoil(n)/tpi
  180 continue
      psix(i) = sum
   75 continue
!...........................................................
!
      val = -(psix(ixplas+1)-psix(ixplas-1))/(2._R8*deex*xary(ixplas))
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
