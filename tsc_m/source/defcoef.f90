!#include "library_names.h"
      subroutine defcoef
!......number 8,70
!
!.....this subroutine uses the adiabatic variables adi,ade,adn,adp
!.....and the geometrical variables xmja2,vp2 to define local coefficien
!.....of a cubic interpolating polynomial.  these coefficients asv,bsv,c
!.....are then used in the subroutines geval,peval,eeval,reval to define
!.....toroidal field function,pressure,electron pressure, and density
!.....as functions of the poloidal flux
!
!
      USE CLINAM
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER jj,j,m,l,k,ier,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 atemp,btemp,ctemp,wkspace,dis,sine,cose,cr1
      REAL*8 cy1,am2,cr2,cr4,cy4,am1,denom
!============
!     dimension yy(ppsi,4,pneq),xx(ppsi,4),gsq(ppsi)
      dimension atemp(4,4),btemp(4),ctemp(4),wkspace(4)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: yy
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xx
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gsq
!============      
      INTERFACE
      subroutine f04aae(a1,n1,a2,n2,n3,n4,a3,n5,a4,n6)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n1,n2,n3,n4,n5,n6
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 a1(n1,*),a2(n2,*),a3(n5,*),a4(*)
      INTEGER ipiv(n1)
      INTEGER i, j
      REAL*8 d
      END subroutine f04aae

      END INTERFACE
!============
      IF(.not.ALLOCATED(yy)) ALLOCATE( yy(ppsi,4,pneq), STAT=istat)
      IF(.not.ALLOCATED(xx)) ALLOCATE( xx(ppsi,4), STAT=istat)
      IF(.not.ALLOCATED(gsq)) ALLOCATE( gsq(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : defcoef  ' 
!============      
!
!......integrate differential equation to determine g
      gsq(npsit+1) = gzero**2
      gsq(npsit+2) = gzero**2
      do 60 jj=2,npsit
      j = npsit+2-jj
      if(j.le.npsit+1) go to 55
      gsq(j) = gzero**2
      go to 60
   55 continue
      gsq(j) = gsq(j+1) + 2._R8/xmja2(j) * ( (gxmja(j+1)-gxmja(j))       &  
     &   + tpi*qprof2(j)*vp2(j)*(adp(j+1)/vpg(j+1)-adp(j)/vpg(j) ) )
      if(gsq(j) .gt. 0) go to 60
      gsq(j) = gsq(j+1)
!     ineg=15
      write(nterm,3111) jj,j,qprof2(j),vp2(j),gsq(j),gsq(j+1),           &  
     &                  xmja2(j),gxmja(j+1),gxmja(j),vp2(j),             &  
     &                  adp(j+1),vpg(j+1),adp(j),vpg(j)
 3111 format("jj,j=",2i3,"   qprof2,vp2,gsq,gsq",1p4e12.4,               &  
     &                 /,1p8e12.4)
   60 continue
      gsq(1) = gsq(2)
!
      do 100 m=1,4
      do 100 l=2,npsit
      yy(l,m,1) = sqrt(gsq(l-2+m))
      yy(l,m,2) = adp(l-2+m)/vpg(l-2+m)
      yy(l,m,3) = ade(l-2+m)/vpg(l-2+m)
      yy(l,m,4) = adn(l-2+m)/vp(l-2+m)
      xx(l,m) = xsv(l-2+m)
  100 continue
      do 70 k=1,pneq
      do 200 l=3,npsit
      dis = sqrt((yy(l,2,k)-yy(l,3,k))**2 + (xx(l,2)-xx(l,3))**2)
      sine = (yy(l,2,k) - yy(l,3,k))/dis
      cose = (xx(l,2)    - xx(l,3)   )/dis
!
      cr1 = cose*(xx(l,1)-xx(l,3)) + sine*(yy(l,1,k)-yy(l,3,k))
      cy1 =-sine*(xx(l,1)-xx(l,3)) + cose*(yy(l,1,k)-yy(l,3,k))
      am2 = cy1/cr1
!
      cr2 = cose*(xx(l,2)-xx(l,3)) + sine*(yy(l,2,k)-yy(l,3,k))
      cr4 = cose*(xx(l,4)-xx(l,3)) + sine*(yy(l,4,k)-yy(l,3,k))
      cy4 =-sine*(xx(l,4)-xx(l,3)) + cose*(yy(l,4,k)-yy(l,3,k))
      am1 = -cy4/(cr2-cr4)
      denom = (xx(l,2)-xx(l,3))**2
      asv(l,k) = (xx(l,2)**2*(yy(l,3,k)-am1*xx(l,3))                     &  
     & + xx(l,3)**2*(yy(l,2,k)-am2*xx(l,2))                              &  
     &         - xx(l,2)*xx(l,3)*(yy(l,2,k)+yy(l,3,k)) )/denom
      bsv(l,k) = ((yy(l,2,k) - yy(l,3,k))*(xx(l,2)-xx(l,3))              &  
     & + xx(l,2)**2*am1 + xx(l,3)**2*am2                                 &  
     &         + 2._R8*xx(l,2)*xx(l,3)*(am1+am2)  )/denom
      csv(l,k) =-(xx(l,2)*(2._R8*am1+am2) + xx(l,3)*(am1+2._R8*am2))/    &  
     & denom
      dsv(l,k) = (am1+am2)/denom
  200 continue
!
!.....force fit to be quadratic near origin
      l=2
      dis = sqrt((yy(l,2,k)-yy(l,3,k))**2 + (xx(l,2)-xx(l,3))**2)
      sine = (yy(l,2,k) - yy(l,3,k))/dis
      cose = (xx(l,2)   - xx(l,3)  )/dis
!
      cr2 = cose*(xx(l,2)-xx(l,3)) + sine*(yy(l,2,k)-yy(l,3,k))
      cr4 = cose*(xx(l,4)-xx(l,3)) + sine*(yy(l,4,k)-yy(l,3,k))
      cy4 =-sine*(xx(l,4)-xx(l,3)) + cose*(yy(l,4,k)-yy(l,3,k))
      am1 = -cy4/(cr2-cr4)
      denom = (xx(l,2)-xx(l,3))**2
      asv(l,k) = (xx(l,2)**2*(yy(l,3,k)-am1*xx(l,3))                     &  
     & + xx(l,3)**2*(yy(l,2,k)+am1*xx(l,2))                              &  
     &         - xx(l,2)*xx(l,3)*(yy(l,2,k)+yy(l,3,k)) )/denom
      bsv(l,k) = ((yy(l,2,k) - yy(l,3,k))*(xx(l,2)-xx(l,3))              &  
     & + am1*(xx(l,2)**2 - xx(l,3)**2) )/denom
      csv(l,k) =-am1*(xx(l,2) - xx(l,3))/denom
      dsv(l,k) = 0._R8
!...added 3/21/96 during killer pellet work
      atemp(1,1) = 0._R8
      atemp(1,2) = 1._R8
!     atemp(1,3) = 2._R8*xsv2(1)
!     atemp(1,4) = 3._R8*xsv2(1)**2
      atemp(1,3) = 2.*xx(l,1)
      atemp(1,4) = 3.*xx(l,1)**2
      atemp(2,1) = 0._R8
      atemp(2,2) = 1._R8
      atemp(2,3) = 2._R8*xx(l,3)
      atemp(2,4) = 3._R8*xx(l,3)**2
      atemp(3,1) = 1._R8
      atemp(3,2) = xx(l,2)
      atemp(3,3) = xx(l,2)**2
      atemp(3,4) = xx(l,2)**3
      atemp(4,1) = 1._R8
      atemp(4,2) = xx(l,3)
      atemp(4,3) = xx(l,3)**2
      atemp(4,4) = xx(l,3)**3
      btemp(1) = 0._R8
      btemp(2) = (yy(l,4,k) - yy(l,2,k))/(xx(l,4) - xx(l,2))
      btemp(3) = yy(l,2,k)
      btemp(4) = yy(l,3,k)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      ier    = 0
      call f04aae(atemp,4,btemp,4,4,1,ctemp,4,wkspace,ier)
      if(ier.ne.0) then
      do j=1,4
        write(nout,1001) j,(atemp(j,i),i=1,4),btemp(j),ctemp(j)
      enddo
 1001 format("error in defcoef",i4,1p6e12.4)
      endif
      if(ier.ne.0) ineg=15
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      asv(l,k) = ctemp(1)
      bsv(l,k) = ctemp(2)
      csv(l,k) = ctemp(3)
      dsv(l,k) = ctemp(4)
   70 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
