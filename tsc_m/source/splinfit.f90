      subroutine splinfit
!
!
!.....special subroutine for loading splin data for ifunc=6
!
      USE CLINAM
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ifield,nrecsf,ii,ll,l,m,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ane1,anh1,ani1,zk,zc
      REAL*8 ank1,anc1,pvalmin,evalmin,dis,sine,cose,cr1,cy1,am2
      REAL*8 cr2,cr4,cy4,am1,denom
!============
!     dimension yy(ppsi,4,pneq),xx(ppsi,4)
!     dimension pinpt(ppsi),einpt(ppsi),rinpt(ppsi)
!     dimension yinpt(ppsi,5)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: yy
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xx
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pinpt
      REAL*8, ALLOCATABLE, DIMENSION(:) :: einpt
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rinpt
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: yinpt
!============      
      IF(.not.ALLOCATED(yy)) ALLOCATE( yy(ppsi,4,pneq), STAT=istat)
      IF(.not.ALLOCATED(xx)) ALLOCATE( xx(ppsi,4), STAT=istat)
      IF(.not.ALLOCATED(pinpt)) ALLOCATE( pinpt(ppsi), STAT=istat)
      IF(.not.ALLOCATED(einpt)) ALLOCATE( einpt(ppsi), STAT=istat)
      IF(.not.ALLOCATED(rinpt)) ALLOCATE( rinpt(ppsi), STAT=istat)
      IF(.not.ALLOCATED(yinpt)) ALLOCATE( yinpt(ppsi,5), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : splinfit  ' 
!============      
      ifield = 4000
!
      nrecsf = acoef(ifield)
      npsit6 = nrecsf
      do 89 ii=1,nrecsf
      xs6(ii+1) = acoef(ii+ifield)
   89 continue
      ifield = ifield + nrecsf
      do  95 ll=1,5
      do 90 ii=1,nrecsf
      yinpt(ii,ll) = acoef(ii+ifield)
   90 continue
      ifield = ifield + nrecsf
   95 continue
!
!....redefine initial impurity densities
      ane1 = yinpt(1,1)
      anh1 = yinpt(1,2)
      ani1 = yinpt(1,3)
!
!....Krypton and Carbon charge
      zk = 36._R8
      zc = 6._R8
!....Krypton relative density
      ank1 = (ane1 - anh1 - zc*ani1)/((anh1+ani1)*(zk-zc))
      anc1 = (ane1 - anh1 - zk*ani1)/((anh1+ani1)*(zc-zk))
!
      write(nterm,7111) fraci(2),anc1
 7111 format(" Carbon  fraction changed from",1pe12.4," to",1pe12.4)
      fraci(2) = anc1
      write(nterm,7112) fraci(6),ank1
 7112 format(" Krypton fraction changed from",1pe12.4," to",1pe12.4)
      fraci(6) = ank1
!
      rinpt(nrecsf+2) = fracn0
      pvalmin = fracn0*acoef(882)*r0*smallt
      pinpt(nrecsf+2) = pvalmin
      evalmin = fracn0*r0*smallt
      einpt(nrecsf+2) = evalmin
      xs6(1) = 0._R8
      xs6(nrecsf+2) = 2._R8*xs6(nrecsf+1) - xs6(nrecsf)
!
      do 98 l=2,nrecsf+1
      rinpt(l) = yinpt(l-1,1)*1.E6_R8*usdd
      einpt(l) = yinpt(l-1,1)*yinpt(l-1,4)*1.E6_R8*usdh
      pinpt(l) =(yinpt(l-1,1)*yinpt(l-1,4) +                             &  
     &          (yinpt(l-1,2)+yinpt(l-1,3))*yinpt(l-1,5))*1.E6_R8*usdh
   98 continue
      rinpt(1) = rinpt(2)
      einpt(1) = einpt(2)
      pinpt(1) = pinpt(2)
!     write(nterm,4444)
!4444 format(" rinpt, einpt, pinpt follow")
!c    write(nterm,4445) (pinpt(i),i=1,nrecsf+2)
!     write(nterm,4445) (einpt(i),i=1,nrecsf+2)
!     write(nterm,4445) (rinpt(i),i=1,nrecsf+2)
!     write(nterm,4445) (xs6(i),i=1,nrecsf+2)
!4445 format(1p10e12.4)
!
!
!
      do 100 m=1,4
      do 100 l=2,npsit6
      yy(l,m,2) = pinpt(l-2+m)
      yy(l,m,3) = einpt(l-2+m)
      yy(l,m,4) = rinpt(l-2+m)
      xx(l,m) = xs6(l-2+m)
  100 continue
!
      do 70 k=2,pneq
      do 200 l=3,npsit6
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
      as6(l,k) = (xx(l,2)**2*(yy(l,3,k)-am1*xx(l,3))                     &  
     & + xx(l,3)**2*(yy(l,2,k)-am2*xx(l,2))                              &  
     &         - xx(l,2)*xx(l,3)*(yy(l,2,k)+yy(l,3,k)) )/denom
      bs6(l,k) = ((yy(l,2,k) - yy(l,3,k))*(xx(l,2)-xx(l,3))              &  
     & + xx(l,2)**2*am1 + xx(l,3)**2*am2                                 &  
     &         + 2._R8*xx(l,2)*xx(l,3)*(am1+am2)  )/denom
      cs6(l,k) =-(xx(l,2)*(2._R8*am1+am2) + xx(l,3)*(am1+2._R8*am2))/    &  
     & denom
      ds6(l,k) = (am1+am2)/denom
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
      as6(l,k) = (xx(l,2)**2*(yy(l,3,k)-am1*xx(l,3))                     &  
     & + xx(l,3)**2*(yy(l,2,k)+am1*xx(l,2))                              &  
     &         - xx(l,2)*xx(l,3)*(yy(l,2,k)+yy(l,3,k)) )/denom
      bs6(l,k) = ((yy(l,2,k) - yy(l,3,k))*(xx(l,2)-xx(l,3))              &  
     & + am1*(xx(l,2)**2 - xx(l,3)**2) )/denom
      cs6(l,k) =-am1*(xx(l,2) - xx(l,3))/denom
      ds6(l,k) = 0._R8
   70 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
