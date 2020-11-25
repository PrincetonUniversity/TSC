!#include "f77_dcomplx.h"
      subroutine fixshape(iter)
!
!.....06/17/96    scj&np
!     called during the equilibrium iteration to define feedback
!     group currents gcurfb to fix the plasma shape.
!.....called for acoef(901) > 0.
!
!.....shape points are input on type 62 and 63 time point cards
!62 xcon   i    xcon(i,1)   xcon(i,2)   xcon(i,3)  ...  xcon(i,6)
!63 zcon   i    zcon(i,1)   zcon(i,2)   zcon(i,3)  ...  zcon(i,6)
!
!.for acoef(901)=1   only shape points are used
!.for acoef(901)=2   shape points + flux linkage (acoef(902)) at xplas,zplas
!.for acoef(901)=3   shape points, + x-point (r_x=acoef(903),z_x=acoef(904))
!.for acoef(901)=4   shape points + flux linkage + x-point
!
!.acoef(905) specifies max number of coil group currents to calculate
!.(actually set this equal to the total number of groups. To fix any
!.coil current, set the desired value in gcur(2) and set the corresponding
!.gcur(3) value to 1.0)
!.acoef(906) is the relative error tolorance...[1.e-3]
!.acoef(907) is the iteration number when shape feedback starts
!.acoef(908) is the iteration number when type 19 feedback ends
!.acoef(909) is the relaxation factor for equilibrium shape feedback
!.acoef(910) is number of iterations between resetting sigmax and relaxation
!            factors
!.acoef(911) is number of iterations to full implementation of vsec constraint
!
!..must use istart=1 for this option
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE CLINAM
      USE CBSVD

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iter,itestdiag,jtime,i,j,null,l,n,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tau,wmax,wmin,wrat,rsq,rmsnorm,gradsq,dpsidx
      REAL*8 dpsidz,gsval,psval,psixz,psixx,psizz,sigmin,add
      REAL*8 sum, AREAL
!============
!     common/cbsvd/sigmax,nranksc,ngscfb
!     dimension amatsc(pngroup,pngroup),bvecsc(pngroup)
!     dimension gsumg(pngroup),grsumg(pngroup),gzsumg(pngroup),
!    1 grrsumg(pngroup),grho(pngroup),grho2(pngroup),
!    2 xvecsc(pngroup),xvecsco(pngroup),
!    3 bvecsco(pngroup)
!     dimension uuscfb(pngroup,pngroup),vvscfb(pngroup,pngroup),
!    1          rv1scfb(pngroup)
!     dimension groupwt(pngroup)
      data itestdiag/0/
      data jtime/0/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: amatsc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: bvecsc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gzsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grrsumg
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grho2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xvecsc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xvecsco
      REAL*8, ALLOCATABLE, DIMENSION(:) :: bvecsco
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: uuscfb
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: vvscfb
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rv1scfb
      REAL*8, ALLOCATABLE, DIMENSION(:) :: groupwt
!============      
      IF(.not.ALLOCATED(amatsc)) ALLOCATE( amatsc(pngroup,pngroup),      &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(bvecsc)) ALLOCATE( bvecsc(pngroup), STAT=istat)
      IF(.not.ALLOCATED(gsumg)) ALLOCATE( gsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grsumg)) ALLOCATE( grsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(gzsumg)) ALLOCATE( gzsumg(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grrsumg)) ALLOCATE( grrsumg(pngroup),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(grho)) ALLOCATE( grho(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grho2)) ALLOCATE( grho2(pngroup), STAT=istat)
      IF(.not.ALLOCATED(xvecsc)) ALLOCATE( xvecsc(pngroup), STAT=istat)
      IF(.not.ALLOCATED(xvecsco)) ALLOCATE( xvecsco(pngroup),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(bvecsco)) ALLOCATE( bvecsco(pngroup),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(uuscfb)) ALLOCATE( uuscfb(pngroup,pngroup),      &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(vvscfb)) ALLOCATE( vvscfb(pngroup,pngroup),      &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(rv1scfb)) ALLOCATE( rv1scfb(pngroup),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(groupwt)) ALLOCATE( groupwt(pngroup),            &  
     &                                 STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : fixshape  ' 
!============      
      jtime = jtime + 1
!
      do 7001 i=1,pngroup
      do 7000 j=1,pngroup
      uuscfb(i,j) = 0._R8
      vvscfb(i,j) = 0._R8
      amatsc(i,j) = 0._R8
 7000 continue
      bvecsc(i) = 0._R8
      sigmasvd(i) = 0._R8
      xvecsc(i) = 0._R8
      rv1scfb(i) = 0._R8
      groupwt(i)=gcur(2,i)
      if(i.le.acoef(905).and.groupwt(i).eq.0) then
      print *,'***warning*** gcur(2,l) should be non zero'
      stop
      endif
 7001 continue
      tau = 0._R8
      wmax = 0._R8
      wmin = 0._R8
      wrat = 0._R8
      rsq = 0._R8
      rmsnorm = 0._R8
      rmsnorm = 0._R8
      null = 0
!      ngscfb = ifix(acoef(905))
      ngscfb=ngroupt
      if(acoef(906).le.0) then
      write(nterm,*) 'acoef(906) .le. 0!!!'
      endif
      if(ncnt .gt. 0 .and. ngscfb .gt. 0 .and. ngscfb .le. ngroupt)      &  
     &     go to 100
      write(nout,8888) ncnt,ngscfb
 8888 format(" error in equilibrium shape control feedback",             &  
     &" ncnt,ngscfb =",2i5)
      write(nterm,8888) ncnt,ngscfb
      ineg=2
      return
!
  100 continue
!
! initialize xvecso
!
      if(jtime.eq.1) then
! if the x-point constraint is applied, add the x-point position to the
! set of position constraints
      if(acoef(901).eq.3._R8.or. acoef(901).eq.4._R8) then
      xcon0(istart,ncnt+1)=acoef(903)
      zcon0(istart,ncnt+1)=acoef(904)
      ncnt=ncnt+1
      endif
      do 897 l=1,ngscfb
      xvecsco(l) = groupwt(l)*gcur(istart,l)*usdi
  897 continue
      endif
!
!.....shape constraint
      nranksc = ncnt
      do 200 n=1,ncnt
      call fieldg(xcon0(istart,n),zcon0(istart,n),gsumg,grsumg,          &  
     &            gzsumg,grrsumg,grho,grho2)
      call grap(1,zcon0(istart,n),xcon0(istart,n),gradsq,dpsidx,dpsidz,  &  
     &     gsval,psval,psixz,psixx,psizz,1)
!
      do 199 l=1,ngscfb
! subtract off coil contrib to flux
      psval = psval - gsumg(l)*(xvecsco(l)/groupwt(l)+gcurfbo(l))
  199 amatsc(n,l) = gsumg(l)/groupwt(l)
      bvecsc(n) = +(psilimo - psval)
!      bvecsc(n)=+(psilim - (psval - bvecsco(n)))
      bvecsco(n)=bvecsc(n)
      if(itestdiag.eq.1 .or. lrswtch.gt.0) bvecsc(n) = -psval
  200 continue
!      write(nterm,4439) (bvecsc(l),l=1,ncnt)
! 4439 format(" bvecsc",1p10e12.4)
      if(acoef(901) .ne. 2._R8.and. acoef(901) .ne. 4._R8) go to 300
!
!.....flux linkage constraint
      nranksc = nranksc + 1
      call fieldg(xplas,zplas,gsumg,grsumg,                              &  
     &            gzsumg,grrsumg,grho,grho2)
      do 299 l=1,ngscfb
  299 amatsc(nranksc,l) = gsumg(l)*tpi/groupwt(l)
      if(iter.lt.(acoef(911)+1)) then
      bvecsc(nranksc) = acoef(902)*(iter-1)/acoef(911)
      endif
      if(iter.ge.(acoef(911)+1)) bvecsc(nranksc) = acoef(902)
      bvecsco(nranksc) = bvecsc(nranksc)
  300 continue
      if(acoef(901) .ne. 3._R8.and. acoef(901) .ne. 4._R8) go to 400
!
!.....x-point constraint
!      print *,'x-point constraint'
      call fieldg(acoef(903),acoef(904),gsumg,grsumg,                    &  
     &            gzsumg,grrsumg,grho,grho2)
      call grap(2,acoef(904),acoef(903),gradsq,dpsidx,dpsidz,            &  
     &     gsval,psval,psixz,psixx,psizz,1)
!
      nranksc = nranksc + 1
      do 398 l=1,ngscfb
! subtract off coil contribution to flux
      dpsidx = dpsidx - grsumg(l)*(xvecsco(l)/groupwt(l)+gcurfbo(l))
  398 amatsc(nranksc,l) = grsumg(l)/groupwt(l)
      bvecsc(nranksc) = - dpsidx
      bvecsco(nranksc) = bvecsc(nranksc)
!
      nranksc = nranksc + 1
      do 399 l=1,ngscfb
! subtract off coil contribution to flux
      dpsidz = dpsidz - gzsumg(l)*(xvecsco(l)/groupwt(l)+gcurfbo(l))
  399 amatsc(nranksc,l) = gzsumg(l)/groupwt(l)
      bvecsc(nranksc) = - dpsidz
      bvecsco(nranksc) = bvecsc(nranksc)
!
  400 continue
!
!...fixed coil current constraint
! if gcur(3,i)=1, that coil current is to be constrained equal
! to its value set in gcur(1,i).
      do i=1,pngroup
      if(gcur(3,i).eq.1.E3_R8) then
      nranksc=nranksc+1
      do 491 l=1,ngscfb
      amatsc(nranksc,l)=0._R8
      if(l.eq.i) amatsc(nranksc,l)=1._R8/groupwt(l)
 491  continue
      bvecsc(nranksc)=gcur(1,i)*usdi
      endif
      enddo
!
      if(itestdiag.eq.1) go to 799
!
!......invert matrix and solve for the gcurfb
!
! copy a into u
      do k=1,nranksc
         do l=1,ngscfb
            uuscfb(k,l)=amatsc(k,l)
         enddo
      enddo
      call svdcmp(uuscfb,nranksc,ngscfb,pngroup,pngroup,sigmasvd,        &  
     & vvscfb)
!     find maximum singular value
      sigmax=0.0_R8
      do 13 k=1,ngscfb
         if (sigmasvd(k).gt.sigmax) sigmax=sigmasvd(k)
 13   continue
      if(jtime.eq.1) then
      write(nterm,6679)
 6679 format('singular values for A matrix are:')
      write(nterm,6680) (sigmasvd(k),k=1,ngscfb)
 6680 format(1p7e12.4)
      write(nterm,6681)
 6681 format('values scaled to sigmax=1. are:')
      write(nterm,6680) (sigmasvd(k)/sigmax,k=1,ngscfb)
      write(nterm,6682)
 6682 format(' ')
      write(nterm,6683)
 6683 format('enter sigma cutoff ')
      read(nterm,*) acoef(906)
      endif
!     define "small"
      sigmin=sigmax*acoef(906)
!     zero the "small" singular values
      do 14 k=1,ngscfb
         if (sigmasvd(k).lt.sigmin) sigmasvd(k)=0.0_R8
 14   continue
      call svbksb(uuscfb,sigmasvd,vvscfb,nranksc,ngscfb,                 &  
     &   pngroup,pngroup,bvecsc,xvecsc)
! test how accurately the svd was done
      if(iter/10*10.eq.iter) then
      sum=0._R8
      do 17 k=1,nranksc
         bvecsc(k)=0.0_R8
         do 16 j=1,ngscfb
            bvecsc(k)=bvecsc(k)+amatsc(k,j)*xvecsc(j)
 16      continue
      add=(bvecsc(k)-bvecsco(k))**2
      if(bvecsco(k).ne.0._R8) add=add/bvecsco(k)**2
      sum=sum+add
 17   continue
      resid=sqrt(sum)/AREAL(nranksc)
      write(nterm,6684) resid
 6684 format('svd residual = ',1pe12.4)
      endif
!
      fluxlink=0._R8
      if(jtime.le.1) go to 501
      call fieldg(xplas,zplas,gsumg,grsumg,                              &  
     &            gzsumg,grrsumg,grho,grho2)
      do 499 l=1,ngscfb
! subtract off preprogrammed part and apply relaxation factor
      gcurfb(l) = gcurfb(l) +                                            &  
     &            xvecsco(l)/groupwt(l) - gcur(istart,l)*usdi            &  
     &          + acoef(909)*(xvecsc(l)-xvecsco(l))/groupwt(l)
      xvecsco(l) = xvecsco(l) + acoef(909)*(xvecsc(l)-xvecsco(l))
      fluxlink=fluxlink+gsumg(l)*tpi*(gcurfb(l)+gcur(istart,l)*usdi)
  499 continue
      if(iter/10*10.eq.iter) then
      write(nterm,5535) (gcurfb(l),l=1,ngscfb)
 5535 format(" gcurfb",1p10e12.4)
      endif
  501 psilimo = psilim
      return
!
!.....coding for diagnostic test for itestdiag=1
!
  799 continue
      write(nterm,1800)
 1800 format(" diagnostic test of equilibrium shape control")
      do 801 n=1,nranksc
      sum = 0._R8
      do 800 l=1,ngscfb
  800 sum = sum + amatsc(n,l)*gcur(istart,l)*usdi
      write(nterm,1801) n,sum,bvecsc(n)
 801  continue
 1801 format(i5,1p2e12.4)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
