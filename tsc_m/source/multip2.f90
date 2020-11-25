      subroutine multip2(anull,adipol,aquad,ahex,aoct,adec,              &  
     &          anull0,adipol0,aquad0,ahex0,aoct0,adec0)
!......5.74 multip2o
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,iii,ii,l
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 adipol,aquad,ahex,aoct,adec,anull0,adipol0,aquad0
      REAL*8 ahex0,aoct0,adec0,anull,gsumg,grrsumg,grho,grho2,grho4
      REAL*8 dsq,xminsq,xsq,xpos,zpos,curr
!============
      dimension gsumg(9),grrsumg(9),grho(9),grho2(9),grho4(9)
!     dimension gsumgg(pngroup,9) , grrsumgg(pngroup,9),
!    1          grhog(pngroup,9),grho2g(pngroup,9),grho4g(pngroup,9),
!    2  grsumgg(pngroup,9), gzsumgg(pngroup,9)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: gsumgg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grrsumgg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grhog
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grho2g
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grho4g
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: grsumgg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: gzsumgg
!============      
      IF(.not.ALLOCATED(gsumgg)) ALLOCATE( gsumgg(pngroup,9),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(grrsumgg)) ALLOCATE( grrsumgg(pngroup,9),        &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(grhog)) ALLOCATE( grhog(pngroup,9), STAT=istat)
      IF(.not.ALLOCATED(grho2g)) ALLOCATE( grho2g(pngroup,9),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(grho4g)) ALLOCATE( grho4g(pngroup,9),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(grsumgg)) ALLOCATE( grsumgg(pngroup,9),          &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(gzsumgg)) ALLOCATE( gzsumgg(pngroup,9),          &  
     &                                 STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : multip2  ' 
!============      
!
      dsq = .5_R8*xplas*deex
      xminsq = xplas**2 - 4._R8*dsq
      do 100 n=1,9
      xsq = xminsq + (n-1)*dsq
      xpos = sqrt(xsq)
      zpos = zplas
      call fieldg(xpos,zpos,gsumgg(1,n),grsumgg(1,n),gzsumgg(1,n),       &  
     &     grrsumgg(1,n),grhog(1,n),                                     &  
     &                      grho2g(1,n))
      call fieldc(xpos,zpos,gsumg(n),grrsumg(n),grho(n),grho2(n))
  100 continue
!
!
      do 260 n=4,6
      grho4(n) = (grho2(n+1)-2._R8*grho2(n)+grho2(n-1))/dsq**2
      do 260 iii=1,ngroupt
      ii = nogroupt(iii)
      grho4g(ii,n) = (grho2g(ii,n+1)-2._R8*grho2g(ii,n)+grho2g(ii,n-1))  &  
     &             /dsq**2
  260 continue
      anull = gsumg(5)/xplas**2
      adipol = grho(5)*2._R8
      aquad = grho2(5)*4._R8*xplas**2
      ahex = (grho2(6)-grho2(4))/(2._R8*dsq)*4._R8*xplas**4
      aoct = (grho2(6)-2._R8*grho2(5)+grho2(4))/(dsq**2)                 &  
     &            *8._R8*xplas**6/3._R8
      adec = (grho4(6)-grho4(4))/(2._R8*dsq)*4._R8*xplas**8/3._R8
!
!.....calculate multipolar moments from preprogrammed currents only
      anull0 = 0
      adipol0 = 0
      aquad0 = 0
      ahex0 = 0
      aoct0 = 0
      adec0 = 0
      do 200 iii=1,ngroupt
      ii = nogroupt(iii)
      do 200 l=1,ntpts
      curr = fact(l)*gcur(l,ii)*usdi
      anull0 = anull0 + curr*gsumgg(ii,5)/xplas**2
      adipol0 = adipol0 + curr*grhog(ii,5)*2._R8
      aquad0 = aquad0 +curr*grho2g(ii,5)*4._R8*xplas**2
      ahex0 = ahex0 +curr*(grho2g(ii,6)-grho2g(ii,4))*4._R8*xplas**4     &  
     &               /(2._R8*dsq)
      aoct0 = aoct0 +curr*(grho2g(ii,6)-2._R8*grho2g(ii,5)+grho2g(ii,4))  &  
!    &                                                                   &  
     &               /(dsq**2)*8._R8*xplas**6/3._R8
      adec0 = adec0 +curr*(grho4g(ii,6)-grho4g(ii,4))*4._R8*xplas**8/    &  
     & 3._R8                                                             &  
     &               /(2._R8*dsq)
  200 continue
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
