!#include "f77_dcomplx.h"
      subroutine impinit
      USE CLINAM
      USE RADTAB
      USE SAPROP
      USE SCR1
      USE SPECIE
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nimp,nchrgs,j,jq,nchrgs1,i,is,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sumfrac,sumzfrc,tekev,ainzrjq,radrjq,recrjq,pimlt
      REAL*8 AREAL, sum
!============
!     dimension snp(pchrgmx)
!
      logical flgver
      data flgver / .false. /
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: snp
!============      
      IF(.not.ALLOCATED(snp)) ALLOCATE( snp(pchrgmx), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : impinit  ' 
!============      
!
      if(flgver) call impdtst
!
      sumfrac = 0._R8
      sumzfrc = 0._R8
      do 100 nimp = 1,pimp
      sumfrac = sumfrac + fraci(nimp)
      sumzfrc = sumzfrc + fraci(nimp)*AREAL(nchrgsr(nimp)-1)
  100 continue
      facimp = 1._R8/(sumzfrc + 1._R8-sumfrac)
!
      do 500 nimp = 1,pimp
      if(fraci(nimp).le.0.0_R8) go to 500
      nchrgs = nchrgsr(nimp)
!
!  * * impurity ionization
!
      do 200 j = 2,npsit
!     write(6,7777) j
 7777 format( "j=",i3,"in impinit")
!
      do 210 jq = 1,nchrgs
      tekev = te(j)*1.E-3_R8
      call lookup (nimp,jq,tekev,ane(j),ainzrjq,radrjq,                  &  
     &recrjq)
      ainz(jq,j) = ane(j)*ainzrjq
      if(ainz(jq,j) .lt. 1.E-8_R8) ainz(jq,j) = 1.E-8_R8
      rad(nimp,jq,j) = ane(j)*radrjq
      rec(jq,j) = ane(j)*recrjq
  210 continue
 
!
! * * LTE model
!
      sum = 0.0_R8
      nchrgs1 = nchrgs-1
      do 220 i = 1,nchrgs
      pimlt = 1.0_R8
      is = i
      if(i.eq.nchrgs) go to 230
      do 240 k = is,nchrgs1
!     write(6,8111) k,j,pimlt,rec(k+1,j),ainz(k,j)
      pimlt = pimlt*rec(k+1,j)/ainz(k,j)
!     if(pimlt.gt.1.e40) pimlt = 1.e40
 8111 format(2i3,1p3e12.4)
  240 continue
!     if(sum.gt.1.e60) go to 220
  230 sum = sum+pimlt
  220 continue
!     write(6,8777) sum
 8777 format(" sum = ",1pe12.4)
 
      snp(nchrgs) = facimp*fraci(nimp)*ane(j)*vp(j)/sum
      nn = nchrgs1
      do 250 i = 1,nchrgs1
!     write(6,2222) i,nimp,j,rec(nn+1,j),snp(nn+1),ainz(nn,j)
 2222 format(3i3,1p3e16.8)
      snp(nn) = rec(nn+1,j)*snp(nn+1)/ainz(nn,j)
      nn = nn-1
  250 continue
      do 410 i = 1,nchrgs
      nq(i,nimp,j) = snp(i)
  410 continue
!
  200 continue
!     stop
!
      do 600 j = npsit+1,npsi
      nq(1,nimp,j) = facimp*fraci(nimp)*ane(j)*vp(j)
      do 610 i = 2,nchrgs
      nq(i,nimp,j) = 0.0_R8
  610 continue
  600 continue
  500 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
