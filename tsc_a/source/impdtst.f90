      subroutine impdtst
!
      USE CLINAM
      USE RADTAB
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nimp,indx,j,ii,nchrgs,jq,nchrgs1,i,is,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 anemks,animks,te
      REAL*8 radrmx,radrmn,delte,ainzrjq,radrjq,recrjq,pimlt,sumrt
      REAL*8 sum
!============
!     dimension radx(100),recx(100),ainzx(100),snp(100),aryte(100)
!    .,         radrt(100)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: radx
      REAL*8, ALLOCATABLE, DIMENSION(:) :: recx
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ainzx
      REAL*8, ALLOCATABLE, DIMENSION(:) :: snp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aryte
      REAL*8, ALLOCATABLE, DIMENSION(:) :: radrt
!============      
      IF(.not.ALLOCATED(radx)) ALLOCATE( radx(100), STAT=istat)
      IF(.not.ALLOCATED(recx)) ALLOCATE( recx(100), STAT=istat)
      IF(.not.ALLOCATED(ainzx)) ALLOCATE( ainzx(100), STAT=istat)
      IF(.not.ALLOCATED(snp)) ALLOCATE( snp(100), STAT=istat)
      IF(.not.ALLOCATED(aryte)) ALLOCATE( aryte(100), STAT=istat)
      IF(.not.ALLOCATED(radrt)) ALLOCATE( radrt(100), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : impdtst  ' 
!============      
!
!
!
! * * * plot comparison plots
!
!tss  call destroy('f4jdel0x')
      anemks = 1.E20_R8
      animks = 1.E18_R8
!tss  call fr80id('film-onl',1,2)
      call setpch(0,1,0,100)
!tss  call keep80('jdel')
!
      do 500 nimp = 1,pimp
!
      te  = 1.E-3_R8
      indx = 0
      radrmx = -1.E30_R8
      radrmn =  1.E30_R8
      do 300 j = 1,5
      delte = 10**(j-1)*1.E-3_R8
      do 310 ii = 1,9
      te = te + delte
      indx = indx + 1
      aryte(indx) = te
!
      nchrgs = nchrgsr(nimp)
      do 325 jq = 1,nchrgs
!
      call lookup (nimp,jq,te,anemks,ainzrjq,radrjq,recrjq)
!
      ainzx(jq) = anemks*ainzrjq
      radx(jq) = anemks*radrjq
      recx(jq) = anemks*recrjq
!
  325 continue
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
      pimlt = pimlt*recx(k+1)/ainzx(k)
      if(pimlt.gt.1.E40_R8) pimlt = 1.E40_R8
  240 continue
  230 sum = sum+pimlt
  220 continue
 
      snp(nchrgs) = animks/sum
      nn = nchrgs1
      do 250 i = 1,nchrgs1
      snp(nn) = recx(nn+1)*snp(nn+1)/ainzx(nn)
      nn = nn-1
  250 continue
!
      sumrt = 0.0_R8
      do 330 jq = 1,nchrgs
      radx(jq) = radx(jq)*snp(jq)
      sumrt = sumrt+radx(jq)
  330 continue
      radrt(indx) = 1.E13_R8*sumrt/(anemks*animks)
      radrmx = max(radrt(indx),radrmx)
      radrmn = min(radrt(indx),radrmn)
  310 continue
  300 continue
      call mapsll(0.001_R8, 100._R8, radrmn, radrmx, 0.2_R8, 0.8_R8,     &  
     & 0.3_R8, 0.7_R8)
      call tracec(1ht,aryte,radrt,indx,-1,-1,0._R8,0._R8)
      call setld(3._R8,20._R8,1,0,2,1)
      write(s100,1989)
      call gtext(s100,80,0)
 1989 format('cooling rate (ergs cm3/s)')
      call frame(0)
!
  500 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
