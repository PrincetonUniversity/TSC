        subroutine  tquit (iquit)
!..rxw/end
!-----------------------------------------------------------------------
!.....3.95.1 fluxmod
      USE CLINAM
      USE FVV1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iquit,iqtest,iabs
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 q95,vtest
      REAL*8 ptotmw,ohcur
!============
!       equivalence (q95, global(61))
!     dimension grsum(pngroup),grsum0(pngroup),gvsum(pngroup),
!    1          gvsum0(pngroup),gcur0ka(pngroup),gcurfka(pngroup)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grsum
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grsum0
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gvsum
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gvsum0
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gcur0ka
      REAL*8, ALLOCATABLE, DIMENSION(:) :: gcurfka
!============      
      IF(.not.ALLOCATED(grsum)) ALLOCATE( grsum(pngroup), STAT=istat)
      IF(.not.ALLOCATED(grsum0)) ALLOCATE( grsum0(pngroup), STAT=istat)
      IF(.not.ALLOCATED(gvsum)) ALLOCATE( gvsum(pngroup), STAT=istat)
      IF(.not.ALLOCATED(gvsum0)) ALLOCATE( gvsum0(pngroup), STAT=istat)
      IF(.not.ALLOCATED(gcur0ka)) ALLOCATE( gcur0ka(pngroup),            &  
     &                                 STAT=istat)
      IF(.not.ALLOCATED(gcurfka)) ALLOCATE( gcurfka(pngroup),            &  
     &                                 STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : tquit  ' 
!============      
!
!                       Test for plasma current and shape quantities
!                       outside range specified by acoef(102)
        q95 = global(61)
        iquit  = 0
        iqtest = int(acoef(101))
        if (iqtest.eq.0)   return
!
      if(igone.eq.1) return
!
        iquit = acoef(101)
        vtest  = acoef(102)
        if (iqtest.eq.-1 .and.    apl.lt.vtest)   return
        if (iqtest.eq. 1 .and.    apl.gt.vtest)   return
        if (iqtest.eq.-2 .and.   zmag.lt.vtest)   return
        if (iqtest.eq. 2 .and.   zmag.gt.vtest)   return
        if (iqtest.eq.-3 .and.  dipdt.lt.vtest)   return
        if (iqtest.eq. 3 .and.  dipdt.gt.vtest)   return
        if (iqtest.eq.-4 .and.   xmag.lt.vtest)   return
        if (iqtest.eq. 4 .and.   xmag.gt.vtest)   return
        if (iqtest.eq.-5 .and.    q95.lt.vtest)   return
        if (iqtest.eq. 5 .and.    q95.gt.vtest)   return
        if (iqtest.eq.-6 .and. shape5.lt.vtest)   return
        if (iqtest.eq. 6 .and. shape5.gt.vtest)   return
        if (iqtest.eq.-7 .and. shape3.lt.vtest)   return
        if (iqtest.eq. 7 .and. shape3.gt.vtest)   return
        ptotmw = (bpowers+pohmic+palpha)*1.E-6_R8
        if (iqtest.eq.-8 .and. ptotmw.lt.vtest)   return
        if (iqtest.eq. 8 .and. ptotmw.gt.vtest)   return
!
      iquit = 0
      if(iabs(iqtest) .lt. 9) return
      iquit = acoef(101)
      call groupcur(grsum,grsum0,gvsum,gvsum0,gcur0ka,gcurfka)
      ohcur = grsum(1)
      if (iqtest.eq.-9 .and. ohcur.lt.vtest)   return
      if (iqtest.eq. 9 .and. ohcur.gt.vtest)   return
!
        iquit = 0
        return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
