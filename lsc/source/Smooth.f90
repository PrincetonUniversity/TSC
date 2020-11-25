!                                                                      |
!     power.F ---------------------------------------------------------+
 
 
!     ql.f(or)  begins                      ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE Smooth(rmat, n1, n1dim, n2, nsmoo, smvec)
      USE params
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ipsi, nsm2, n1, n2, nsmoo, n1dim
      REAL*8 rmat(n1dim, n2), smvec(n1dim)
!     smooth quasilinear diffusion coefficient
!     VecDmp in matr.F puts the Dql values for each v into wkv.
!
!     The calling logic is tricky in that it depends on storage order in f77.
!     convolve takes wkv, outputs Dql.  If qlsm is a delta fn, no change.
!     qlsm is integer vector gaussian of width nsmw, normalized to add to 1.
 
      nsm2 = (nsmoo - 1) / 2
      do 10 ipsi=1, n2
         call VecDmp(rmat(1,ipsi), wkv, n1, 1)
         call convolve(n1, nsm2, smvec, rmat(1, ipsi), wkv)
 10   continue
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
