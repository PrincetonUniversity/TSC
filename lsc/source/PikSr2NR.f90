!
!     ------------------------------------------------------------------
!
      SUBROUTINE PikSr2NR ( n, arr, brr )
!     Sorts an array ARR of length N into acending numerical order,
!     by straight instertion --- while making the corresponding
!     rearrangement of the array BRR.
!     N is input; ARR is replaced on output
!     by its sorted rearrangement.  So is BRR.
!     Numerical Recipes in Fortran by W. H. Press, et al., p 228
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, j, i
      REAL*8    arr(n), brr(n), a, b
      do 12 j = 2, n
        a = arr(j)
        b = brr(j)
        do 11 i = j-1, 1, -1
            if ( arr(i) .le. a ) go to 10
            arr(i+1) = arr(i)
            brr(i+1) = brr(i)
 11     continue
        i = 0
 10     arr(i+1) = a
        brr(i+1) = b
 12   continue
      return
      END
!     -                                                                |
!     -                                                                |
!     Sorting by Straight Insertion   ends ----------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
