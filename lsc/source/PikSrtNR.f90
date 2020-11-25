!
!     Sorting by Straight Insertion   begins --------------------------|
!     -                                                                |
!     -                                                                |
      SUBROUTINE PikSrtNR ( n, arr )
!     Sorts an array ARR of length N into acending numerical order,
!     by straight instertion.  N is input; ARR is replaced on output
!     by its sorted rearrangement.
!     Numerical Recipes in Fortran by W. H. Press, et al., p 226
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, j, i
      REAL*8    arr(n), a
      do 12 j = 2, n
        a = arr(j)
        do 11 i = j-1, 1, -1
            if ( arr(i) .le. a ) go to 10
            arr(i+1) = arr(i)
 11     continue
        i = 0
 10     arr(i+1) = a
 12   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
