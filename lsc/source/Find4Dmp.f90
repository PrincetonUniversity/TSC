 
!     -----------------------------------------------------------------
 
      SUBROUTINE Find4Dmp(y, n, i4, nfound)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, n, i4(4), j, nfound
      REAL*8    lvl(4)
      REAL*8    y(n), dum, MxdPdz, start
      DATA    MxdPdZ, lvl(1), lvl(2), lvl(3), lvl(4)                     &  
     &      / 0.0001_R8,   0.75_R8,  0.50_R8,  0.25_R8, 0.10_R8/
 
      nfound = 1
      start = 1._R8
      j = 1
      do 5 i=1,4
         i4(i) = n
 5    continue
 
      do 10 i=2,n
          dum = 1._R8+ y(i)
!                                        Limit the damping per zone to the
!                                        value MxdPdZ
          if (dum .lt. MxdPdz) then
            dum = MxdPdZ
          endif
          start =  start * dum
          if (start .le. lvl(j) ) then
            i4(j) = i
            j=j+1
          endif
          if (j .eq. 5) go to 11
 10   continue
 
      nfound = j
      return
 
 11   continue
      nfound = 4
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
