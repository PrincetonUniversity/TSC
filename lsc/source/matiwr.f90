!
!     -----------------------------------------------------------------
!
      SUBROUTINE matiwr(title, comat, n1dim, n1, i1l, n2, i2l)
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER *(*) title
      CHARACTER *(30) label
      INTEGER n1dim, n1, i1l, i1u, n2, i2l, i2u, cmaxr, cminr, i1, i2,   &  
     &     comatr
      INTEGER comat(n1dim, n2)
      write(nLSCcom2,41)title
 41   format(/,/,a, '(i1, i2)',/)
      write(nLSCcom2,4)n1dim,n1,n2, i1l, i2l
 4    format(' n1dim = ', i3, 3x,' n1 = ', i3, 3x, ' n2 = ',i3 /         &  
     &       ' i1l = ', i3, 3x, ' i2l = ', i3)
!
!     find extrema and print
!
      cmaxr = comat(1, 1)
      cminr = comat(1, 1)
        i1u = i1l + n1 - 1
        i2u = i2l + n2 - 1
      do 5 i1 = i1l, i1u
      do 5 i2 = i2l, i2u
      comatr=comat(i1,i2)
      if(cmaxr.lt.comatr)cmaxr=comatr
      if(cminr.gt.comatr)cminr=comatr
 5    continue
      write(nLSCcom2,50)cminr,cmaxr
50    format(' min = ', i3, 5x,' max = ', i3)
        do 1 i2 = i2l, i2u
           write(label, 101) i2
 101       format('i2 = ', i2)
           call viwrite(n1, comat(1, i2), label)
1     continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
