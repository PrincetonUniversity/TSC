      subroutine cdiv(ar,ai,br,bi,cr,ci)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 ar,ai,br,bi,cr,ci
 
!     complex division, (cr,ci) = (ar,ai)/(br,bi)
 
      REAL*8 s,ars,ais,brs,bis
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
