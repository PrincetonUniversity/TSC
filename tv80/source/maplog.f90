!#include "f77_dcomplx.h"
        subroutine maplog (xlow, xhi)
 
!*************************************************************************
!
!  maplog  -  Convert log coordinates of a map
!
!  synopsis     call maplog (xlow,xhi)
!
!               float xlow                      Low end of range
!               float xhi                       Hi end of range
!
!  description  Converts the low and high end of a range into a log
!               range that will have integer exponents.  This algorithm
!               was taken directly from tv80lib.
!
!*************************************************************************
 
!
!c purpose: make sure that xlow is less than xhi
!           and generate the log of the values if necessary
!
! if xlow eq xhi, make xhi a smidge bigger
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL log10a
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   xhi,xlow,t,g1
      REAL   REAL
!============
        if (xlow .gt. xhi) then
          t = xlow
          xlow = xhi
          xhi = t
        else if (xlow .eq. xhi) then
          xhi = xlow + 1.0 
        endif
!
! reset xlow and xhi to non-negative values
!
        if (xlow .lt. 0.0 ) xlow = 0.0 
        if (xhi .lt. 0.0 ) xhi = 0.0 
!
! reset xlow and xhi to default values if xhi has an bad value
!
        if (xhi .eq. 0) then
          xlow = .00001 
          xhi = 100000. 
          call tgerror (8)
          go to 110
        endif
 
! reset xlow with respect to xhi
 
        if (xlow .eq. 0.0 ) then
          xlow = xhi*0.0000000001 
          call tgerror (9)
        endif
 
! get log value for xlow
 
110     xlow = log10a(xlow)
 
! get log value for xhi
 
        g1 = alog10(xhi)
        if (g1 .gt. 0.0 ) g1 = g1+0.9999999999 
        xhi = REAL(int(g1))
        xlow = REAL(int(xlow))
 
! reset xlow
 
        if (xlow .eq. xhi) xlow = xhi-1.0 
!
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
