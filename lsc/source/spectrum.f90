 
!     spectrum                              ---------------------------|
!                                                                      |
!                                                                      |
       REAL*8 FUNCTION spectrum (ParIdx)
      USE params
      USE RayWrk
       IMPLICIT NONE
       INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
       REAL*8                                                            &  
     &        ParIdx, a, c
       INTEGER istart
!     ParIdx    Parallel refractive index, n_\parallel
!     a         width parameter
!     c         normalization constant
!     istart    0 on loading, 1 after first call
!     spectrum  function based on f(x) = c * 1/(x^2 + a^2)
!               where c= a/PI if the range were +/- infinity
       REAL*8 spwidth, spcentr
       DATA      istart/0/
       DATA spwidth /0.0_R8/, spcentr/0.0_R8/
 
       if (istart .eq. 0 ) then
         istart = 1
         a = spwidth/2.
         c = a/( atan((nparmax-spcentr)/a)                               &  
     &        - atan((nparmin-spcentr)/a) )
       endif
 
       spectrum = c/( (ParIdx-spcentr)**2 + a*a )
       return
 
       END
!c                                                                      |
!c                                                                      |
!c     spectrum ends                         ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
