      subroutine q8plot
 
!**************************************************************************
!
!  q8plot
!
!  calls        trace
!               tracep
!
!**************************************************************************
 
!
!c author: f. n. fritsch
!
!c variable declarations:
!
      USE q7quad
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,j1,isplat,knumch,lpts
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   t,xi,yi,tt,ave
!============
!     common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
!    * , int, k01, k02, c01, cdf, czmin, czmax
!
! coding for quad begins here...
      dimension ba(5)
      dimension bb(5)
      dimension t(5)
      REAL*8 ba, bb
      save
!
!c program statements:
!
! load t
      do 1 i = 1, 4
    1 t(i) = czc(i)-ctest
      t(5) = t(1)
      j1 = 1
      isplat = 0
!
! begin loop
      do 9 knumch = 1, 4
    2 if (t(knumch)) 3, 4, 5
!
    3 if (t(knumch+1)) 9, 9, 6
!
    4 xi = cxc(knumch)
      yi = cyc(knumch)
      isplat = isplat+1
      go to 7
!
    5 if (t(knumch+1)) 6, 9, 9
!
! interpolate
    6 tt =   t(knumch)/(t(knumch)-t(knumch+1))
      xi = (cxc(knumch+1)-cxc(knumch))*tt+cxc(knumch)
      yi = (cyc(knumch+1)-cyc(knumch))*tt+cyc(knumch)
    7 ba(j1) = xi
      bb(j1) = yi
      j1 = j1+1
    9 continue
!
! switch for closure and correct number of lines..
   12 go to (9999, 9999, 14, 15, 17), j1
!
   14 lpts = 2
      go to 30
!
   15 lpts = j1
      ba(j1) = ba(1)
      bb(j1) = bb(1)
      go to 30
!
! decide which two lines to plot..
! jump if this quadrilateral is a plateau
   17 if (isplat .eq. 4) go to 9999
      ave = 0.25 *(t(1)+t(2)+t(3)+t(4))
      if (ave*t(1)) 20, 15, 22
!
   20 ba(4) = ba(2)
      bb(4) = bb(2)
      ba(2) = xi
      bb(2) = yi
   22 lpts = 2
!
! end of subroutine quad
!
! now plot the contours that have been found...
      if (int .eq. 0) go to 23
      call trace (ba(3), bb(3), 2, -1, -1, 0.0 , 0.0 )
      go to 31
!
   23 call tracep (ba(3), bb(3), 2, 3, -1, -1)
      go to 32
!
   30 if (int .eq. 0) go to 32
   31 call trace (ba, bb, lpts, -1, -1, 0.0 , 0.0 )
      go to 9999
!
   32 call tracep (ba, bb, lpts, 3, -1, -1)
!
9999  return
 
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
