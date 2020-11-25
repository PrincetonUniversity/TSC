 
!
!     -----------------------------------------------------------------|------
!
 
      SUBROUTINE WsloPrm ( uGiven , MuGiven , ZGiven,                    &  
     &                     dWsduou, iWhichWay)
 
!     Provides for a linear interpolation of Zeff in the region 1.0 <Zeff< 10.
!     The interpolation is calculated with nodes at Z = 1, 2, 5, 10 as given
!     in the Karney and Fisch paper.
 
!     Written by D. Enright, June 1992.
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iWhichWay
      REAL*8                                                             &  
     &        Z , uGiven, MuGiven, ZGiven,                               &  
     &        dWsduou, dWsduou1, dWsduou2, dWsduou5, dWsduou10
      REAL*8    ONE, TWO, FIVE, TEN
      DATA    ONE, TWO, FIVE, TEN/                                       &  
     &        1.0_R8, 2.0_R8,  5.0_R8,10.0_R8/
!
      Z    = ZGiven
!
!
      if (Z .le. ONE) then
         call WsloPrmZ(uGiven, MuGiven, ONE, dWsduou, iWhichWay)
      else if (Z .le. 2.0_R8) then
         call WsloPrmZ(uGiven, MuGiven, ONE, dWsduou1, iWhichWay)
         call WsloPrmZ(uGiven, MuGiven, TWO, dWsduou2, iWhichWay)
         dWsduou = (dWsduou2 - dWsduou1)/(TWO-ONE)*(Z-ONE) + dWsduou1
      else if (Z .le. 5.0_R8) then
         call WsloPrmZ(uGiven, MuGiven, TWO, dWsduou2, iWhichWay)
         call WsloPrmZ(uGiven, MuGiven,FIVE, dWsduou5, iWhichWay)
         dWsduou = (dWsduou5 - dWsduou2)/(FIVE-TWO)*(Z-TWO)+ dWsduou2
      else if (Z .le. TEN) then
         call WsloPrmZ(uGiven, MuGiven,FIVE, dWsduou5, iWhichWay)
         call WsloPrmZ(uGiven, MuGiven, TEN, dWsduou10, iWhichWay)
         dWsduou = (dWsduou10 - dWsduou5)/(TEN-FIVE)*(Z-FIVE)+ dWsduou5
      else
         call WsloPrmZ(uGiven, MuGiven, TEN , dWsduou, iWhichWay)
      endif
 
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
