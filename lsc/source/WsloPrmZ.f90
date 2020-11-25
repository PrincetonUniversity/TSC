 
!
!     -----------------------------------------------------------------|-------
!
      SUBROUTINE WsloPrmZ( uGiven , MuGiven , ZGiven,                    &  
     &                     dWsduou, iWhichWay)
!     Ratio of power coupled form the rf source into electromagnetic energy to
!     the rf power absorbed by electrons; P_el / P_in.
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iWhichWay
      REAL*8                                                             &  
     &        u , Mu , Z , uGiven, MuGiven, ZGiven,                      &  
     &        x, dWsduou
      REAL*8    e01, e02, e03, e04, e1, e2, A0,  A1,  A2
      REAL*8                                                             &  
     &    WPpZ01a1, WPpZ01a2, WPpZ01a3, WPpZ01b1, WPpZ01b2, WPpZ01b3 ,   &  
     &    WPpZ02a1, WPpZ02a2, WPpZ02a3, WPpZ02b1, WPpZ02b2, WPpZ02b3 ,   &  
     &    WPpZ05a1, WPpZ05a2, WPpZ05a3, WPpZ05b1, WPpZ05b2, WPpZ05b3 ,   &  
     &    WPpZ10a1, WPpZ10a2, WPpZ10a3, WPpZ10b1, WPpZ10b2, WPpZ10b3
      DATA                                                               &  
     &    WPpZ01a1, WPpZ01a2, WPpZ01a3, WPpZ01b1, WPpZ01b2, WPpZ01b3 /   &  
     &     0.66445_R8, -0.36032_R8,  0.07328_R8,  0.17769_R8,            &  
     &     -0.25452_R8,  0.07278_R8/
      DATA                                                               &  
     &    WPpZ02a1, WPpZ02a2, WPpZ02a3, WPpZ02b1, WPpZ02b2, WPpZ02b3 /   &  
     &     0.56760_R8, -0.38984_R8,  0.08634_R8, -0.04019_R8,            &  
     &     -0.24673_R8,  0.08508_R8/
      DATA                                                               &  
     &    WPpZ05a1, WPpZ05a2, WPpZ05a3, WPpZ05b1, WPpZ05b2, WPpZ05b3 /   &  
     &     0.39906_R8, -0.32879_R8,  0.07670_R8, -0.28281_R8,            &  
     &     -0.16275_R8,  0.07436_R8/
      DATA                                                               &  
     &    WPpZ10a1, WPpZ10a2, WPpZ10a3, WPpZ10b1, WPpZ10b2, WPpZ10b3 /   &  
     &     0.27028_R8, -0.23261_R8,  0.05272_R8, -0.39140_R8,            &  
     &     -0.07526_R8,  0.04981_R8/
 
      REAL*8                                                             &  
     &    WPmZ01a1, WPmZ01a2, WPmZ01a3, WPmZ01a4 ,                       &  
     &    WPmZ02a1, WPmZ02a2, WPmZ02a3, WPmZ02a4 ,                       &  
     &    WPmZ05a1, WPmZ05a2, WPmZ05a3, WPmZ05a4 ,                       &  
     &    WPmZ10a1, WPmZ10a2, WPmZ10a3, WPmZ10a4
 
      DATA                                                               &  
     &    WPmZ01a1, WPmZ01a2, WPmZ01a3, WPmZ01a4 /                       &  
     &    -0.63673_R8, -1.39960_R8,  3.37662_R8, -4.23684_R8/
      DATA                                                               &  
     &    WPmZ02a1, WPmZ02a2, WPmZ02a3, WPmZ02a4 /                       &  
     &    -0.55777_R8, -0.80763_R8,  1.43144_R8, -2.03866_R8/
      DATA                                                               &  
     &    WPmZ05a1, WPmZ05a2, WPmZ05a3, WPmZ05a4 /                       &  
     &    -0.39704_R8, -0.33811_R8,  0.23607_R8, -0.51011_R8/
      DATA                                                               &  
     &    WPmZ10a1, WPmZ10a2, WPmZ10a3, WPmZ10a4 /                       &  
     &    -0.26600_R8, -0.17342_R8,  0.01896_R8, -0.13349_R8/
!
!
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
      iWhichWay = 0
!
!     .....
!
!     For Mu = 1 and 0 < u < 5 use Table IV. on page 191.
      if ( Mu .ge. 0.90_R8) then
 
             if ( u .gt. 5.00_R8) then
!     .                                 If the velocity is too large, limit it
!     .                                 to a value covered by the table, and
!     .                                 report WhichWay the velocity is large.
               u  = 5.00_R8
               iWhichWay = +1
             endif
!
             x = u*u
!
             if ( Z  .le. 1.5_R8)  then
                dWsduou =                                                &  
     &           (   WPpZ01a1*x    + WPpZ01a2*x**2 + WPpZ01a3*x**3 ) /   &  
     &       ( 1._R8+  WPpZ01b1*x    + WPpZ01b2*x**2 + WPpZ01b3*x**3 )
 
        else if ( Z  .le. 3.0_R8)  then
                dWsduou =                                                &  
     &           (   WPpZ02a1*x    + WPpZ02a2*x**2 + WPpZ02a3*x**3 ) /   &  
     &       ( 1._R8+  WPpZ02b1*x    + WPpZ02b2*x**2 + WPpZ02b3*x**3 )
 
        else if ( Z  .le. 7.0_R8)  then
                dWsduou =                                                &  
     &           (   WPpZ05a1*x    + WPpZ05a2*x**2 + WPpZ05a3*x**3 ) /   &  
     &       ( 1._R8+  WPpZ05b1*x    + WPpZ05b2*x**2 + WPpZ05b3*x**3 )
 
        else
                dWsduou =                                                &  
     &           (   WPpZ10a1*x    + WPpZ10a2*x**2 + WPpZ10a3*x**3 ) /   &  
     &       ( 1._R8+  WPpZ10b1*x    + WPpZ10b2*x**2 + WPpZ10b3*x**3 )
 
        endif
      return
      endif
 
!     For Mu = -1 and 0 < u < 1   use Table V. on page 192.
      if ( Mu .le.-0.90_R8.and.  u .le. 1.00_R8) then
!
             x = u*u
!
             if ( Z  .le. 1.5_R8)  then
                dWsduou = WPmZ01a1*x                                     &  
     &                 + WPmZ01a2*x**2 + WPmZ01a3*x**3 + WPmZ01a4*x**4
 
        else if ( Z  .le. 3.0_R8)  then
                dWsduou = WPmZ02a1*x                                     &  
     &                 + WPmZ02a2*x**2 + WPmZ02a3*x**3 + WPmZ02a4*x**4
 
        else if ( Z  .le. 7.0_R8)  then
                dWsduou = WPmZ05a1*x                                     &  
     &                 + WPmZ05a2*x**2 + WPmZ05a3*x**3 + WPmZ05a4*x**4
 
        else
                dWsduou = WPmZ10a1*x                                     &  
     &                 + WPmZ10a2*x**2 + WPmZ10a3*x**3 + WPmZ10a4*x**4
 
        endif
      return
      endif
 
 
!     For Mu = -1 and 0 < u < 1   use Table V. on page 192.
      if ( Mu .le.-0.90_R8.and.  u .gt. 1.00_R8) then
!     .                                 If the velocity is too large,
!     .                                 damp the results at high u in such
!     .                                 a way that the value and derivative
!     .                                 are continuous at u=e01 (1), but the
!     .                                 values returned at u > e01 (1) are not
!     .                                 too large.  The table is not valid.
!     .                                 The hope is to induce E field
!     .                                 in TSC.
!
 
!     e01 etc:   power 01 etc of the u^2 (E-like) at which we expand
!     e1, e2:    first and second power of the expansion parameter u^2-e01
!     A0,01,02:  expansion coeficients around u^2=e01
             x = u*u
             iWhichWay = -1
             e01 =1._R8
             e02=e01*e01
             e03=e01*e02
             e04=e01*e03
             e1 =(x-e01)
             e2 =(x-e01)**2
!
             if ( Z  .le. 1.5_R8)  then
                A0 =    WPmZ01a1*e01 +    WPmZ01a2*e02                   &  
     &             +    WPmZ01a3*e03 +    WPmZ01a4*e04
                A1 =    WPmZ01a1     + 2._R8*WPmZ01a2*e01                &  
     &             + 3._R8*WPmZ01a3*e02 + 4._R8*WPmZ01a4*e03
                A2 =                      WPmZ01a2                       &  
     &             + 3._R8*WPmZ01a3*e01 + 6._R8*WPmZ01a4*e02
 
        else if ( Z  .le. 3.0_R8)  then
                A0 =    WPmZ02a1*e01 +    WPmZ02a2*e02                   &  
     &             +    WPmZ02a3*e03 +    WPmZ02a4*e04
                A1 =    WPmZ02a1     + 2._R8*WPmZ02a2*e01                &  
     &             + 3._R8*WPmZ02a3*e02 + 4._R8*WPmZ02a4*e03
                A2 =                      WPmZ02a2                       &  
     &             + 3._R8*WPmZ02a3*e01 + 6._R8*WPmZ02a4*e02
 
        else if ( Z  .le. 7.0_R8)  then
                A0 =    WPmZ05a1*e01 +    WPmZ05a2*e02                   &  
     &             +    WPmZ05a3*e03 +    WPmZ05a4*e04
                A1 =    WPmZ05a1     + 2._R8*WPmZ05a2*e01                &  
     &             + 3._R8*WPmZ05a3*e02 + 4._R8*WPmZ05a4*e03
                A2 =                      WPmZ05a2                       &  
     &             + 3._R8*WPmZ05a3*e01 + 6._R8*WPmZ05a4*e02
 
        else
                A0 =    WPmZ10a1*e01 +    WPmZ10a2*e02                   &  
     &             +    WPmZ10a3*e03 +    WPmZ10a4*e04
                A1 =    WPmZ10a1     + 2._R8*WPmZ10a2*e01                &  
     &             + 3._R8*WPmZ10a3*e02 + 4._R8*WPmZ10a4*e03
                A2 =                      WPmZ10a2                       &  
     &             + 3._R8*WPmZ10a3*e01 + 6._R8*WPmZ10a4*e02
 
 
        endif
        dWsduou = A0 + A1*e1 + A2*e2
      return
      endif
 
 
!     Unanticipated parameter range, set return to 0.00.
      dWsduou = 0.00_R8
      return
      END
!                                                                      |
!     Fits.F ends       -----------------------------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
