!
!     -----------------------------------------------------------------|-------
!
      REAL*8 FUNCTION WsloDwn ( uGiven , MuGiven , ZGiven )
!     WsloDwn   is the inergy in units of m v_runaway^2 imparted to the
!               electric field by an electron as it slows down.
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8                                                             &  
     &        u , Mu , Z , uGiven, MuGiven, ZGiven,                      &  
     &        x
      REAL*8                                                             &  
     &      Mu0 , Mu1 , Mu2 , Mu3 , Mu4 , U02, U04, U06, U08, U10
 
      REAL*8                                                             &  
     &     WpZ01a2,  WpZ01a3,  WpZ01a4,  WpZ01b1,  WpZ01b2,  WpZ01b3 ,   &  
     &     WpZ02a2,  WpZ02a3,  WpZ02a4,  WpZ02b1,  WpZ02b2,  WpZ02b3 ,   &  
     &     WpZ05a2,  WpZ05a3,  WpZ05a4,  WpZ05b1,  WpZ05b2,  WpZ05b3 ,   &  
     &     WpZ10a2,  WpZ10a3,  WpZ10a4,  WpZ10b1,  WpZ10b2,  WpZ10b3
      DATA                                                               &  
     &     WpZ01a2,  WpZ01a3,  WpZ01a4,  WpZ01b1,  WpZ01b2,  WpZ01b3 /   &  
     &     0.16612_R8, -0.01495_R8,  0.00775_R8,  0.37136_R8,            &  
     &     0.02240_R8,  0.01645_R8/
      DATA                                                               &  
     &     WpZ02a2,  WpZ02a3,  WpZ02a4,  WpZ02b1,  WpZ02b2,  WpZ02b3 /   &  
     &     0.14200_R8, -0.04048_R8,  0.01145_R8,  0.12253_R8,            &  
     &     0.00384_R8,  0.02440_R8/
      DATA                                                               &  
     &     WpZ05a2,  WpZ05a3,  WpZ05a4,  WpZ05b1,  WpZ05b2,  WpZ05b3 /   &  
     &     0.09880_R8, -0.05152_R8,  0.01113_R8, -0.19484_R8,            &  
     &     0.00559_R8,  0.02362_R8/
      DATA                                                               &  
     &     WpZ10a2,  WpZ10a3,  WpZ10a4,  WpZ10b1,  WpZ10b2,  WpZ10b3 /   &  
     &     0.06537_R8, -0.03895_R8,  0.00738_R8, -0.32456_R8,            &  
     &     0.02797_R8,  0.01526_R8/
 
 
      REAL*8                                                             &  
     &     WmZ01a2,  WmZ01a3,  WmZ01a4,  WmZ01a5 ,                       &  
     &     WmZ02a2,  WmZ02a3,  WmZ02a4,  WmZ02a5 ,                       &  
     &     WmZ05a2,  WmZ05a3,  WmZ05a4,  WmZ05a5 ,                       &  
     &     WmZ10a2,  WmZ10a3,  WmZ10a4,  WmZ10a5
      DATA                                                               &  
     &     WmZ01a2,  WmZ01a3,  WmZ01a4,  WmZ01a5 /                       &  
     &    -0.16483_R8, -0.13420_R8,  0.15346_R8, -0.24314_R8/
      DATA                                                               &  
     &     WmZ02a2,  WmZ02a3,  WmZ02a4,  WmZ02a5 /                       &  
     &    -0.14186_R8, -0.09297_R8,  0.06661_R8, -0.12870_R8/
      DATA                                                               &  
     &     WmZ05a2,  WmZ05a3,  WmZ05a4,  WmZ05a5 /                       &  
     &    -0.09975_R8, -0.04781_R8,  0.00606_R8, -0.03545_R8/
      DATA                                                               &  
     &     WmZ10a2,  WmZ10a3,  WmZ10a4,  WmZ10a5 /                       &  
     &    -0.06651_R8, -0.02797_R8, -0.00247_R8, -0.00934_R8/
!
!
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
!
!
!     If Mu < 0.5, then use the MACSYMA-derived formula of page 191.
      U02 = u*u
      if ( U02 .le. 0.25_R8) then
 
        Mu0 = 1._R8
        Mu1 = Mu
        Mu2 = Mu*Mu
        Mu3 = Mu*Mu2
        Mu4 = Mu2*Mu2
 
        U04 = U02*U02
        U06 = U02*U04
        U08 = U04*U04
        U10 = U04*U06
!       --------------------------------------------------------------|
        WsloDwn = Mu1 * U04 / (5._R8+Z)                                  &  
     &   -    (2._R8+ Z + 3._R8*Mu2)      * U06 /                        &  
     &                                     ( 3._R8*(3._R8+Z)*(5._R8+Z) )  &  
     &   + 2._R8*( (24._R8+19._R8*Z+3._R8*Z**2)*Mu1 + (9._R8+Z)*Mu3 )    &  
     &        * U08 /                                                    &  
     &                       ( (3._R8+Z)*(5._R8+Z)*(7._R8+3._R8*Z)       &  
     &   *(9._R8+Z) )                                                    &  
     &  -(                                                               &  
     &    (1041._R8+1864._R8*Z + 1189._R8*Z**2 + 316._R8*Z**3 +          &  
     &    30._R8*Z**4 ) * Mu0                                            &  
     &   +( 417._R8+ 497._R8*Z +  181._R8*Z**2 +  21._R8*Z**3) *         &  
     &   10._R8* Mu2                                                     &  
     &   +  (9._R8+Z)*(13._R8+3._R8*Z)                  *  5._R8* Mu4    &  
     &   )                                                     * U10 /   &  
     &      ( 5._R8*(2._R8+Z)*(3._R8+Z)*(5._R8+Z)*(7._R8+3._R8*Z)*       &  
     &   (9._R8+Z)*(13._R8+3._R8*Z) )
 
!       --------------------------------------------------------------|
 
        return
      endif
 
!     So if u > 0.5 but -1 < Mu < 1 , then set WsloDwn to 0.00
!     The paper gives no fits for this region, but Fig. 4 does graph it.
      if ( Mu .gt. -0.90_R8.and. Mu .lt. 0.90_R8) then
        WsloDwn = 0.00_R8
        return
      endif
 
 
      x = u*u
 
!     Mu is +1, so use the Table II. on page 190.
!     This is valid for u < 5.
      if ( Mu .ge. 0.90_R8.and. u .le. 5.00_R8) then
             if ( Z  .le. 1.5_R8)  then
                WsloDwn =                                                &  
     &                 ( WpZ01a2*x**2 + WpZ01a3*x**3 + Wpz01a4*x**4 ) /  &  
     &            ( 1._R8+ WpZ01b1*x    + WpZ01b2*x**2 + WpZ01b3*x**3 )
        else if ( Z  .le. 3.0_R8)  then
                WsloDwn =                                                &  
     &                 ( WpZ02a2*x**2 + WpZ02a3*x**3 + Wpz02a4*x**4 ) /  &  
     &            ( 1._R8+ WpZ02b1*x    + WpZ02b2*x**2 + WpZ02b3*x**3 )
 
        else if ( Z  .le. 7.0_R8)  then
                WsloDwn =                                                &  
     &                 ( WpZ05a2*x**2 + WpZ05a3*x**3 + Wpz05a4*x**4 ) /  &  
     &            ( 1._R8+ WpZ05b1*x    + WpZ05b2*x**2 + WpZ05b3*x**3 )
 
        else
                WsloDwn =                                                &  
     &                 ( WpZ10a2*x**2 + WpZ10a3*x**3 + Wpz10a4*x**4 ) /  &  
     &            ( 1._R8+ WpZ10b1*x    + WpZ10b2*x**2 + WpZ10b3*x**3 )
        endif
      return
      endif
 
!     If Mu = -1 and 0 < u < 1 then use the fit of Table III. on page 191.
      if ( Mu .le.-0.90_R8.and. u*u .le. 1.00_R8) then
             if ( Z  .le. 1.5_R8)  then
                WsloDwn = WmZ01a2*x**2                                   &  
     &                 + WmZ01a3*x**3 + WmZ01a4*x**4 + WmZ01a5*x**5
        else if ( Z  .le. 3.0_R8)  then
                WsloDwn = WmZ02a2*x**2                                   &  
     &                 + WmZ02a3*x**3 + WmZ02a4*x**4 + WmZ02a5*x**5
        else if ( Z  .le. 7.0_R8)  then
                WsloDwn = WmZ05a2*x**2                                   &  
     &                 + WmZ05a3*x**3 + WmZ05a4*x**4 + WmZ05a5*x**5
        else
                WsloDwn = WmZ10a2*x**2                                   &  
     &                 + WmZ10a3*x**3 + WmZ10a4*x**4 + WmZ10a5*x**5
        endif
      return
      endif
 
!     Unanticipated input parameter, no data; set return to 0.00.
      WsloDwn = 0.00_R8
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
