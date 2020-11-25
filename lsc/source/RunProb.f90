!     Fits.F begins     -----------------------------------------------|
!                                                                      |
!
!     Reference: ``Current in wave-driven plasmas,'' by
!     Charles F. F. Karney and Nathaniel J. Fisch,
!     Phys. Fluids {\bf 29} 180-192 (1986)
!
!     Fits.F written March 1991 by D. W. Ignat.  All DATA statements
!     and logic based on the reference.
!     Copyright D. W. Ignat, E. J. Valeo, C. F. F. Karney, N. J. Fisch,
!     Princeton University, Plasma Physics Laboratory
!
!     These are memonics for the numerical constants:
!     R         Run-away probability
!     W         Energy imparted to the electric field by stopped electrons
!     WP        dW/du /u
!     p         mu is plus:  +1
!     m         mu is minus: -1
!     Z         effective charge
!
!     These are other usages:
!     RunProb   Run-away probability, or the R(u,mu) function of reference.
!     WsloDwn   Energy imparted to electric field by an electron as it
!               slows down, or the W_s(u,mu) function of reference.
!     WsloPrm   \partial W_s/\partial u / u function of reference.
!               This is the ratio P_{\rm el} to P_{\rm in} , or ratio of
!               power coupled from the rf source into electromagnetic energy
!               to the rf power absobed by the electrons.
!               Very roughly, Pel/Pin rises from 0 to 1 as u_par goes from
!               0 to 5, with most rise as u_par goes from 0 to 3.
!               See Fig. 7a of the reference.
!
!     u         velocity normalized to run-away velocity (magnitude!)
!     Mu        direction cosine, v-par/v-total
!     Z         Z_effective
!     x         u*u, or (u-1) in RunProb
!
!
!     Review of definitions and conventions of the reference:
!     q         carries sign of electron charge == - e
!     E         is parallel to B, and in the positive direction
!     v_par     is positive in the direction of E and B
!     \Gamma == \ln \Lambda \frac{n q^4}{4\pi \epsilon_0^2 m^2}
!     v_r    == -sign(qE) \sqrt{m\Gamma / \abs{qE} } (a positive number)
!     Dreicer Velocity v_D = - \sqrt{2 + Z} v_r
!     u         v /(\abs{ v_r})
!     u_\par    v_\par / v_r
!     u_\perp   v_\perp / \abs{v_r}
!
!     The following numbered lines are a copy of the MACSYMA fit equation,
!     kept here for redundancy.
!01  WsloDwn = Mu1 * U04 / (5.+Z)
!02 ^   -    (2. + Z + 3.*Mu2)                                * U06 /
!03 ^                                     ( 3.*(3.+Z)*(5.+Z) )
!04 ^   + 2.*( (24.+19.*Z+3.*Z**2)*Mu1 + (9.+Z)*Mu3 )         * U08 /
!05 ^                       ( (3.+Z)*(5.+Z)*(7.+3.*Z)*(9.+Z) )
!06 ^  -(
!07 ^    (1041.+1864.*Z + 1189.*Z**2 + 316.*Z**3 + 30.*Z**4 ) * Mu0
!08 ^   +( 417.+ 497.*Z +  181.*Z**2 +  21.*Z**3) * 10.       * Mu2
!09 ^   +  (9.+Z)*(13.+3.*Z)                      *  5.       * Mu4
!10 ^   )                                                     * U10 /
!11 ^      ( 5.*(2.+Z)*(3.+Z)*(5.+Z)*(7.+3.*Z)*(9.+Z)*(13.+3.*Z) )
!
!
!     -----------------------------------------------------------------|-------
!
      REAL*8 FUNCTION RunProb ( uGiven , MuGiven , ZGiven )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8                                                             &  
     &        u , Mu , Z , uGiven, MuGiven, ZGiven,                      &  
     &        x
      REAL*8                                                             &  
     &     RpZ01a0,  RpZ01a1,  RpZ01a2,  RpZ01a3,  RpZ01b2,  RpZ01b3 ,   &  
     &     RpZ02a0,  RpZ02a1,  RpZ02a2,  RpZ02a3,  RpZ02b2,  RpZ02b3 ,   &  
     &     RpZ05a0,  RpZ05a1,  RpZ05a2,  RpZ05a3,  RpZ05b2,  RpZ05b3 ,   &  
     &     RpZ10a0,  RpZ10a1,  RpZ10a2,  RpZ10a3,  RpZ10b2,  RpZ10b3
 
      DATA                                                               &  
     &     RpZ01a0,  RpZ01a1,  RpZ01a2,  RpZ01a3,  RpZ01b2,  RpZ01b3 /   &  
     &    -3.68063_R8,  4.23913_R8, -4.55894_R8, -0.39755_R8,            &  
     &    -1.22774_R8,  1.41450_R8/
      DATA                                                               &  
     &     RpZ02a0,  RpZ02a1,  RpZ02a2,  RpZ02a3,  RpZ02b2,  RpZ02b3 /   &  
     &    -4.97636_R8,-16.09015_R8,  0.83188_R8,  0.21737_R8,            &  
     &    6.84615_R8, -0.98649_R8/
      DATA                                                               &  
     &     RpZ05a0,  RpZ05a1,  RpZ05a2,  RpZ05a3,  RpZ05b2,  RpZ05b3 /   &  
     &    -4.27687_R8, -4.33629_R8,  0.30338_R8,  0.05697_R8,            &  
     &    3.21315_R8, -0.47749_R8/
      DATA                                                               &  
     &     RpZ10a0,  RpZ10a1,  RpZ10a2,  RpZ10a3,  RpZ10b2,  RpZ10b3 /   &  
     &    -4.94597_R8, -1.53482_R8,  0.10112_R8,  0.03087_R8,            &  
     &     2.45288_R8, -0.36896_R8/
 
!
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
!
!     .                                 Electrons of velocity less than
!     .                                 the runaway velocity simply do not
!     .                                 run away:
      if ( u*u .le. 1._R8) then
        RunProb = 0._R8
        return
      endif
!
!     Backward-running electrons run away easily once u > 1, and it
!     does not depend much on Z.  However, there is no fit given in the
!     paper, even though a graph is given.  This is MY crude fit to
!     the Fig. 2 of the reference.
      if ( Mu .le. -0.90_R8) then
        x = 0.85_R8* (10._R8/Z)**2 * ( (u*u-1._R8)/3._R8)**2
        if (x .gt. 0.85_R8)                                              &  
     &                  x = 1.00_R8- 0.15_R8*exp(-(x-0.85_R8)**2)
        RunProb = x
        return
      endif
 
!     No fit given, no graph given for -1 < Mu < 1 , so I set RunProb to 0.0.
      if ( Mu .gt. -0.90_R8.and. Mu .lt. 0.90_R8) then
        RunProb = 0.00_R8
        return
      endif
 
!     So we are left with u > 1, and Mu = 1. (or at least Mu > 0.90)
!     This uses the fit data from Table I., page 190.
      x = (u - 1._R8)
      if        ( Z .le. 1.5_R8) then
        RunProb = exp (                                                  &  
     &  ( RpZ01a0 + RpZ01a1*x + RpZ01a2*x**2 + RpZ01a3*x**3 ) /          &  
     &  (               1.0_R8*x + RpZ01b2*x**2 + RpZ01b3*x**3 )         &  
     &                )
        else if ( Z .le. 3.0_R8) then
        RunProb = exp (                                                  &  
     &  ( RpZ02a0 + RpZ02a1*x + RpZ02a2*x**2 + RpZ02a3*x**3 ) /          &  
     &  (               1.0_R8*x + RpZ02b2*x**2 + RpZ02b3*x**3 )         &  
     &                )
        else if ( Z .le. 7.0_R8) then
        RunProb = exp (                                                  &  
     &  ( RpZ05a0 + RpZ05a1*x + RpZ05a2*x**2 + RpZ05a3*x**3 ) /          &  
     &  (               1.0_R8*x + RpZ05b2*x**2 + RpZ05b3*x**3 )         &  
     &                )
        else if ( Z .le. 10._R8) then
        RunProb = exp (                                                  &  
     &  ( RpZ10a0 + RpZ10a1*x + RpZ10a2*x**2 + RpZ10a3*x**3 ) /          &  
     &  (               1.0_R8*x + RpZ10b2*x**2 + RpZ10b3*x**3 )         &  
     &                )
      endif
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
