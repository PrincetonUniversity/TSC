!#include "f77_dcomplx.h"
!     Obviously, any fudge factor can be thrown in on the \nu_slow or on
!     the DiffCoef owing to the heuristic derivation of this treatment.
!     The extention to spatially varying DiffCoef is straightforward.
!
!     A modification for an input Jrf based on a radius-like grid
!     is very simple:
!
!     Radius^2  *  \nu_slow
!     ----------------------  *  r * (Jrf - Jdriven) ==
!     1 DiffCoef (NumJs+1)^2
!
!    (rp+r)/2 * (Jrfp-Jrf) - (r+rm)/2 * (Jrf-Jrfm)
!
!     where rp is r-plus rm is r-minus etc and \delta r = 1/(NumJs+1)
!
!     This spirit of this approach is similar to that of
!     V. Fuchs, I. P. Shkarofsky, R. A. Cairns, and P. T. Bonoli,
!     ``Simulations of lower hybrid current drive and ohmic transformer
!     recharge,'' Nuclear Fusion {\bf 29} 1479 1989.
!     It appears that the practical difference is that
!     in the Fuchs work a loss-rate at the edge is specified,
!     rather than the DiffCoef, and then the DiffCoef
!     is found by a shooting method.
!
!     In contrast, the present routine specifies the DiffCoef and finds
!     the modified (reduced) current.
!     (and the current loss by implication)
!
!     The following is a heuristic derivation of the formula,
!     loosely based on the Fuchs paper:                  (rectangular coords!)
!
!    df/dt = d/dv ( Dq df/dv )  + d/dv(Dc df/dv + \nu v f) + d/dr(D df/dr)
!            RF (q-l) souce term  Collisional terms        Spatial  Diffusion
!
!   multiply by e v and integrate over v,
!   integrating by parts in first 3 terms on RHS
!
!    dJ/dt = \int e(-Dq df/dv) dv + \int e(-Dc df/dv) dv - \int \nu e v f dv
!            RF q-l source term    Collisional diffusion  Dynamical friction
!
!           +                                                d/dr(D dJ/dr)
!                                                           Spatial Diffusion
!
!   In the dynamical friction term, pull \nu out of the integral
!   making it \nu_eff \cdot {fudge-factor}, but call it \nu still.
!   Then you have - \nu J for that term
!
!   Tell yourself that if RF is there, the q-l source overwhelms the
!   collisional velocity diffusion, which we now ignore.  Then
!
!    dJ/dt = \int e(-Dq df/dv) dv  -\nu J + d/dr( D dJ/dr )
!
!   Tell yourself that the rf souce term can be written as
!   \nu' times the RF current found in the normal LSC calculation J_o;
!   and then tell yourself that we might as well treat \nu' == \nu
!
!   Then
!
!     dJ/dt = \nu (J_o - J) +  d/dr ( D dJ/dr )
!
!   This has some desirable properties, like if
!   D is small then the J moves to J_o on the \nu time scale
!   and if D is large then J and J_o can be quite different.
!   Also if \nu is space-dependent and D not then the J will
!   tend to move out (as Giruzzi finds) because there are
!   fewer collisions there.
!
!   It does not do anything like move the fast particles around
!   and change the damping because of moved fast particles.
!
!   The Edc treatment based on Karney Fisch is a confused
!   thing after this.
!
!
      SUBROUTINE SmoothJrf(Radius, DiffCoef, NumJs, idxFor0,             &  
     &                     NuSlow, Jdriven, Jdiffus,                     &  
     &                     Awk, Bwk, Cwk, Jwk, Xwk,                      &  
     &                     IndeptVar, BoundCond)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER NumJs, idxFor0, idxFor0m1
      REAL*8    Radius, DiffCoef
      REAL*8    NuSlow(NumJs), Jdriven(NumJs), Jdiffus(NumJs)
      REAL*8    Awk(NumJs), Bwk(NumJs), Cwk(NumJs), Jwk(NumJs)
      REAL*8    Xwk(NumJs), xwk0
      REAL*8    TimeDum, Jp1
      INTEGER i, ierror
      CHARACTER*(*) IndeptVar, BoundCond
      REAL*8    BC1, BC2
      REAL*8 AREAL
!
      Jp1 = AREAL(NumJs+1)
 
      if (idxFor0 .lt. NumJs) then
         idxFor0m1 = idxFor0 - 1
      else
         idxFor0m1 = NumJs
      endif
!
!     Default Radial Independent Variable Type is:  PolFlux_Like
!     This is Area_Like.
!
      if (IndeptVar .eq. 'RadDist_Like') then
         TimeDum = Radius**2/(1._R8*DiffCoef) / Jp1**2
      else if (IndeptVar .eq. 'PolFlux_Like')then
         TimeDum = Radius**2/(4._R8*DiffCoef) / Jp1**2
      else
         TimeDum = Radius**2/(4._R8*DiffCoef) / Jp1**2
      endif
!
!     Zero slope boundary condition at i=0;i=1 or
!     Zero second deriv boundary condition at i=0;i=1;i=2
      if(BoundCond .eq. 'Zero1stDeriv') then
         BC1 = 1.0_R8
         BC2 = 0.0_R8
      else if(BoundCond .eq. 'Zero2ndDeriv') then
         BC1 = 2.0_R8
         BC2 =-1.0_R8
      else
         BC1 = 1.0_R8
         BC2 = 0.0_R8
      endif
!
!     Xwk is filled with the average value of x(i) and x(i+1);
!     a more important purpose is workspace for the matrix inverter.
!     Xwk does not have to be filled above index idxFor0m1, but we
!     fill it anyway to emphasize that xwk goes with the geometry of
!     the problem, not with where we think the current must go to zero.
!
      do i=1,NumJs
         xwk(i) = (AREAL(i)+0.5_R8)/ Jp1
      enddo
         xwk0 = 0.5_R8/Jp1
!
      Awk(1)     = 0.0_R8
      Bwk(1)     = TimeDum*NuSlow(1)
      if(IndeptVar .eq. 'RadDist_Like') Bwk(1) = Bwk(1)*1._R8/Jp1
      Bwk(1)     = Bwk(1) + xwk0*(1._R8-BC1) + xwk(1)
      Cwk(1)     = - (xwk(1) + BC2*xwk0 )
!
!     Awk(1) and Cwk(NumJs) are never used;  setting to zero is cosmetic
!
      do i=2,idxFor0m1
         Awk(i) = - xwk(i-1)
         Bwk(i) = TimeDum*NuSlow(i)
         if(IndeptVar .eq. 'RadDist_Like') Bwk(i) = Bwk(i)*AREAL(i)/Jp1
         Bwk(i) = Bwk(i) + xwk(i-1) + xwk(i)
         Cwk(i) = - xwk(i)
      enddo
!
!     Zero value boundary condition at i=idxFor0
      Cwk(idxFor0m1) = 0.0_R8
!
!     Awk(1) and Cwk(idxFor0-1) are never used;  setting to zero is cosmetic
!
      do i=1,idxFor0m1
         Jwk(i) = TimeDum*NuSlow(i)*Jdriven(i)
         if(IndeptVar .eq. 'RadDist_Like') Jwk(i) = Jwk(i)*AREAL(i)/Jp1
      enddo
!
!
!     Find the 'diffused' current;
!     Note that we are done with Xwk, and now pass it to tridiaNR for
!     it to use as workspace.
!
      call tridiaNR(Awk,Bwk,Cwk,Jwk,Jdiffus,idxFor0m1,                   &  
     &                                        Xwk,ierror)
!     Fill remaining Jdiffused locations with zeros, if necessary.
 
      if (idxFor0 .le. NumJs) then
         do i=idxFor0, NumJs
            Jdiffus(i) = 0.00_R8
         enddo
      endif
 
      if (ierror .ne. 0) then
         do i=1,NumJs
            Jdiffus(i) = Jdriven(i)
         enddo
      endif
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
