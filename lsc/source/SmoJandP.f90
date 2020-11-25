!
!     ------------------------------------------------------------------
!
      SUBROUTINE SmoJandP(Radius)
      USE Jrf
      USE MKSetc
      USE params
      USE PIetc
      USE power_mod
      USE ProfBody
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip, idxFor0
      REAL*8    Radius
      REAL*8    NuSlow(NPSIDIM), Jdriven(NPSIDIM), Jdiffus(NPSIDIM),     &  
     &        Awk(NPSIDIM),    Bwk(NPSIDIM),     Cwk(NPSIDIM),           &  
     &        Jwk(NPSIDIM),    Xwk(NPSIDIM),     NparCube
      REAL*8    PrfNorm,         PrfOrig(NPSIDIM), PrfSmoo(NPSIDIM)
      EQUIVALENCE (PrfOrig(1),Jdriven(1))
      EQUIVALENCE (PrfSmoo(1),Jdiffus(1))
      PARAMETER (NparCube = 8.0_R8)
!
!     Do nothing if there is no diffusion; or no sensible radius
      if (DiffuJrf .le. 0.00_R8.or. Radius .le. 0.00_R8) return
 
!     Form the 'slowing down frequency'
!     \nu_slow = \ln \lambda n_e e^4 /(4 \pi \epsilon_0^2 m_e^2 c^3)
!                \times  n_\parallel^3
!
!     The assumed boudary conditions in SmoothJrf are that there is
!     a guard point on the low(central) side of the passed array:
!         > zero slope (in r-space) on inner side
!         > zero value on the outer side at idxFor0
!     However, idxFor0 must not be more than npsi+1
!
      do ip=1,npsi
        NuSlow(ip) =  NeAry(ip) * LnlAry(ip) *                           &  
     &               ECOULB**4 / ELECMS**2 * CLIGHT *                    &  
     &               4._R8*PI * 1.0E-14_R8* NparCube
      enddo
!
!     Smooth the current computed with the given field:
      do ip=1,npsi
         Jdriven(ip) = js(ip)
         if(Jdriven(ip) .ne. 0.00_R8) idxFor0=ip
      enddo
         if(idxFor0 .lt. npsi-1) idxFor0 = (npsi + idxFor0)/2
      call SmoothJrf      (Radius, DiffuJrf, npsi, idxFor0,              &  
     &                     NuSlow, Jdriven, Jdiffus,                     &  
     &                     Awk, Bwk, Cwk, Jwk, Xwk,                      &  
     &                    'PolFlux_Like', 'Zero2ndDeriv')
      do ip=1,npsi
         js(ip) = Jdiffus(ip)
      enddo
!
!     Smooth the current computed with the given field plus delta-field
      do ip=1,npsi
         Jdriven(ip) = jsp(ip)
      enddo
!
      call SmoothJrf      (Radius, DiffuJrf, npsi, idxFor0,              &  
     &                     NuSlow, Jdriven, Jdiffus,                     &  
     &                     Awk, Bwk, Cwk, Jwk, Xwk,                      &  
     &                    'PolFlux_Like', 'Zero2ndDeriv')
 
      do ip=1,npsi
         jsp(ip) = Jdiffus(ip)
      enddo
!
      if(PrfSpred .le. 0.00_R8.or. PrfSpred .gt. 1.00_R8) return
!     Smooth Power according to
!         Praytot(psi) = Prf-total n(psi) J-diffused(psi) dV(psi) /
!                        [ \Sigma  n(psi) J-diffused(psi) dV(psi) ] *
!                                                 PrfSpred     +
!                  Prf-ray-undiffused(psi)* (1. - PrfSpred)
!
!     Calculate normalization factor
         PrfNorm = 0.00_R8
      do ip=1,npsi
         PrfNorm     = abs(js(ip))*NeAry(ip)*dVol(ip) + PrfNorm
      enddo
         PrfNorm     = PraySum / PrfNorm
!
!     Calculate raw smoothed power
      do ip=1,npsi
         PrfOrig(ip) = PRaytot(ip)
         PrfSmoo(ip) = abs(js(ip))*NeAry(ip)*dVol(ip) * PrfNorm
      enddo
!
!     Fold raw smooth power with unsmoothed power
      do ip=1,npsi
         PrfSmoo(ip) = PrfSpred*PrfSmoo(ip) +                            &  
     &   (1._R8-PrfSpred)*PrfOrig(ip)
         PrayTot(ip) = PrfSmoo(ip)
      enddo
 
      return
      END
!
!     -----------------------------------------------------------------|
!
!     SUBROUTINE SmoothJrf(Radius, DiffCoef, NumJs, ....
!
!     Routine to Smooth Jrf
!     Copyright 1994 by
!     D. W. Ignat and S. C. Jardin
!     Plasma Physics Laboratory
!     Box 451
!     Princeton, New Jersey 08543
!
!     Smooth RF-driven current with a diffusion-like equation.
!     The idea is that the actual RF-driven current is modified from
!     the computed RF-driven current by slowing-down, characterized by
!     an inverse-time \nu_slow, and a cross-field diffusion,
!     characterized by DiffCoef.
!
!     Jdriven is given on points 1,2,...NumJs; regularly spaced
!     in poloidal flux. We add guard points at 0 and NumJs + 1.
!     We assume that poloidal flux is like radius squared, so that
!     the diffusion-like equation in radius
!
!     dJrf/dt = \nu_{slow} (Jdriven-Jrf) + 1/r d/dr( r DiffCoef dJrf/dr)
!
!     becomes
!
!     dJrf/dt = \nu_{slow} (Jdriven-Jrf) + d/dx(4 x DiffCoef dJrf/dx)/Radius^2
!
!     where  x  is normalized poloidal flux ranging from 0 to 1.
!
!     More correct use of generalised coordinates would follow this:
!         df/dt = 1/Jac  d/dp[ Jac (grad p)**2 D df/dp ]
!         where Jac is the Jacobian, or volume element factor.
!         Jac == (1/grad p) (\int^p dp/grad p)
!     In the case of p ~ r^2, grad p = 2 sqrt(p) and
!     \int^p dp/grad p== sqrt(p) so the J is 1/2, and cancels out.
!
!     We solve the equation for the steady state, from a difference
!     equation
!
!     Radius^2  *  \nu_slow
!     ----------------------  *  (Jrf - Jdriven) ==
!     4 DiffCoef (NumJs+1)^2
!
!    (xp+x)/2 * (Jrfp-Jrf) - (x+xm)/2 * (Jrf-Jrfm)
!
!     where xp is x-plus xm is x-minus etc and \delta x = 1/(NumJs+1)
!
!     The boundary condtions at index 0,1,2 are :
!        (1) that Jrf is constant (Jdiff_0 = Jdiff_1)              <or>
!        (2) that Jrf slope is constant (Jdiff_0 = 2 Jdiff_1 - Jdiff_2)
!
!     and at index idxFor0 that Jrf is zero.  Normally, one expects that
!     idxFor0 is NumJs+1, but maybe the current should not be zero at the edge-
!     it should be zero at idxFor0.
!
!     Boundary condition (1) goes naturally with radius-like coordinates;
!     Boundary condition (2) goes naturally with flux-like coordinates.
!
!     \nu_slow = \ln \lambda n_e e^4 /(4 \pi \epsilon_0^2 m_e^2 c^3)
!                \times  n_\parallel^3
!     where n_e is position dependent...but does not have to be.
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
