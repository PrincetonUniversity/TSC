      MODULE Jrf
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
      INTEGER, PRIVATE, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     Jrf.inc ---------------------------------------------------------
!     -                                                               |
!     -                                                               |
      INTEGER vnormNOK(NPSIDIM), ivrun
      REAL*8    vnorm, nuRuna, vnmax, gmrun, muminus, muplus, dEdcAmnt
      REAL*8    DiffuJrf, PrfSpred
      REAL*8    js(NPSIDIM), jsp(NPSIDIM), Jray(NVELDIM,NPSIDIM)
      REAL*8    nRunDot(NPSIDIM), jRunDot(NPSIDIM), vRunIdx(NPSIDIM)
      REAL*8    IrIntgrl(NPSIDIM), IpIntgrl(NPSIDIM)
      REAL*8    ugr(NVELDIM)
      REAL*8    vnormPos(NPSIDIM), vnormNeg(NPSIDIM)
      REAL*8    VparMaxP(NPSIDIM), VparMaxN(NPSIDIM)
      DATA    vnmax, dEdcAmnt / 0.99_R8, 0.0001_R8/
      DATA    DiffuJrf, PrfSpred / 0.000_R8, 0.000_R8/
!
!     COMMON /JrfCom0 /
!    ^        vnorm, nuRuna, vnmax, gmrun, muminus, muplus, dEdcAmnt,
!    ^        DiffuJrf, PrfSpred,
!    ^        js,          jsp,        Jray,
!    ^        nRunDot,     jRunDot,    vRunIdx,
!    ^        IrIntgrl,          IpIntgrl,
!    ^        ugr,
!    ^        vnormPos,          vnormNeg,
!    ^        VparMaxP,          VparMaxN
!
!     COMMON /JrfCom1 / vnormNOK, ivrun
!
!     js        j_stopped.  See Karney and Fisch paper
!               current density deposited by the rf power
!     jsp       j_stopped if Edc were Edc + dE [= Edc + dEdcAmnt]
!               the p is for plus
!     Jray      js resolved by velocity bin and by psi bin;
!               needed for a graph requested by BERNABEI October 94
!               which is implemented in JrfDiagn
!
!     IrIntgrl  I from rays inside indexed psi surface
!               found by Intgrl js
!     IpIntgrl  I as above but at increased E; Intgrl jsp
!     DiffuJrf  Arbitrary diffusion (smoothing) of Jrf in m^2/sec
!     PrfSpred  in range 0.0 --> 1.0 controls spreding of Prf to match
!               shape of n_e(psi) * J_{rf-smoothed}(psi)
!     ivrun     index of the velocity at which the v=v_runaway for
!               cooperative current drive; if zero, then does not happen
!     nRunDot   time rate of change of runaway population cm-3/sec
!     jRunDot   time rate of change of runaway current, from density
!               increase
!               both above evaluated at the runaway velocity
!     vRunIdx   floating value of the integer iv index of the runaway
!               velocity
!     -                                                               |
!     -                                                               |
!     -----------------------------------------------------------------
 
 
 
 
 
 
 
 
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE Jrf
