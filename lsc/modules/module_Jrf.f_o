      MODULE Jrf
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
      INTEGER, PRIVATE, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
c     Jrf.inc ---------------------------------------------------------
c     -                                                               |
c     -                                                               |
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
c
!     COMMON /JrfCom0 /
!    ^        vnorm, nuRuna, vnmax, gmrun, muminus, muplus, dEdcAmnt,
!    ^        DiffuJrf, PrfSpred,
!    ^        js,          jsp,        Jray,
!    ^        nRunDot,     jRunDot,    vRunIdx,
!    ^        IrIntgrl,          IpIntgrl,
!    ^        ugr,
!    ^        vnormPos,          vnormNeg,
!    ^        VparMaxP,          VparMaxN
c
!     COMMON /JrfCom1 / vnormNOK, ivrun
c
c     js        j_stopped.  See Karney and Fisch paper
c               current density deposited by the rf power
c     jsp       j_stopped if Edc were Edc + dE [= Edc + dEdcAmnt]
c               the p is for plus
c     Jray      js resolved by velocity bin and by psi bin;
c               needed for a graph requested by BERNABEI October 94
c               which is implemented in JrfDiagn
c
c     IrIntgrl  I from rays inside indexed psi surface
c               found by Intgrl js
c     IpIntgrl  I as above but at increased E; Intgrl jsp
c     DiffuJrf  Arbitrary diffusion (smoothing) of Jrf in m^2/sec
c     PrfSpred  in range 0.0 --> 1.0 controls spreding of Prf to match
c               shape of n_e(psi) * J_{rf-smoothed}(psi)
c     ivrun     index of the velocity at which the v=v_runaway for
c               cooperative current drive; if zero, then does not happen
c     nRunDot   time rate of change of runaway population cm-3/sec
c     jRunDot   time rate of change of runaway current, from density
c               increase
c               both above evaluated at the runaway velocity
c     vRunIdx   floating value of the integer iv index of the runaway
c               velocity
c     -                                                               |
c     -                                                               |
c     -----------------------------------------------------------------
 
 
 
 
 
 
 
 
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE Jrf
