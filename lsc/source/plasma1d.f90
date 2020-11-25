 
!
      SUBROUTINE plasma1d (psi, psiold, RBphi,Tee,pe2,pi2,aio,ael)
      USE params
      USE ProfBody
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8    psderiv(0:2, 0:2), psimat(0:3, 0:3)
      COMMON / plxxr / psderiv, psimat
      INTEGER jRBphi, jpe2, jTee, jpi2, jAio, jAel, jVpr,                &  
     &     RECALC, first_call
      PARAMETER(RECALC = -1)
      DATA jRBphi, jpe2, jTee, jpi2, jAio, jAel, jVpr  /                 &  
     &     RECALC, RECALC, RECALC,                                       &  
     &     RECALC, RECALC, RECALC, RECALC /
      DATA first_call / TRUE /
      REAL*8                                                             &  
     &     RBphipr, pe2pr, Teepr, pi2pr, Aiopr, Aelpr,                   &  
     & psl,psi, psiold, RBphi, Tee, pe2, pi2, aio, ael,                  &  
     &     RBphio, pe2old, Teeold, pi2old, Aioold, Aelold
      DATA RBphio, pe2old, Teeold, pi2old, Aioold, Aelold /              &  
     &     0,       0,      0,      0,       0,     0    /
      if(first_call .eq. TRUE)then
         first_call = FALSE
         Psiold = 0._R8
         RBphio = 0._R8
         pe2old = 0._R8
         Teeold = 0._R8
         pi2old = 0._R8
         Aioold = 0._R8
         Aelold = 0._R8
      endif
 
      if(psi .eq. psiold)then
         RBphi  = RBphio
         pe2    = pe2old
         Tee    = Teeold
         pi2    = pi2old
         Aio    = Aioold
         Ael    = Aelold
         return
      endif
 
      call grnu1d(NpsiJ, PsiVec, RBpVec, jRBphi, RBPhiCoefs,             &  
     &     psi, RBphi, RBphipr)
!
!     enforce bounds on psi !!!!!!!!!!!!!!
      psi = max(psiminx, psi)
!     psi = min(psimaxx, psi)
      psl = min(psimaxx, psi)
!     stop forcing bounds on psi on the maximum
!     enforce bounds on psi !!!!!!!!!!!!!!
!
      call grnu1d(NpsiJ, PsiVec, pe2Vec, jpe2  , pe2Coefs,               &  
     &     psl, pe2  , pe2pr)
      call grnu1d(NpsiJ, PsiVec, TeeVec, jTee  , TeeCoefs,               &  
     &     psl, Tee  , Teepr)
      call grnu1d(NpsiJ, PsiVec, pi2Vec, jpi2  , pi2Coefs,               &  
     &     psl, pi2  , pi2pr)
      call grnu1d(NpsiJ, PsiVec, AioVec, jAio  , AioCoefs,               &  
     &     psl, Aio  , Aiopr)
      call grnu1d(NpsiJ, PsiVec, AelVec, jAel  , AelCoefs,               &  
     &     psl, Ael  , Aelpr)
 
      Psiold = Psi
      RBphio = RBphi
      pe2old = pe2
      Teeold = Tee
      pi2old = pi2
      Aioold = Aio
      Aelold = Ael
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
