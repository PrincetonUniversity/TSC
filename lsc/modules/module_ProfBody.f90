      MODULE ProfBody
      USE PARAMS
      USE EMPARAMS
      IMPLICIT NONE
!     File: ProfBody.inc                        starts
!     Commentary on --Vec --Ary dVlVec dVlAry dVol:
!     --Vec refers to a vector on the TSC psi grid
!     --Ary refers to a vector on the LSC psi grid
!     dVlVec(j) is the volume between j-1 and j because vptemp
!     is a backward-calculated quantity.  iVlVec is the sum of dVlVec
!     dVol(j) is the volume centered on j
      REAL*8                                                             &  
     &        NeAry(NPSIDIM), PsiAry(NPSIDIM),  TeAry(NPSIDIM),          &  
     &        ZbrAry(NPSIDIM),iVlAry(NPSIDIM),                           &  
     &        dVol(NPSIDIM),  EdcAry(NPSIDIM), LnlAry(NPSIDIM),          &  
     &        BetZAry(NPSIDIM), MidAry(NPSIDIM)
 
!     COMMON /prcom/
!    ^        NeAry, PsiAry, TeAry,
!    ^        ZbrAry,iVlAry,
!    ^        dVol, EdcAry, LnlAry, BetZAry, MidAry
 
      REAL*8                                                             &  
     &        psiminx, psimaxx, delpsi, Te, Ne, Zbr, pe2min
!     COMMON / plascom /
!    ^        psiminx, psimaxx, delpsi, Te, Ne, Zbr, pe2min
 
!     File: ProfBody.inc                        ends
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE ProfBody
