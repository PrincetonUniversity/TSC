 
      SUBROUTINE GePlPar(psi, Negl, Tegl, Zbrgl, Lnlgl, BetZgl)
!                                       get plasma parameters
      USE params
      USE ProfBody
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8                                                             &  
     &        psi, RBphi,      pe2, pi2,                                 &  
     &        aio, ael, Negl, Tegl, Zbrgl, Lnlgl, BetZgl,  psiold
      DATA psiold / -100._R8/
      call plasma1d (psi, psiold, RBphi, Tegl, pe2, pi2, aio, ael)
      Negl = 1.E14_R8* pe2 / pe2Fac
      Zbrgl = pe2 / pi2 * (pi2Fac / pe2Fac)
      if (Tegl .le. .01_R8) then
         Lnlgl = 23._R8- log(sqrt(Negl) * 1._R8/(Tegl * sqrt(Tegl) *     &  
     &                                     1000._R8**(3._R8/2._R8)))
      else
         Lnlgl = 24._R8- log(sqrt(Negl) * 1._R8/(Tegl * 1000._R8))
      endif
      BetZgl = (1._R8+ Zbrgl)/5._R8
      return
      END
!                                                                      |
!                                                                      |
!     PlasEJV   ends                        ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
