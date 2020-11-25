 
!     PlasEJV   begins                      ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE ProfInit
      USE params
      USE ProfBody
      USE RayBins
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER                                                            &  
     &        ip
      do 10 ip = 1, npsi
         call GePlPar(PsiAry(ip), NeAry(ip), TeAry(ip), ZbrAry(ip),      &  
     &                LnlAry(ip), BetZAry(ip))
!                                       get plasma parameters
!                                       given PsiAry, returns :
!                                       TeAry   electron temperature in Kev
!                                       NeAry   electron density in cm^-3
!                                       ZbrAry  Z bar for e-i collison rate
!                                       LnlAry  Coulombic logarithm
!                                       BetZAry Beta value of Eder & Valeo
!                                       Lnl BetZ added by D. P. Enright
 
 10   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
