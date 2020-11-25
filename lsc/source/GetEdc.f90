!
!     ------------------------------------------------------------------
!
      SUBROUTINE GetEdc (DifAmt)
      USE params
      USE PIetc
      USE ProfBody
      USE RayBins
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER idxPsi
      REAL*8    DifAmt, psi, Edc
      REAL*8    MinEdc
      DATA    MinEdc    / 0.0001_R8/
!     DifAmt is intended to be 0.0 when the E_dc is desired as given
!     and intended to be something like 0.0001 when we are trying to form
!     the derivative   d ln J / d ln E
!
!
      do 10 idxPsi = 1, npsi
        psi = PsiAry(idxPsi)
        call linr1d(NpsiJ, PsiVec, EdcVec,                               &  
     &       psi, Edc)
        Edc = Edc + DifAmt
        if (abs(Edc) .le. MinEdc) then
          if (Edc .ge. 0._R8) then
            Edc = +MinEdc
          else
            Edc = -MinEdc
          endif
        endif
        EdcAry(idxPsi) =  Edc
 10   continue
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
