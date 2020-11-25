!
!     ValInit begins  -------------------------------------------------|
!     -                                                                |
!     -                                                                |
!
      SUBROUTINE ValInit
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      call MiscInit
!                                       doesn't fit elsewhere
      call MkGrids
!                                       set up velocity, psi grids
!                                       needed for quasilinear calculations
!     call WrGrids
!                                       output same ! not after June 2000
      call ProfInit
!                                       interpolate Ne, Ni, Te, Ti
!                                       to PsiAry from TSC data points
      call FeInit
!                                       compute constants needed for
!                                       solution of Fe ....
      call DqlInit
!                                       and for solution of Dql
!     call WrDcoll
!                                       output same ! not after June 2000
      call DqlClear
!                                       set Dql == 0
      return
      END
!
!                                                                      |
!     Openings, Initializations end         ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
