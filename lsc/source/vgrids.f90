!     ------------------------------------------------------------------
      SUBROUTINE vgrids
      USE FeBins
      USE params
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL ugrid
!     generate parallel velocity grid
      if (2 * (nv / 2) .eq. nv) then
         call LSCwarn (' require ODD nv ')
         nv = nv - 1
      endif
!                                       set UniformGRID into Vpar
         call ugrid(Vpar, nv, Vmin, Vmax)
!                                       set EXPonentialGRID into Vpar
!        call egrid(Vpar, nv, Vmin, Vmax, StpRange)
!
      call mkdelv(nv, dvsym, Vpar)
!                                       compute forward difference on Vpar,
!                                       and place into dvsym
      call mkdvp(nv, dvplus, Vpar)
      call mkvinv(nv, dvip, dvplus)
      call mkvinv(nv, dvisym, dvsym)
!                                       construct inverse, so dvisym is 1/dv
      IvZero = (nv + 1) / 2
!                                       require nv to be odd !!!!!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
