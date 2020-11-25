!
      SUBROUTINE RyZnInit
      USE params
      USE RayBins
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,j
!
!     Initialize  ray zones                     -----------------------|
!                                                                      |
 
!
!                                       `Zero' arrays.
      j = iray
                do 10 i=1,NZONDIM
                    RofRay (i) = Rmaj
                    ZofRay (i) = 0.0_R8
                    PofRay (i) = 0.0_R8
                    NperRy (i) = 1.0_R8
                    rtPsRy (i) = 1.0_R8
                    PowrRy (i) = 1.0_R8
                    TimeRy (i) = 0.0_R8
                    DistRy (i) = 0.0_R8
                    DetrRy (i) = 0.0_R8
!               do 10 j=1,NRAYDIM     !! 9sep 93
                    ezsq (i,j) = 1.0_R8
                  dlnPds (i,j) = 0.0_R8
                  dlnPdsK(i,j) = 0.0_R8
                  dlnPdsX(i,j) = 0.0_R8
!                                      For izind, ivind the value 0 has
!                                      no meaning -- its an incorrect address.
!                                      This is used as a flag that this ray bin
!                                      was not filled in the calculation.
                    izind(i,j) = 0
                    ivind(i,j) = 0
 10             continue
!                                       But do npar separately since the
!                                       first zone location is filled with
!                                       information from launched spectrum.
                do 20 i=2,nzones
!               do 20 j=1,nrays      !! 9sep 93
                    npar (i,j) = 1.0_R8
 20             continue
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
