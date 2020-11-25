      block data impdata
!
!******************************************************************************


!
      USE CLINAM
      USE RADTAB
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
! * * * impurity radiation tables
!           nte = number of electron temperatures in table
!           nne = number of electron densities in table
!           altei = alog10 of interpolation temperatures
!           alnei = alog10 of interpolation densities
!           alinzr = alog10 of ionization rate
!           alrecr = alog10 of recombination rate
!           alradr = alog10 of radiative emission rate
!

!============
      end block data impdata
