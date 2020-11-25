      MODULE SPDMOD
      USE PARAM
      IMPLICIT NONE
!
      INTEGER ncurrentg, numfluxloopsg
      REAL*8 :: atot,csum,vesppe,vescure,totcure,chicur,vescure2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: csuma, ppcur, fbcur
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pflux, pfluxo, eflux
!
      END MODULE SPDMOD
