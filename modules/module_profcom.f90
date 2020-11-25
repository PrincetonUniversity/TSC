      MODULE PROFCOM
      USE PARAM
      IMPLICIT NONE

      INTEGER, PARAMETER :: pnplt=16
      INTEGER :: pnplts
      INTEGER :: np2sav
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  npsitsv
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  ihds

      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xplt, yplt
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  scary
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  psix
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  timesv
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  ymaxsf, yminsf
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xcs, ycs
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  div

      END MODULE PROFCOM


