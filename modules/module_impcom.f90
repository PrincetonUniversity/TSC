      MODULE IMPCOM
      USE PARAM
      IMPLICIT NONE
!......pimp-----no of impurity species
!......pne------no of density locations in radiation arrays
!......pte------no of temp locations in radiation arrays
!......pchrg----charge state in radiation arrays (block data)
!......pchrgmx--maximum charge state for impurity tranpsort (nucz+1)
!......ppsi--maximum number of radial zones
!
!============
! idecl:  explicitize implicit INTEGER declarations:
!============
! idecl:  explicitize implicit REAL declarations:
!============
!     integer pne,pte,pcgOx,pcgC,pcgFe,pcgBe,pcgNe,pcgKr,
!    1        pcgAr,pcgW,pimp,pchrgmx,ppsi
!     parameter(pne=5,pte=25,pcgOx=9,pcgC=7,pcgFe=27,pcgBe=5,pcgNe=11)
!     parameter(pcgKr=37,pcgAr=19,pcgW=75)
!     parameter(pimp=09,pchrgmx=80)
!     parameter(ppsi=500)
!
      REAL*8 :: begcsp
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  nq
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  ainz, rec
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  rad
      REAL*8 :: fincsp
!
      REAL*8 :: nqo
!     ------------------------------------------------------------------
      REAL*8 :: begcrt
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alinzOx
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  alradOx, alrecOx
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alinzC 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  alradC , alrecC
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alinzFe
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  alradFe, alrecFe
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alinzBe
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  alradBe, alrecBe
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alinzNe
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  alradNe, alrecNe
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alinzKr  
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  alradKr, alrecKr
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alinzAr
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  alradAr, alrecAr
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alinzW
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  alradW , alrecW
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  altei
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  alnei
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  nne, nte, nchrgsr
      REAL*8 :: fincrt
!     ------------------------------------------------------------------
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  snp, aimp, bimp, cimp, dimp,  &  
     &                              ework, fwork 
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE IMPCOM
