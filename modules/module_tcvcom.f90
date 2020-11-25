      MODULE TCVCOM
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 :: begtcv
      REAL*8 :: apla,aplp,caltime,dzel,dissi,epsza,epszo,epsqa,epsqo
      INTEGER :: icur,iprtcv,ielch,mmax,nelz,nprob,nprtcv,nelch,nshac,   &  
     &           mshap,nluup,nob,nfob,nbxob,nbzob,nshap,nvvel,nvsob
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  ax
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  aph
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  a7
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  a8
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  a8l
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  berr, bee, coth
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  bxob, bzob
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  coesv
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  cerr, curr
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  cx, cxold, cxave, cxold2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  dex
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  dep
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  dsh
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  ddd, dddd
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  elss
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  elsv
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: elvvi, elvv
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  eee, ffob
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  ferr, flxa, flxo, flxd
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  gx
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  gp
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  hshap
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  hvvel
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  isvf, isvl
!
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  rhoa, rhoo, rho
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  rscurr
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  rshap, wshap, zshap
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  rvvel, rhn, wvvel, zvvel
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  sith, thb, vvv, xb, zb
      REAL*8 :: tmstp,timtcv2,timesnw,timesld,time1,time2,uc,vssh2,vssh,  &  
     &          vssha,vssh0,xip,xop,xvsob,zup,zlp,zvsob
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  volttcv
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  umeg
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  www, xf, zf
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xfob, zfob
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  xbxob, xbzob, zbxob, zbzob
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  xniv
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  zniv,xhnn,zhnn,xhne,xhoo,    &  
     &                                 zhoo,                             &  
     &                            xhee,xhss,zhss,xhse
      REAL*8 :: fintcv
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE TCVCOM
