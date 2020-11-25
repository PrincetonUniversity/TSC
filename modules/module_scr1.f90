      MODULE SCR1
      USE PARAM
      USE SCRATCH
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 :: sumpel
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  piz,wiz,giz,uiz,aiz,pim,wim,  &  
     &                                 gim,uim,                          &  
     &                          aim,pip,wip,gip,uip,abip,biz,bim,bip,    &  
     &                          rim,riz,rip,qim,qiz,qip
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  abouna,abounb,abounc,abound,  &  
     &                          qt2   , qt4  , qt   ,  fac7,             &  
     &                          term3a, term4, term5, term6,             &  
     &                          term8 , rt   , term9, fac7a
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  wf1pa , gf1pa, uf1pa, bf1pa,  &  
     &                          rf1pa , qf1pa
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  wf2pa , gf2pa, uf2pa, bf2pa,  &  
     &                          rf2pa , qf2pa
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  wf1ta , gf1ta, uf1ta, bf1ta,  &  
     &                          rf1ta , qf1ta
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  dpara , dparb, dparc, dpard
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  ajp2,utsave,rjcdg
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  gapr,gapi,gapv
!
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  pling1,pling2,pling3,      &  
     &                                 pling4,                           &  
     &                              gling1,gling2,gling3,pling5,         &  
     &                              ajp3,ajp4
!
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR1
