      MODULE SAPROP
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8 :: begcsa
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  ane,zeffa,avez,simpe,te,ti,  &  
     &                                 sradion,                          &  
     &                           zeffa2,                                 &  
     &                           ajavlh,ajpary,ajavfw,ajavlh2,djdelh2,          &
     &                           djdetsc2,voltlp,anne,animp,sradiono,    &  
     &                           aneimp,                                 &  
     &                           avezimpa,amassimpa,amasshyda,aimassa,   &  
     &                           sumnia,sumnqa,pden,ajavbs,ajavcd,anhy,  &  
     &                           anox,anca,anfe,anbe,powtsc,currtsc,     &  
     &                           djdetsc,djdelh,sravejet, ajavec
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  sinzrec
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  powtsci
      REAL*8 :: tpsir,tpsil,tpsiri,tpsili,sjane0, sumjet
      REAL*8 :: fincsa
!     ------------------------------------------------------------------
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SAPROP
