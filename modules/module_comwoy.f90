      MODULE COMWOY
      USE PARAM
      IMPLICIT NONE
!  woyke common blocks used in subroutines:
!      regler, conag, sort, cirequ1, efield und main
 
!       common block zur variation der reglerparameter
!               und zur begrenzung der integrierer (anti - reset - windup)
!
!       gain    reglerverstaerkungen und entkopplung k, k', e
!
!       nsz     nachstellzeiten und vorhaltzeiten tn, tv
!
!       arw     minimal- und maximalwerte zur begrenzung der integrierer
!
!       dtabt   abtastzeit der echtzeitregelung
!
!       wlwdoc, wlwnext beschraenkung der ausgabe auf 200 punkte
        REAL*8  fb1vk1 (6), fb1vk2 (6)
        REAL*8  fb1ck1 (2), fb1ck2 (2)
        REAL*8  fb1ohk1, fb1ohk2
        REAL*8  fb2ck1 (2)
        REAL*8  fb1ve (6,6), fb1ce (2,2)
        REAL*8  dtabt, wlwdoc, wlwnext
        REAL*8  uohkp
 
 
        REAL*8  fb1vtn (6), fb1ctn (2), fb1ohtn, fb2ctn (2)
        REAL*8  fb1vtv (6), fb1ctv (2), fb1ohtv
 
 
        REAL*8  ivmax (6), ivmin (6)
        REAL*8  icmax (2), icmin (2)
        REAL*8  ucmax (2), ucmin (2)
        REAL*8  iohmax, iohmin
        REAL*8  uohmax, uohmin
        integer n35
 
 
!   common block zur steuerung mittels "preprogrammed waveforms"
!
!   wavef   sollgroessen, welche von der experimentsteuerung definiert
!           werden.
!
!
      REAL*8 gsoll (6)
      REAL*8 ivff  (10), iffint (10)
      REAL*8 icff  (2)
      REAL*8 ipsoll
      REAL*8 iohff
 
!     alte werte und ableitungen fuer tracking
      REAL*8 gsollm1 (6), gsollp (6), ipsollm1, ipsollp, icgso (2)
 
!     integer nwoy  !  file nummer fuer test output
!     integer nwoyo !  file nummer fuer restart control output
!     integer nwoyi !  file nummer fuer restart control input
      integer nwoy,nwoyo,nwoyi
      REAL*8 gist(6), ipist, icist (2), ivist (6), iohist
      REAL*8 ivsoll(7),iohsoll,icsoll (2), ucsoll(2)
!
!        ipow    stromprogramm der poloidalfeldspulen
!        ipowm1  zeitpunkt n - 1
!        ipowm2  zeitpunkt n - 2
!        ipowp   ableitung zum zeitpunkt n - 1
!        flxp    induzierter fluss innerer kreise
!        vegint  zustaende der integrierer fuer stromregler der vertikalfeldspul
 
!       REAL*8 ipow (pngroup), ipowm1 (pngroup), ipowm2 (pngroup),
!    1     ipowp (pngroup),
!    &         flxp (pngroup), vegint (pngroup),
        REAL*8                                                           &  
     &     fb2cei (2), fb2cep (2),                                       &  
     &         fb2cd (2), fb2cx (2)
       REAL*8, ALLOCATABLE, DIMENSION(:) :: ipow, ipowm1, ipowm2,        &  
     &                                      ipowp, flxp, vegint
 
!======================================================================
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE COMWOY
