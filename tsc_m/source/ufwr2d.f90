!#include "f77_dcomplx.h"
      subroutine ufwr2d
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!...writes 2d (y vs. rho at time t) UFILES        C. Kessel 10/2000
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER kv,numone,numv,i,nuf2d,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tfluxb,gps,sqrgps
      REAL*8 tflux,ellip,delta,x1,x2,z1,z2,position
      REAL*8 AREAL
!============
      dimension kv(50)
!     dimension usave(50,ppsi),urho(ppsi),kv(50),ugrho(ppsi),
!    +ugrho2(ppsi),urmaj(ppsi),urmin(ppsi),ukappa(ppsi),udelta(ppsi),
!    +usurf(ppsi),uprad(ppsi),ujtot(ppsi),ubpol(ppsi),upoh(ppsi),
!    +uindent(ppsi)
      character*8 ulab(pglobs)
!
      data numone/1/
      data numv/41/
!     data kv(1)/ 1 /
      data (kv(i),i=1,41)/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,        &  
     & 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,           &  
     & 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,           &  
     & 41/
      data (ulab(i),i=1,41)/                                             &  
     &'RMAJOR  ','RMINOR  ','KAPPAR  ','DELTAR  ','INDENTR ',            &  
     &'VOLUME  ','GRHO1   ','GRHO2   ','SURF    ','TE      ',            &  
     &'TI      ','NE      ','SNBIE   ','SNBII   ','QNBIE   ',            &  
     &'QNBII   ','CURNBI  ','QICRHE  ','QICRHI  ','CURICRH ',            &  
     &'QLHE    ','QLHI    ','CURLH   ','QFUSE   ','QFUSI   ',            &  
     &'QOHM    ','QRAD    ','QEI     ','CURTOT  ','NIMP    ',            &  
     &'NM1     ','NM2     ','NM3     ','ZEFFR   ','Q       ',            &  
     &'BPOL    ','CHIE    ','CHII    ','QECHE   ','QECHI   ',            &  
     &'CURECH  '/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: usave
      REAL*8, ALLOCATABLE, DIMENSION(:) :: urho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ugrho
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ugrho2
      REAL*8, ALLOCATABLE, DIMENSION(:) :: urmaj
      REAL*8, ALLOCATABLE, DIMENSION(:) :: urmin
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ukappa
      REAL*8, ALLOCATABLE, DIMENSION(:) :: udelta
      REAL*8, ALLOCATABLE, DIMENSION(:) :: usurf
      REAL*8, ALLOCATABLE, DIMENSION(:) :: uprad
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ujtot
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ubpol
      REAL*8, ALLOCATABLE, DIMENSION(:) :: upoh
      REAL*8, ALLOCATABLE, DIMENSION(:) :: uindent

      IF(.not.ALLOCATED(usave)) ALLOCATE( usave(50,ppsi), STAT=istat)
      IF(.not.ALLOCATED(urho)) ALLOCATE( urho(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ugrho)) ALLOCATE( ugrho(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ugrho2)) ALLOCATE( ugrho2(ppsi), STAT=istat)
      IF(.not.ALLOCATED(urmaj)) ALLOCATE( urmaj(ppsi), STAT=istat)
      IF(.not.ALLOCATED(urmin)) ALLOCATE( urmin(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ukappa)) ALLOCATE( ukappa(ppsi), STAT=istat)
      IF(.not.ALLOCATED(udelta)) ALLOCATE( udelta(ppsi), STAT=istat)
      IF(.not.ALLOCATED(usurf)) ALLOCATE( usurf(ppsi), STAT=istat)
      IF(.not.ALLOCATED(uprad)) ALLOCATE( uprad(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ujtot)) ALLOCATE( ujtot(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ubpol)) ALLOCATE( ubpol(ppsi), STAT=istat)
      IF(.not.ALLOCATED(upoh)) ALLOCATE( upoh(ppsi), STAT=istat)
      IF(.not.ALLOCATED(uindent)) ALLOCATE( uindent(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : ufwr2d  ' 
!============      
!
!     call flxvol(2)
!
      tfluxb = AREAL(npsit-1)*dpsi
      do 1000 i=2,npsit
      gps = gxmja2(i)/xmja2(i)*(tpi*qprof2(i))**2
      sqrgps = sqrt(gps)
      tflux = AREAL(i-1)*dpsi
      urho(i) = sqrt(tflux/tfluxb)
      ugrho(i) = sqrgps/(2._R8*sqrt(tflux*tfluxb))
      ugrho2(i) = ugrho(i)**2
      call shapex(ellip,delta,x1,x2,z1,z2,xsv2(i))
      urmaj(i)=0.5_R8*(x1+x2)
      urmin(i)=0.5_R8*(x1-x2)
      ukappa(i)=ellip
      udelta(i)=delta
      uindent(i)=0._R8
      usurf(i)=2._R8*pi*urmaj(i)*2._R8*pi*urmin(i)                       &  
     &*sqrt((1._R8+ukappa(i)**2)/2._R8)
      uprad(i)=(sradion(i)*(usdp/usdt)+savebre(i)+savecyc(i)+            &  
     &saveimp(i))
      ujtot(i)=(gxmja2(i)-gxmja2(i-1))*udsi*xplas/vp(i)
      upoh(i)=(0.5_R8*(as(i-1)+as(i))+vlooph(i))/tpi                     &  
     &*(gxmja2(i)-gxmja2(i-1))
 1000 continue
      urmaj(1)=xmag
      urmin(1)=0.0_R8
      ukappa(1)=ukappa(2)+(ukappa(3)-ukappa(2))/(xsv2(3)-xsv2(2))        &  
     &*(psimin-xsv2(2))
      udelta(1)=0.0_R8
      usurf(1)=0.0_R8
      uprad(1)=uprad(2)
      ujtot(1)=ujtot(2)
      upoh(1)=upoh(2)
!
      do 2000 i=1,npsit
      usave(1,i) = urmaj(i)
      usave(2,i) = urmin(i)
      usave(3,i) = ukappa(i)
      usave(4,i) = udelta(i)
      usave(5,i) = uindent(i)
      usave(6,i) = vary(i)
      usave(7,i) = ugrho(i)
      usave(8,i) = ugrho2(i)
      usave(9,i) = usurf(i)
      usave(10,i) = te(i)
      usave(11,i) = ti(i)
      usave(12,i) = ane(i)
!..beam elec fueling
      usave(13,i) = 0.0_R8
!..beam ion fueling
      usave(14,i) = 0.0_R8
!..beam elec heating
      usave(15,i) = 0.0_R8
!..beam ion heating...savei
      usave(16,i) = 0.0_R8
!..beam current drive
      usave(17,i) = 0.0_R8
!..ICRH elec heating...savefw
      usave(18,i) = 0.0_R8
!..ICRH ion heating
      usave(19,i) = savei(i)*udsp/udst
!..ICRH current drive
      usave(20,i) = 0.0_R8
!..LH elec heating...savelh
      usave(21,i) = 0.0_R8
!..LH ion heating...savilh
      usave(22,i) = 0.0_R8
!..LH current drive
      usave(23,i) = 0.0_R8
!..alpha elec heating
      usave(24,i) = savea(i)*udsp/udst
!..alpha ion heating
      usave(25,i) = savia(i)*udsp/udst
!..ohmic elec heating
      usave(26,i) = upoh(i)*udsp/udst
!..radiated power
      usave(27,i) = uprad(i)*udsp/udst
!..equipartition term
      usave(28,i) = 0.0_R8
!..total toroidal current density
      usave(29,i) = ujtot(i)
!..impurity density
      usave(30,i) = 0.0_R8
!..additional species
      usave(31,i) = 0.0_R8
      usave(32,i) = 0.0_R8
      usave(33,i) = 0.0_R8
!..Zeff
      usave(34,i) = zeffa(i)
!..safety factor
      usave(35,i) = qprof2(i)
!..poloidal field
      usave(36,i) = 0.0_R8
!..electron thermal diffusivity
      usave(37,i) = chiesec(i)
!..ion thermal diffusivity
      usave(38,i) = chiisec(i)
!..ECH elec heating
      usave(39,i) = 0.0_R8
!..ECH ion heating
      usave(40,i) = 0.0_R8
!..ECH current drive
      usave(41,i) = 0.0_R8
 2000 continue
!
      nuf2d=9
      if( numargs .lt. 1 ) then
         filename = 'tok_99999TSC_2d.dat'
      else
         filename = 'tok_99999TSC_2d.dat' // '.' // trim(suffix)
      end if
      open(nuf2d,file=trim(filename),status='unknown',                   &  
     &position='append')
      call timedate(itime,idate,imach,isuffix)
!
      do 5000 k=1,numv
!     write(nuf2d,949)
      write(nuf2d,948) ulab(kv(k))
      write(nuf2d,950) idate
      write(nuf2d,951)
      write(nuf2d,952)
      write(nuf2d,953)
      write(nuf2d,954)
      write(nuf2d,955) ulab(kv(k))
      write(nuf2d,956)
      write(nuf2d,957) npsit
      write(nuf2d,958) numone
      write(nuf2d,1001) (urho(i),i=1,npsit)
      write(nuf2d,1002) tpest(itpest)
      write(nuf2d,1001) (usave(k,i),i=1,npsit)
      write(nuf2d,959)
      write(nuf2d,949)
 5000 continue
!
 1001 format(6(1pe13.6))
 1002 format(1pe13.6)
 949  format('**************************************************')
 948  format(';-FILE NAME : UF99999.TSC#',a8)
 950  format('  99999TSC  2 0 6',14x,';-SHOT #- F(X,Y) DATA -UF2DWR- ',  &  
     &a8)
 951  format(' TSC Simulation',16x,';-SHOT DATE-  UFILES ASCII FILE',    &  
     &1x,'SYSTEM')
 952  format('   0',27x,';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-')
 953  format(' RHO                   ',8x,';-INDEPENDENT VARIABLE        &  
     & LABEL:                                                            &  
     & X-')
 954  format(' TIME               SEC',8x,';-INDEPENDENT VARIABLE        &  
     & LABEL:                                                            &  
     & Y-')
 955  format(1x,a8,22x,';-DEPENDENT VARIABLE LABEL-')
 956  format(' 0',29x,';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM.')
 957  format(7x,i3,21x,';-# OF X PTS-')
 958  format(7x,i3,21x,';-# OF Y PTS-  X,Y,F(X,Y) DATA FOLLOW:')
 959  format(' ;----END-OF-DATA-----------------COMMENTS:-----------')
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
