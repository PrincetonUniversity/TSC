!#include "f77_dcomplx.h"
      subroutine ufwr1d
!
!
      USE CLINAM
      USE SCR14
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!...writes 1D (y vs. time) UFILES keeping total number
!...of points around 100        C. Kessel 10/2000
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER kv,numv,i,k,nuf1d,njump,ncou
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ajump
      REAL*8 AREAL
!============
      dimension kv(250)
      character*8 ulab(pglobs)
!
!...prescribe the number of variables to be written to UFILES in NUMV
!...and put number in KV array that corresponds to number of
!...desired variable in ULAB array
!
      data numv/17/
      data (kv(i),i=1,17)/ 9 , 35 , 34 , 37 , 36 , 48 , 119 ,            &  
     & 135 , 77 ,173 , 80 , 27 , 28 , 61 , 47 , 141 , 46/
!
!...those variable labels in all CAPITALS correspond to standard
!...UFILE variables
!
      data (ulab(i),i=1,40)/                                             &  
     &'  amach ','   ekin ','   ekinp','   ekint','  iplim ',            &  
     &'ZMAG    ','RMAG    ','   cur  ','IP      ','delp*tpi',            &  
     &'pmin*tpi',' diamag ','VSURF   ','  qzero ','  qedge ',            &  
     &'    dt  ','beta-tor','<n>/nmur','1/q nr/b','li/2+bp ',            &  
     &'BEPMHD  ',' li/2   ',' li vs q','taue-kg ','taue(ms)',' taue2  ',  &  
!    &                                                                   &  
     &'TE0     ','TI0     ','te/te-av','ti/ti-av',' chiauxs',            &  
     &'chiohms ','hflux-mw','RGEO    ','AMIN    ',                       &  
     &'DELTA   ','KAPPA   ',' xsep(1)',' zsep(1)','vsec-int'/
      data (ulab(k),k=41,96)/                                            &  
     &'resv-sec','vsec-abs','loopv-oh','vsec-oh ','ptot(mw)',            &  
     &'PFUSION ','POHM    ','PICRH   ','nullapol','nullapol',            &  
     &'  dipole','  dipole','quadrupl','quadrupl','hexapole',            &  
     &'hexapole','octapole','octapole','decapole','decapole',            &  
     &'Q95     ','  qstar ','  qcyl  ',' delta95',' ellip95',            &  
     &' xsepcal',' zsepcal','  psep  ',' psepcal',' n/ncrit',            &  
     &' dt1s   ',' dt2s   ',' dt3s   ','1/qcylin',' beta-b0',            &  
     &' beta-0 ','NEL     ','n-vol av','density ','WTH     ',            &  
     &'ellip-90',' delt-90','VOL     ',' no part','  eps1a ',            &  
     &'  eps1c ','  eps10 ','  eps2a ','  eps2c ','  eps20 ',            &  
     &'  eps3a ','  eps3c ','  eps30 ','  eps4a ','  eps4c ',            &  
     &'  eps40 '/
      data (ulab(k),k=97,109)/                                           &  
     &'thm ener','pf ener ','tf ener ','hflux en','aux ener',            &  
     &'int c en','pf bn en','tf bn en','total en','  tevv  ',            &  
     &'  ffac  ',' npsit  ',' resid  '/
!
      data ulab(119)/'PRAD    '/
      data ulab(122)/'ploss   '/
      data (ulab(k),k=125,126)/'cpu(min)','cpu-rate'/
      data (ulab(k),k=129,130)/'alphdens','alphbeta'/
      data (ulab(k),k=134,139)/'alphfrac','ZEFF    ',                    &  
     &                         ' ni-ave ',' nhe-ave',                    &  
     &                         ' W-dot  ',' Taue/G '/
      data ulab(141)/'IBOOT   '/
      data (ulab(k),k=147,148)/' xsep(2)',' zsep(2)'/
      data (ulab(k),k=151,152) / ' r(q=1) ',' li-(ga)'/
      data (ulab(k),k=164,165) / ' thalo  ',' whalo  '/
      data (ulab(k),k=166,167) / ' Pohm-H ',' Pohm-P '/
      data (ulab(k),k=168,170) / 'VS-Poynt','        ','Ejima Co'/
      data (ulab(k),k=173,174) /'LI      ','   Ct   '/
      data (ulab(k),k=175,177) / 'Ptot    ',' Palpha ','frac-T  '/
      data (ulab(k),k=200,204) / 'dndtpel ','xpel    ',                  &  
     & 'tepel   ','totradp ','Radpel  '/
      data (ulab(k),k=205,206) /'pellatom','tot atom'/
      data (ulab(k),k=207,212) / 'iqtrubmx','nskipsf ',                  &  
     & 'runawayc','plas_cur','avalanch','resource' /
      data (ulab(k),k=240,241)/' pecrh(2)',' eccd(2)'/
!
      nuf1d=10
      if ( numargs .lt. 1 ) then
         filename = 'tok_99999TSC_1d.dat'
      else
         filename = 'tok_99999TSC_1d.dat' // '.' // trim(suffix)
      end if
      open(nuf1d,file=trim(filename),status='unknown')
      call timedate(itime,idate,imach,isuffix)
!
      if(nrecord .le. 100) then
      njump=1
      go to 998
      endif
      ajump = AREAL(nrecord)/100._R8
      njump = int(ajump)
      if( (ajump-AREAL(njump)) .gt. 0.5_R8) then
      njump = njump+1
      endif
!
 998  continue
      ncou = 0
      do 999 i=1,nrecord,njump
      ncou = ncou+1
 999  continue
!
      do 5000 k=1,numv
!     write(nuf1d,949)
      write(nuf1d,948) ulab(kv(k))
      write(nuf1d,950) idate
      write(nuf1d,951)
      write(nuf1d,952)
      write(nuf1d,953)
      write(nuf1d,954) ulab(kv(k))
      write(nuf1d,955)
      write(nuf1d,956) ncou
      call disk(kv(k)+1,1,big1,big2)
      write(nuf1d,1001) (big2(i),i=1,nrecord,njump)
      if(kv(k) .eq. 46 .or. kv(k) .eq. 47 .or. kv(k) .eq. 48             &  
     &.or. kv(k) .eq. 80) then
      write(nuf1d,1001) (big1(i)*1.E6_R8,i=1,nrecord,njump)
      go to 1111
      endif
      if(kv(k) .eq. 77) then
      write(nuf1d,1001) (big1(i)*1.E20_R8,i=1,nrecord,njump)
      go to 1111
      endif
      write(nuf1d,1001) (big1(i),i=1,nrecord,njump)
 1111 write(nuf1d,957)
      write(nuf1d,949)
 5000 continue
!
 1001 format(6(1pe13.6))
 949  format('**************************************************')
 948  format(';-FILE NAME : UF99999.TSC#',a8)
 950  format('  99999TSC  1 0 6',14x,';-SHOT #- F(X,Y) DATA -UF1DWR- ',  &  
     &a8)
 951  format(' TSC Simulation',16x,';-SHOT DATE-  UFILES ASCII FILE',    &  
     &1x,'SYSTEM')
 952  format('   0',27x,';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-')
 953  format(' TIME               SEC',8x,';-INDEPENDENT VARIABLE        &  
     & LABEL:                                                            &  
     & X-')
 954  format(1x,a8,22x,';-DEPENDENT VARIABLE LABEL-')
 955  format(' 0',29x,';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM.')
 956  format(7x,i3,21x,';-# OF PTS-  X, F(X) DATA FOLLOW:')
 957  format(' ;----END-OF-DATA-----------------COMMENTS:-----------')
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
