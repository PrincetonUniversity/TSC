      MODULE SCR14
      USE PARAM
      IMPLICIT NONE
!============
! idecl:  explicitize implicit REAL declarations:
!============
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  big1,big2,big3,big4,big5
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  temp
      CHARACTER*8, DIMENSION(pglobs) ::  nlab
!     INTEGER, ALLOCATABLE, DIMENSION(:) ::  ibig1
!     equivalence (big1(1),ibig1(1))

      INTEGER, PRIVATE :: i, k

      data (nlab(i),i=1,40)/                                             &  
     &'  amach ','   ekin ','   ekinp','   ekint','  iplim ',            &  
     &'  zmag  ','  xmag  ','   cur  ','   cur  ','delp*tpi',            &  
     &'pmin*tpi',' diamag ','surfvolt','  qzero ','  qedge ',            &  
     &'    dt  ','beta-tor','<n>/nmur','1/q nr/b','li/2+bp ',            &  
     &'betapol ',' li/2   ',' li vs q','taue-kg ','taue(ms)',            &  
     &' taue2  ',' te(0)  ',' ti(0)  ','te/te-av','ti/ti-av',            &  
     &' chiauxs','chiohms ','hflux-mw','   r0   ','minorrad',            &  
     &'delt-tri','ellip   ',' xsep(1)',' zsep(1)','vsec-int'/
      data (nlab(k),k=41,96)/                                            &  
     &'resv-sec','vsec-abs','loopv-oh','vsec-oh ','ptot(mw)',            &  
     &'palpha  ','pohmic  ','paux    ','nullapol','nullapol',            &  
     &'  dipole','  dipole','quadrupl','quadrupl','hexapole',            &  
     &'hexapole','octapole','octapole','decapole','decapole',            &  
     &'  q95   ','  qstar ','  qcyl  ',' delta95',' ellip95',            &  
     &' xsepcal',' zsepcal','  psep  ',' psepcal',' n/ncrit',            &  
     &' dt1s   ',' dt2s   ',' dt3s   ','1/qcylin',' beta-b0',            &  
     &' beta-0 ','n-lineav','n-vol av','density ','int ener',            &  
     &'ellip-90',' delt-90',' volume ',' no part','  eps1a ',            &  
     &'  eps1c ','  eps10 ','  eps2a ','  eps2c ','  eps20 ',            &  
     &'  eps3a ','  eps3c ','  eps30 ','  eps4a ','  eps4c ',            &  
     &'  eps40 '/
      data (nlab(k),k=97,109)/                                           &  
     &'thm ener','pf ener ','tf ener ','hflux en','aux ener',            &  
     &'int c en','pf bn en','tf bn en','total en','  tevv  ',            &  
     &'  ffac  ',' npsit  ',' resid  '/
!
      data (nlab(k),k=125,126)/'cpu(min)','cpu-rate'/
      data (nlab(k),k=129,130)/'alphdens','alphbeta'/
      data (nlab(k),k=134,139)/'alphfrac',' zeff   ',                    &  
     &                         ' ni-ave ',' nhe-ave',                    &  
     &                         ' W-dot  ','  H_98  '/
      data (nlab(k),k=147,148)/' xsep(2)',' zsep(2)'/
      data (nlab(k),k=151,152) / ' r(q=1) ',' li-(ga)'/
      data (nlab(k),k=164,165) / ' thalo  ',' whalo  '/
      data (nlab(k),k=166,167) / ' Pohm-H ',' Pohm-P '/
      data (nlab(k),k=168,170) / 'VS-Poynt','        ','Ejima Co'/
      data (nlab(k),k=173,174) /'   li   ','   Ct   '/
      data (nlab(k),k=175,177) / 'Ptot    ',' Palpha ','frac-T  '/
      data (nlab(k),k=200,204) / 'dndtpel ','xpel    ',                  &  
     & 'tepel   ','totradp ','Radpel  '/
      data (nlab(k),k=205,206) /'pellatom','tot atom'/
      data (nlab(k),k=207,212) / 'iqtrubmx','nskipsf ',                  &  
     & 'runawayc','plas_cur','avalanch','resource' /
      data (nlab(k),k=218,220) / 'Te(ped) ','Ti(ped) ','Inj Curr'/
      data (nlab(k),k=221,221) / 'Injcur 2'/
      data (nlab(k),k=222,222) / '-dW core'/
      data (nlab(k),k=224,224) / '-dW tot '/
      data (nlab(k),k=237,239) / 'P&H Area','P+Halo A','HALOAREA' /

! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE SCR14
