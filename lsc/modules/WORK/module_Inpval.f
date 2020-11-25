C     NameList begins                       ---------------------------|
C                                                                      |
C                                                                      |
!============
!============
      NAMELIST /inpvalue/
     ^            fghz,   ntors,    npols,     nrays,    nGrps,
     ^         nparmax,  nparmin,  npolmin,  npolmax,
     ^         centers, couplers,   widths,   powers, phaseDeg,
     ^          DoBram,  nslices,
     ^        DiffuJrf, PrfSpred
c
      NAMELIST /inpexprt/
     ^          HstpLH,    nstep,    nfreq,     npsi,   nzones,
     ^              nv,    nsmoo,     nsmw,
     ^         nRampUp,    nFlat, WeghtItr,
     ^           PlFlg,   PrFlg ,    idiag,
     ^          DoXcam,   Do1Rpr,   Do0Edc,
     ^        TailTeps, TailPeps, TailNeps,
     ^        ScatKdeg, TurnNegs,   DoTRAN,
     ^           thet0
C                                                                      |
C                                                                      |
C     NameList ends                         ---------------------------|
 
c
c EXAMPLE FILE
c  &inpvalue
c  fghz = 4.6,
c  nGrps = 1,
c  powers(1) = 1.00, powers(2) = 1.00,
c  phaseDeg(1) = 130., phaseDeg(2) = 130., nslices=301,
c  ntors=15, npols=1,
c  DoBram = 1,
c  couplers(1)='TFTRLHCD',
c  couplers(2)='TFTRLHCD',
c  DiffuJrf = 0.00,
c  /
c
c  &inpexprt
c  HstpLH = .01,  nstep = 20000,  nfreq = 50,  npsi = 100,
c  nzones = 2000,
c  nv = 401, nsmoo = 9, nsmw = 3,  weghtitr = .2,
c  plflg( 1) = 0, 1, 0,
c  plflg( 4) = 1, 0, 0, 0,
c  prflg( 1) = 0, 0, 0,
c  /
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
