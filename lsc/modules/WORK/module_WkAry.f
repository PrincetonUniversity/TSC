c     File: Wkary.inc                            starts
c     scratch, miscelaneous arrays
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER npltdim,nzondim
!============
      INTEGER MxPLTZON
      REAL*8
     ^     wkzr(NZRDIM), ilist(NVELDIM),
     ^     wkv(NWKVDIM), wkpsi(NPSIDIM), wkzx(ZXDIM)
 
      PARAMETER (MxPLTZON=NPLTDIM+NZONDIM)
      REAL*8
     ^     xp(MxPLTZON),  yp(MxPLTZON), xxp(MxPLTZON), yyp(MxPLTZON)
      COMMON / wkcom0 / wkzr, wkzx, wkv, wkpsi, ilist, xp, yp,
     ^     xxp, yyp
c
c     note;  After April 16, 1993 the nwkvdim parameter is active
c            and this line is replaced:
c    ^     wkv(NVELDIM), wkpsi(NPSIDIM), wkzx(ZXDIM)
c     note;  Also, the parameter MxPLTZON introduced for plotting arrays
c     File: Wkary.inc                            ends
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
