c     tscunits.inc      ------------------------------------------------
c
c     nTSCwrit          Unit number to which TSC writes equilibrium data
c     nTSCread          Unit number from which TSC reads Prf, Jrf , etc
c     nTSCscrn          Unit number TSC uses for writing to the screen
c     nTSCgraf          Unit number to which LSC writes gnuplot stuff
c                         for each graph
c     nLSCcomm          Unit number to which LSC sends 'communication' output
c     nTSCunus          Unit number not used by TSC; used internally by LSC
c  ORIGINAL USAGE of iRayTrsi:
c     iRayTrsi     0 use old ray data; recalculate f_e(v)
c                    taking account of new n_e and T_e,
c                    use new E_dc for the current
c	           1 calculate new rays, and f_e(v), from new equilibrium
c  NEW USAGE of iRayTrsi:
c     iRayTrsi     0 use old ray data, and old f_e(v);
c                    use new E_dc for the current
c	           1 calculate new rays, and f_e(v), from new equilibrium
c	           2 use old ray data, but calculate new f_e(v)
c                    taking account of new n_e and T_e,
c                    use new E_dc for the current
c
c     iPlotsig          0 do not make plots; 1 make plots as in input.lhh
c     iError            0 LSC finishes without errors; 1 (or more) errors found
c                      -1 LSC found an error; calling program can keep going
c     iEndRy            Number of aborted rays on this call; failed restart of
c                       error encountered in path, such as too short wavelength
c     nLSCcom2          Alternate for LSC communication
      INTEGER         nTSCwrit, nTSCread, nTSCscrn,
     ^                nTSCgraf, nLSCcomm, nTSCunus,
     ^                iRayTrsi, iPlotsig, iXraysi, iError, iEndRy
      INTEGER         nLSCcom2
      COMMON/tscunits/nTSCwrit, nTSCread, nTSCscrn,
     ^                nTSCgraf, nLSCcomm, nTSCunus,
     ^                iRayTrsi, iPlotsig, iXraysi, iError, iEndRy,
     ^                nLSCcom2
c
c     tscunits.inc      ------------------------------------------------
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
