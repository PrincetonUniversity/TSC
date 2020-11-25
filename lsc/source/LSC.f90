      SUBROUTINE LSC(LhPwrMWx,                                           &  
     &               nTSCwrix, nTSCreax, nTSCgrax, nLSCcomx, nTSCunux,   &  
     &               nTSCscrx,                                           &  
     &               iRayTrsx, iPlotsix, iErrorx)
!     LSC  main control module              ---------------------------|
!                                                                      |
!                                                                      |
!     LSC  -- Lower hybrid Simulation Code
!     AcDc -- A Current Drive Code
!     Copyright D. W. Ignat, E. J. Valeo,
!               N. J. Fisch, C. F. F. Karney,
!               S. C. Jardin, S. Bernabei, D. P. Enright.
!     Princeton University
!     1992, 1993, 1994, 1995, 1996, 2000
!
!     References:
!       D. W. Ignat, Phys. Fluids, <24>, 1110 (81);
!       E. J. Valeo and D. C. Eder, J. Comp. Physics, <69>, 341 (87).
!       C. F. F. Karney and N. J. Fisch, Phys. Fluids <29>, 180 (86).
!       D. W. Ignat, E. J. Valeo, and S. C. Jardin,
!       "Dynamic Modeling of Lower Hybrid Current Drive,"
!       Nuclear Fusion <34>, 837-852 (94)
!       D. W. Ignat,  B. P. LeBlanc, C. K. Phillips,
!       J. R. Wilson, and R. V. Budny,
!       "Computational Model for Fast Wave Current Drive,
!       Proc. 11th Conf. on Radio Frequency Power in Plasmas,
!       (Palm Springs, CA, May 1995)
!       D. W. Ignat, R. Kaita, S. C. Jardin, and M. Okabayashi
!       "Spreading of Lower Hybrid Wave Driven Currents in PBX-M,"
!       Nucl. Fusion <36> 1733-1742 (96)
!
!            Warning: Local variables are assumed to be SAVEd by the
!            compiler.  While this is a standard implementation, it is
!            not ANSI standard.
!
!     LSC_16Oct00
!        1. Make changes found necessary for running with TSC
!           set nLSCcom2 = nTSCscrn
!           recycle ns30c1 to be nTSCgraf, the unit for gnuplot stuff
!           abandon writes of GrfTitle to ns30c1, or whatever
!        2. Use "call date_and_time(zzdate,zztime,zzzone,ivals)" in gnup.F
!           instead of "call date", "call itime" (or "call clock" for Cray)
!           "date_and_time" is a standard intrinsic in Fortran 90.   But,
!           under Fujitsu Linux an -X9 compiler switch is needed.
!           g77 knows about "date_and_time" .
!        3. Remove spurious "fdate" in XrayCam2.F; add "date_and_time".
!           Remove spurious mention of ezlilo, not used any more.
!     LSC_25Sep00
!        1. Remove unused *.inc : Labels, Diags, local1d, plotf .
!        2. Re-arrange inputs and parameters `graph' for more logical order
!        3. Remove disagreement in VecMnMx over real*4 / real*8
!           (Glitch not previously noticed while ridding code of
!            dependence on the -r8 type compiler switch.)
!        4. Search/fix other similar.  Found in ezsg, agraphdi; Raytrace, GetKscat;
!           jrf, SmoothJrf; Raytrace, ran3
!     LSC_03Aug00
!        0. Pull out Minneapolis/Boulder graphs, testSmooth, SmoothInt
!        1. Correct errors in Jray _graphs_ which eliminated negative J's
!        2. Remove Rayout, mktrajb, mktrajc.  These are ray plots in
!           the vpar-psi plane.  (I did not grasp the difference, surely
!           not obvious on looking at the output.)  These did not use EZ
!           calls, so were getting hard to support.
!        3. Adjust placement of n_\parallel, n_\perp vs rt psi;
!           3d isometric of ray, poloidal cross section of ray;
!           Adjust placement of n_\parallel enhancement and accessibility
!        4. Rename iFastPiece to be iOrigMaxwl in Giruzzi of pbxio.F
!           to make code more understandable.
!        5. Recast and simplify all Plot Print flags; see PlPr.inc
!        6. Add GNplot for rays, Jrf, Prf.  Revise SG-only wave
!           scattering plots to go in the middle of ray plots
!        7. Insert three `dummy' couplers USR1 USR2 USR3 and set them
!           equal to TFTR.  The intent is for these to be modified in the
!           source.
!        8. Fiddle to get gnuplot filename to take date/time from
!           both Fujitsu f90 and G77.
!        9. Many small changes aimed at getting things to work in
!           double precision without "-r8" switches.
!     LSC_27Jun00
!        0. Work on writes to the main log file, improve readability.
!           Small non-important bug in forgetting path length in ray log,
!           but have not found it yet.  (I remember once understanding the issue.)
!     LSC_21Jun00
!        0. Take out unit 6 writes; DoBern, DoComm,  DoEVsm, DoAdFs
!           DoEVsm is still available through iSMO, see below, but recompile.
!        1. Take out iRememFe nRampOld
!     LSC_12Jun00
!        0. Write Pql Pray and integrals for gnuplot
!     LSC_05Jun00
!        0. Remove phony symmetric hi n_parallel spectrum if DoBram = -1
!           (See 30Apr93)
!        1. Write spectrum to data file for gnuplot (SM ?) post processing
!     LSC_19May00
!        0. Take out DoRFHe, DoRFDa, DoPaus, DoScrn, DoDeBg, DoGrzi, DoHarv
!        1. Take out epol Err Ert Err2 Ert2 (NJ Fisch used these once)
!     LSC_17May00
!        0. Change setpwrl to do by default linear ramp up assuming about
!           100 ramp steps, with flat at the top, but keep in the geometric
!           ramp code, both original (2 * per) and trial (1.eps * per).
!           Invent parameter for NFLATDEF; use NRAMPDIM for default nrampup
!        1. Take out graph of Er Ntht Nrho Eth Err
!        2. Fix bug that diffused power can be negative where current is
!           negative.  Just take abs(Jrf) in the weighting formula.
!     LSC_15May00  -- simplify input and remove features not, or seldom, used
!        0. Split NAMELIST files in two pieces,
!           input.lhh containing  `inpvalue'  --- the most  basic  parameters
!           input.xpt containing  `inpexprt'  --- the experts-only parameters
!        1. Take out the stuff put in for Norton L Bretz ( NLB , bret )
!        2. Take out the stuff put in for Nathaniel J Fisch  (Drt DoFish) .
!        3. Take out the stuff put in for Bob Harvey, but save the code
!        4. Take out the stuff put in for Franco Paoletti, but save the code
!        5. Take out non uniform velocity grids
!     LSC_21Apr00  -- interim fixes, not a stable version
!        0.  Redefine rFudgDmp to be the fractional power lost in a zone
!            in the ray picture, and use that fraction constructing Dql
!            Abandon iFudgDmp, since rFudgDmp needed for all zones, rays.
!        1.  Greatly increase the allowed number of iterations on f_e to
!            200, and set up nflat with kludge to select a build up method:
!              nflat = 1,2 ... nramup-1  ==> old method
!              nflat = nrampup           ==> 1.15 increase each cycle
!              nflat > nrampup           ==> 1/nrampup increase each cycle
!        2.  Change the test for stopping the ray calculation so that
!            the calculation runs longer, and therefore can accomodate
!            large amounts of QL burn through.  It used to be that if
!            Maxwellian damping was ``strong'' then the ray stopped, but
!            now we say that Maxwellian damping cannot be more than 10%
!            per region, so it takes hundreds of strong damping zone crossings
!            to absorb the ray.
!     LSC_03Jan00
!        0.  Add '#include "PlPr.inc"' in SUBROUTINE GrafParm to support
!            printing out of "idiag".  See 07Feb96.  g77 found this for me.
!            Gasp.  Thank gnu.
!     LSC_07Feb96
!        0.  Add idiag to plot of paramters.  Missing all this time!
!     LSC_16Jan96
!        0.  Take out references to dpsi.  It was not actually
!        1.  Install the EhstKarney formula.
!     LSC_20Sep95
!        0.  Allow diffusing the rf power with PrfSpred if current has
!            been diffused with DiffuJrf such that
!            Praytot(psi) = Prf-total n(psi) J-diffused(psi) dV(psi) /
!                           [ \Sigma  n(psi) J-diffused(psi) dV(psi) ] *
!                                                 PrfSpred     +
!                  Prf-ray-undiffused(psi)* (1. - PrfSpred)
!        1.  Add phony listing of power to ions in transfer to calling
!            code; this to help with compatibility with TSC/FCD/TRANSP.
!     LSC_14Jul95
!        0.  Warning: Local variables are assumed to be SAVEd by the
!            compiler.  While this is a standard implementation, it is
!            not ANSI standard.  In particular, the NERSC A machine
!            compiler must be forced to do this with an option switch.
!        1.  call EZfini(0,0) consistently; remove 4-transition-point
!            grid stretcher with calls mkgrid eegrd vecadd and vars
!            frv1minus, frv2minus, v1minus, v2minus,
!            frv1plus,  frv2plus,  v1plus,  v2plus,  epsvgr
!
!     LSC_10Feb95
!        0.  Add coupler 'TFTRLHCD' which is 32 0.55cm wgs, 0.15cm septum
!            Move EXCLUDED to 1.2 from 1.5; better for TFTR
!        1.  Repair bug in the error trap on nrays, ntors, npols;
!            add trap of error nparmin>nparmax
!        2.  Replaced "100." with REXPMIN in a couple of  places in FeMaxw;
!            this avoids compiler problem under osf/1 v3.0, f77 v3.6:
!            evaluation of exp( -100. ) gave an arithmetic error.
!            In the TRANSP environment, however, it should be noted that
!            exp( -REXPMIN ) evaluates by silent underflow to zero. --DMC
!        3.  Add time of equilibrium to certain graphs; see EquTitle(81:96)
!     LSC_03Dec94
!        0.  Diffuse Jrf accorting to a method like the Fuchs method.
!            See extensive comments under SmoothJrf.
!            DiffuJrf in the namelist governs this; REPLACES JrfEnhc
!        1.  Add graphs of Jrf by velocity and psi for Bernabei 94 APS meeting
!            see JrfDiagn
!        2.  Move SUBROUTINE LSC to top of file...ahead of comments
!     LSC_30Sep94
!        0.  Insert conditional compilation for TRANSP
!        1.  Split new modules on the /soft/TRANSP-LSC with fsplit
!     LSC_30Aug94
!        0.  Enter option for JET_LHCD coupler
!        1.  Initialize iLastRy in Block Data as Fortran requires
!        2.  Remove infinite loop potential in RayIni....
!        3.  Add nslices to parameter output graph
!        4.  Trap read errors in RdTSC
!        5.  Trap powers(i) = 0.0; reset to 0.001
!     LSC_23May94
!        0.  Install Do0Edc so the Edc in lhcdoua/jardin.d etc can be ignored
!            and treated as zero.  This to investigate possible stability
!            problems with T/LSC combined on some shots like 313258
!        1.  Add graph of n// enhancement on midplane to pitch(r),
!            Ne(r), Pray(r) attached to PlFlg(26)
!     LSC_28Mar94
!        0.  Change the Pitch of the field line graph to degrees; range: +/- 6
!        1.  Add plot of Pray(v) for various psi. Requested by Stefano for the
!            Boulder meeting....to show that PBX is relevant to TPX. Flag 19
!        2.  Also, add plot of Pray(iv) summed over all positions and
!            Plaunch at the edge in v, nparallel, and 1/nparallel(?) space.
!     LSC_10Dec93
!        0.  Allow the use of phony foils to see the Photon Temperature
!            as determined at lower energy might be higher.  Foil code
!            'PF' sets up only Aluminum with the thickness given
!        1.  Allow the introduction of a fast particle population such that
!            T/Tfast = epsT, nfast/n = epsN,  nfast Tfast /(nT) = epsP
!            using variables Tail?eps where ? = N, T, P
!        2.  Apply trapezoid rule for the integration over
!            \nuColl / (Dc + Dql) dv uniformly for +- v and starting and
!            subsequent iterations.
!        3.  Allow for operation on non-symmetric equilibria
!     LSC_18Oct93
!        0.  Allow re-calcluations of rays to take place 1 ray at a time
!            with Do1Rpr=1 (Do 1 Ray per retrace call)  See iLastRy, iDoRay
!        1.  Added entry points to help with TSC restarts:
!            ENTRY PutLSCrestart
!            ENTRY GetLSCrestart
!        2.  Add fast electron population with Doflag DoAdFs.
!            The npsi*nv numbers come from another code GGfe.F
!     LSC_03Sep93
!        0.  Give TSC the power and current !!centered!! on the psi values.
!            Due to a previous misunderstanding this had been given on
!            volume shells ending at the psi values passed.
!        1.  Allow TSC to input directly the new loop voltage and call
!            directly for the new current, and dJ/dE
!              SUBROUTINE JdepFromTSC(curtsc,djdets2,npsiTake) in jrf.F
!              ENTRY PutVinLSC(voltPut,npsiPut)                in grapgrd.F
!              note: location#1 is ignored---first data at second location
!              in LSC  input arrays are 1-->NpsiJ;   in TSC 1 -->npsiPut
!              in LSC output arrays are 1-->NpsiJ-1; in TSC 1 -->npsiTake
!     LSC_31Aug93
!        0.  Add EZdscw EZdscr EZopen EZclos to help capture numbers
!            for later manipulation; employ on vertical xray slice
!        1.  Fix error in paramater-list-graph for nrays
!        2.  work xray units: inten in w/cm^2 on camera; w/cm^3 in plasma;
!            including solid angle effect, see "solAng", and absorbers;
!            make frac = ri - i inst of ri - 1 in interpolator; move Cth0 = mu;
!            fix E-indexing; fix multiplication of normalizations; clamp Emax.
!        3.  Add nRunDot, jRunDot == d/dt {n_e_runaway,j_runaway_from_n_dot}
!            = Dql df/dv * {1; q_e v_parallel} | cooperative runaway velocity
!     LSC_19Jun93
!        0.  Add ability to launch rays in the middle, bypass the ql
!            calculation; for TFTR support as discussed with Norton L. Bretz.
!            DoFlag: DoBret; Starting points rstNLB, zstNLB, npaNLB, degNLB
!            ALL THIS UNDER ``0.'' TAKEN OUT MAY 2000; Sorry, Dr Bretz (NLB)
!        1.  Add ability to write out the psi, r, z, n//, n_rho, etc for
!            Paoletti and Bernabei who are exploring with MatLab the
!            'Volcano Limit' to n//.
!        2.  Also for Paoletti, add r,z,rtpsi,ompe/omce,Bth,Bph
! !      3.  Fixed bugs in Wr2TSC.F and Rayio.F where Npsi was given instead
!            of NpsiJ as size of normalized flux array, for graphs of
!            approximate damping and dlnJdlnE
!     LSC_24May93
!        0.  Repair work on code that writes data for Bob Harvey.
!            sene was wrong; exde etc in doubt; wdnpar erratic
!        1.  Properly put $ at end of pline with call to DolAtEnd(pline)
!        2.  Call FastFrac every time at request of TRANSP
!        3.  Add plot of n// enhancement factor in AccesiSB; puncture plot of
!            actual enhancement
!     LSC_30Apr93
!        0.  Allow for a phony symmetric hi n_parallel spectrum if DoBram = -1
!            This is for an IBW synergy simulation, but not with IBW disp rel.
!        1.  Short circuit warnings if iRayTrsi=0
!     LSC_20Apr93
!        0.  Density fluctuation scattering at bounces at edge; graphs w/ rays
!        1.  Eliminate NparEnhc in namelist; add scatKdeg;
!            add GetKscat and ran3 subroutines; modify BounShft extensively
! !!!    2.  Split the npar array into ntor(NTORDIM) npol(NPOLDIM)
!            These form components of npar(NRAYDIM); add nrays=ntors*npols;
!            ntor replaces npargrid; add npolmax/min;
!        3.  Stopped a ray with more than 4 zones of heavy damping
!            according to a Maxwellian calculation: MxdlnPdz = 0.001
!            This means no penalty for large number of steps.
!        4.  Moved NZONDIM to 2000;
!            removed restriction that NVELDIM should be greater than
!            NPSIDIM + NZONDIM by making WKVDIM = NVELDIM+NPSIDIM+NZONDIM;
!            made plot workspace large enough with MxPLTZON=NPLTDIM+NZONDIM &
!            moving IDM from 200 to NZONDIM.  These requested by Takahashi
!        5.  Changed LSCwarn to be left justified
!     LSC_29Mar93
!        0.  E. Valeo corrected the integral over absorber function,
!            and made other changes to the x ray camera. Constants still to be
!            made correct.  Forced nMUbins to be odd, but not sure this is
!            necessary.  Took out an extra density in the sigma.
!        1.  Small adjustments to absorber code reflecting refinements from
!            S. vonGoeler.  Dimensioning error in GetXmiss was fixed.
!     LSC_15Mar93
!        0.  Improved, added graphs having to do with Xray signal and emission
!        1.  Set ISIZElcfs = 800 because 100 was causing errors
!        2.  Made aspect ratio for coutours of flux and accessibility ok
!     LSC_23Feb93
!        0.  Adjusted Rayio.F on tracing flux surface at 0.3*PsiLim
!     LSC_16Feb93
!        0.  Allowed edge density to be reduced with DoBern=1.  See grapgrd.F
!        1.  Made namelist Inpval an include file at request of TRANSP
!     LSC_06Feb93
!        0.  Improved labeling of launched spectrum graph
!        1.  Trapped and reported possible errors on coupler names
!        2.  Added a very slow coupler, called SLOWSLOW with 16 wgs
!            This not accessible from TRANSP namelist which is integer-driven
! !!!    3.  Corrected error that fghz was not being reported to brambJES !
!     LSC_03Feb93
!        0.  Added plot of power in ray vs root psi, with the first 4 passes
!            at no damping followed in the graph
!        1.  Adjusted step in psi surface plotter to handle D-III-D extremes.
!            Also adjusted start search point to be 0.9 of Rmag
!        2.  Took out DoRFDa and DoRFHe as functioning; left in namelist
!            but took out of include file Doflags.inc
!     LSC_28Jan
!        0.  Made SMALL in fe.F 2.5e-05 for TRANSP use.
!     LSC_22Jan93
!        0.  Allowed Ip as well as ne Te Bo to be changed on DoBern Switch
! !!!    1.  Made adjustment to Dql in zones of high damping so Pql Pray
!            should be pretty much equal.  Key variables: iFudgDmp rFudgDmp
!        2.  Added TORSUPRA coupler option; phaseDeg follows CEN definition
!     LSC_15Jan93
!        0.  Made normalization of parametrized spectrum (DoBram=0) correct.
!        1.  Allowed for Bo to be changed on the DoBern switch
!        2.  Attempted to fix Pray/Pql; by using only discrete n_parallel
!     LSC_08Jan93 :
!        0.  Tokamak de Varennes and PBXM slow couplers selected by character
!            variable couplers(nGrups)
!        1.  dtdV substituted for dVdt in Raystore cause dt can be 0 in TRANSP
!        2.  Ray data written for RW Harvey & CQL3D  with switch DoHarv
!        3.  Grid for reporting Prf and Jrf shifted back to be consistent
!            with comments.  Thanks to Ted.
!        4.  Camera graphs buffer flushing worked on.  Seems ok now.
!        5.  Parameters not graphed if DoTRAN .eq. 1 as requested by TRANSP
!     LSC_13Nov92 :
!        0.  copy of 7Nov TSC code for J/E iteration in jrf.F as comments
!        1.  fix counting in Giruzzi graphs
!     LSC_07Nov92 :
!        0.  SMALL set to 1.5e-5 in fe.F for TRANSP
!        1.  x ray camera w absorbers in XrayCam2 by Valeo, Enright, Ignat
!        2.  last closed flux surface R_, Z_" lcfs " min/max known for graphs;
!        3.  dJ/dE passed to TSC/TRANSP after E/J dJ/dE ;
!        4.  Effective date written with parameters and inputs
!
!      SUBROUTINE LSC(LhPwrMWx,
!     ^               nTSCwrix, nTSCreax, nTSCgrax, nLSCcomx, nTSCunux,
!     ^               nTSCscrx,
!     ^               iRayTrsx, iPlotsix, iErrorx)
!                                       IMPLICIT NONE or equivalent
!                                       IMPLICIT UNDEFINED (A-Z)
!                                       is in file Implic.inc
      USE Doflags
      USE Escan
      USE Implic
      USE params
      USE PlPr
      USE RayBins
      USE RayWrk
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit REAL declarations:
!============
      CHARACTER*8 LSCdate
      COMMON /LSCversn/ LSCdate
      CHARACTER*70 ErrMsg
      INTEGER     nTSCwrix, nTSCreax, nTSCscrx, nTSCgrax, nLSCcomx,      &  
     &            nTSCunux, iRayTrsx, iPlotsix, iErrorx
      REAL*8        LhPwrMWx
      REAL*8    PwrFnd, JrfFnd, PqlFnd
      INTEGER ifirst, iDoRay
!     LhPwrMWx       Lower Hybrid power in MW passed from TSC
!     nTSCwrit       Unit number to which TSC writes equilibrium data
!     nTSCread       Unit number from which TSC reads Prf, Jrf , etc
!     nTSCscrn       Unit number TSC uses for writing to the screen
!     nTSCgraf       Unit number to which LSC writes gnuplot stuff
!                      for each graph
!     nLSCcomm       Unit number to which LSC sends 'screen' output
!     nTSCunused     Unit number not used by TSC; used internally by LSC
!     iRayTrsi     0 use old ray data, and old f_e(v);
!                    use new E_dc for the current
!	           1 calculate new rays, and f_e(v), from new equilibrium
!	           2 use old ray data, but calculate new f_e(v)
!                    taking account of new n_e and T_e,
!                    use new E_dc for the current
!
!     iPlotsig     0 do not make plots; 1 make plot files asked in input.lhh
!     iError       0 LSC finishes without errors; 1 (or more) errors found
!                 -1 LSC found an error; calling program can keep going
!
!     Confusion posible over the suffix Vec vec and Ary ary; explained below:
!     Vec as in PsiVec ... is on the TSC grid, not on the LSC grid: that is Ary
!     Ary as in PsiAry ... in on the LSC grid, not on the TSC grid: that is Vec
!
!     -                                 iLastRy is the index of the last ray
!     -                                 traced; or -1 if no rays have been
!     -                                 traced.  iLastRy controls iDoRay
!                                       iDoRay is a flag in DoRay
!                                       if 1 thru nrays then trace that ray
!                                       if otherwise, then trace all rays
!
 
      DATA ifirst /1/
!     -                                 Copy all the variables and flags
!     -                                 given by TSC into a common block
!     -                                 except iError which is always set to 0.
!     -                                 (value passed from TSC ignored)
!     LSCdate  = 'ddmmmyy '             LSC_date
      LSCdate  = '16Oct00 '
      TotPwr   =  LhPwrMWx * 1.0E+06_R8
      nTSCwrit =  nTSCwrix
      nTSCread =  nTSCreax
      nTSCscrn =  nTSCscrx
      nTSCgraf =  nTSCgrax
      nLSCcomm =  nLSCcomx
      nLSCcom2 =  nTSCscrx
      nTSCunus =  nTSCunux
      iRayTrsi =  iRayTrsx
      iPlotsig =  iPlotsix
      iError   =  0
      iErrorx  =  0
 
!     -                                 If this is the first call,
!     -                                 force things a certain way, such as
!     -                                 make sure the rays are traced, and set
!     -                                 the TSC flag to 1 accordingly
!
!     -                                 Clear warning counter; iEdc flag
      call LSCclear
      call RdVars
           if (iError .gt. 0) then
              call LSCtrace ( ' RdVars')
              iErrorx = iError
              return
           endif
!
      if (Do0Edc .ne. 0) then
         iEdc   = 1
         EdcInp = 0.00_R8
!     .                                 set the iEdc switch which was installed
!     .                                 long ago for debugging the sign of
!     .                                 current driven etc; and set the value
!     .                                 of the EdcInp to 0.00
      endif
!
      call RdTSC
           if (iError .gt. 0) then
              call LSCtrace ( ' RdTSC')
              iErrorx = iError
              return
           endif
!
      if(iRayTrsi .eq. 0 .and. ifirst .eq. 1) then
        iRayTrsi =  2
      endif
!     .                                 if iRayTrsi=0 and firstcall, then
!     .                                 make iRayTrsi=2 so we find f_e etc
 
      if(iPlotsig .eq. 0) call NoPlots
!
 
      if(iRayTrsi .ne. 0 .or. ifirst .eq. 1)  then
        ifirst = 0
        call ValInit
!                                       initialization routines
!                                       for quasilinear calculations
!                                       -- needed after
!                                       each TSC read
        if (iError .gt. 0) then
          call LSCtrace ( ' ValInit')
          iErrorx = iError
          return
        endif
      endif
!
        if(iRayTrsi .eq. 1 .or. iRayTrsi .eq. 2)  then
          call SpecInit
!     .                                 generate launch spectrum with
!     .                                 nGrps humps at centers(nGrps) with
!     .                                 widths(nGrps) & relative powers(nGrps).
!     .                                 But if nGrps = 0, then get a Brambilla
!     .                                 spectrum using the JE Stevens code.
             if (iError .gt. 0) then
                call LSCtrace ( ' SpecInit')
                iErrorx = iError
                return
             endif
        endif
!
!
!
      if (iRayTrsi .ne. 0) then
        call PowrInit
!     .                                 apply the launch spectrum to the
!     .                                 total power.
           if (iError .gt. 0) then
              call LSCtrace ( ' PowrInit')
              iErrorx = iError
              return
           endif
!
!       call RyZnInit     !! 9sep 93
!                                       initialization routines
!                                       for ray calculations
!          if (iError .gt. 0) then       !! 9sep 93
!             call LSCtrace (' RyZnInit') !! 9sep 93
!             iErrorx = iError !! 9sep 93
!             return !! 9sep 93
!          endif !! 9sep 93
      endif
!
      if(iRayTrsi .eq. 1 .and.                                           &  
     &     Do1Rpr .eq. 0)         then
              iDoRay = nrays+1
         call DoRay (iDoRay)
           if (iError .gt. 0) then
              call LSCtrace (' DoRay')
              iErrorx = iError
              return
           endif
         call WrRay
           if (iError .gt. 0) then
              call LSCtrace (' WrRay')
              iErrorx = iError
              return
           endif
      endif
!
      if(iRayTrsi .ne. 1 .and.                                           &  
     &     Do1Rpr .eq. 0)         then
         call RdRay
           if (iError .gt. 0) then
              call LSCtrace (' RdRay')
              iErrorx = iError
              return
           endif
      endif
!
      if(iRayTrsi .eq. 1 .and.                                           &  
     &     Do1Rpr .ne. 0)         then
         iDoRay = iLastRy + 1
         if (iDoRay .gt. nrays) iDoRay=1
         call DoRay (iDoRay)
         iLastRy = iDoRay
           if (iError .gt. 0) then
              call LSCtrace (' DoRay')
              iErrorx = iError
              return
           endif
           if(iLastRy .eq. nrays .or.                                    &  
     &        iLastRy .eq.     0      ) then
             call WrRay
             if (iError .gt. 0) then
               call LSCtrace (' WrRay')
               iErrorx = iError
                return
             endif
           endif
      endif
 
!      if(iRayTrsi .ne. 1 .and.          comment out because of other
!     ^     Do1Rpr .ne. 0)         then  initializations
!         call RdRay
!           if (iError .gt. 0) then
!              call LSCtrace (' RdRay')
!              iErrorx = iError
!              return
!           endif
!      endif
 
!                                       trace/read ray trajectories,
!                                       polarization
!     call ryzcent
!                                       compute zone center ray indices
!                                       as required for power deposition
      if (iRayTrsi .ne. 0) then
 
        call Ramp_Pwr
           if (iError .gt. 0) then
              call LSCtrace (' RampPwr')
              iErrorx = iError
              return
           endif
!                                       ramp up power
!
        call PdepCalc
!                                       compute Prf
      endif
 
      call JdepCalc
          call JrfDiagn
!         BERNABEI request for 1994 APS Meeting in Minneapolis
!                                       compute (and plot) jrf
           if (iError .gt. 0) then
              call LSCtrace(' JdepCalc')
              iErrorx = iError
              return
           endif
      call JsplCalc
!                                       compute Js at E_cd pl(us) dE_dc
           if (iError .gt. 0) then
              call LSCtrace(' JsplCalc')
              iErrorx = iError
              return
           endif
 
      call SmoJandP(Rmax-Rmaj)
!
      call Wr2TSC(PwrFnd, JrfFnd, PqlFnd)
!                                       write deposition profile
!                                       and current drive profile
!                                       and dJ/dE profile to TSC
!
           if (abs(TotPwr-PwrFnd) .ge. 0.2_R8*TotPwr .or.                &  
     &         abs(PqlFnd-PwrFnd) .ge. 0.2_R8*TotPwr )         then
              write(ErrMsg,'('' TotPwr PwrFnd JrfFnd PqlFnd: '',         &  
     &                            4(1pe10.2) )')                         &  
     &                          TotPwr,PwrFnd,JrfFnd,PqlFnd
              call LSCwarn(ErrMsg)
           endif
!     .
!     .                                 write a warning if power is not
!     .                                 well absorbed or if ray/ql answers
!     .                                 are not too close; give current too
!
           if (iError .gt. 0) then
              call LSCtrace(' Wr2TSC')
              iErrorx = iError
              return
           endif
 
!      call ProfOut ! no more writing of profiles July 2000
 
      if(PlFlg(DAMPL) .ge. TRUE) then
         call FirstPas
      endif
 
      if(PlFlg(DQLFEPL) .ge. TRUE) then
         call Giruzzi
      endif
!
!                                      Inserted, August '92 by Doug P. Enright
!                                      for - Xray bremsstrahlung
!                                      routine
!
      if (DoXcam .eq. 1 .and. iPlotsig .eq. 1) then
         call xRayCam2
      endif
 
           if (iEndRy .gt. 0) then
              call LSCwarn ( ' ray stopped early; see neg rc')
              iErrorx =-iEndRy
              return
           endif
 
      return
!
      ENTRY PutLSCrestart
      call WrRay
      return
!
      ENTRY GetLSCrestart
      iError = 0
      call RdRay
      if (iError .gt. 0 ) then
        call LSCwarn(' no LSC restart file; use dead start')
        return
      endif
      ifirst = 0
      return
 
      END
!                                                                      |
!                                                                      |
!     AcDc main control module ends         ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
