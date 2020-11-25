!#include "f77_dcomplx.h"
 
!     ------------------------------------------------------------------
!     Authors:
!     E. J. Valeo and  D. P. Enright, August 1992
!     Modified for loading with LSC by D. W. Ignat, October 1992
!     Modified August 24 1993 by Ignat to get w/cm^2; also fix some
!     indexing errors.
!     Copyright:
!     Princeton University, Plasma Physics Laboratory, 1992, 1993, by
!     E. J. Valeo, D. W. Ignat, S. von Goeler, S. C. Jardin,
!     P. G. Roney, J. E. Stevens, and D. P. Enright
!
!     Units: Energy in MeV mostly; kev in absorber routines;
!     mc^2 in the evaluation of Koch and Motz formula.
!     Bremsstrahlung photon radiation in w/cm^2 on camera; w/cm^3 in plasma
!
!     References:
!     H. W. Koch and J. W. Motz,
!     ``Bremsstrahlung Cross-Section Formulas and Related Data,''
!     \RMP{31}{920}{59}.
!     Formula 2BN on p 924
!     '2' denotes differential in photon energy and angle
!     'B' denotes Born approximation
!     'N' denotes no screening
!     and
!     X ray crosss sections based on W. H. McMaster, el al, UCRL 50174
!     Sec. II rev. 1,  May 1969.
!
!     .         --------------------------------------------------------
!
!#define float dble
!     -----------------------------------------------------------------|
!     -                                                                |
!     -                                                                |
      SUBROUTINE BrambJES (fghz, phaseDeg, nrays, nparalls, amplitds,    &  
     &                     nparaMin, nslices, coupler, turnnegs)
!     J. E. Stevens' Brambilla coupling code
!
!     Copyright 10/22/80     by      J. E. Stevens
!	edit:19-feb-87 [Phyllis Roney] convert graphics to sg for vax
!       (and other fortran 77 modifications).
!
!       Compute spectrum of lower hybrid antenna,
!       power flow, and reflection in guides
!
!     Copyright November 1991 by
!     J. Stevens, S. Bernabei, D. Ignat, E. Valeo, S. Jardin
!     with modifications by Ignat to be called by LSC/TSC.
!     -
!     arguments developed for LSC/TSC
!           phaseDeg    phase in degrees between waveguides               INPUT
!           nrays       number of nparallels to return to make spectrum   INPUT
!           nparalls    array of n_par values from low to high           OUTPUT
!           amplitds    array of amplitudes corresponding to those n_par OUTPUT
!           nparaMin    minimum value of abs(n_par) of interest, say 1.5  INPUT
!           nslices     number of slices taken on spectrum in -8 to 8     INPUT
!           turnnegs    negative going spectral elements are turned to    INPUT
!                       positive if this is 1 and if phaseDeg is different
!                       from 180 or 0 by more than 20 degrees
!
!
!     input variables for original code (JES)
!
!           nguide      no. of waveguides
!                       .lt. 0 => phase of total b field fixed
!           nd          total no. of integration steps
!           delphi      phase diff. in degrees to step
!                       .lt. 14  => use phase in input for guide
!           f0          rf frequency in MHz
!           epsr        relative dielectric constant of material in guides
!                       .lt. 0 => don't include guide impedance in calc.
!           x0          antenna to plasma spacing in cm
!           gradn       density gradient at plasma edge in cm-4
!                       .lt. 0 => compute for +4 orders of magnitude
!           a           waveguide height in cm
!           anmax       max. nz to use in integration
!           anplt       max. nz to use in plot
!           xprobe      rf probe position(cm) at which
!                       total e & b fields are calculated
!           edged       plasma density at the edge in cm-3
!                       .lt. 0 => compute +3 orders of magnatude
!           p0(i)       incident power in i'th guide in kw
!                       .le. 0. => passive guide
!           phi(i)      phase of electric field in i'th guide in degrees
!           az(i)       z distance to  one  edge of i'th guide in cm
!           bz(i)       z distance to other edge of i'th guide in cm
!
!     if not enough lines of az,bz,phi,p0 are specified, then
!     the program assumes az(n,1) is the septum width, bz(1)
!     is the guide to guide spacing, phi(1) is the relative
!     phase, and p0(1) is the power for the remaining guides.
!
!     The NumRay largest components are selected for later
!     ray tracing.
 
! The JET_LHCD launcher consists of 6 rows, by 8 columns of multijunctions.
! Each multijunction, consists at the interface with the plasma of 2 rows by 4
! collumbs of waveguides. The 4 waveguides in a row are phased, by fixed phase
! shifters in the waveguide, giving a 90 deg phase difference between the
! forward travelling wave in adjacent waveguides.
!
! In total the plasma launcher interface consists of 12 rows with
! 4X8 = 32 active waveguides in each. On either side of each row of waveguides,
! are 2 passive waveguides.
!
! The dimmension of the wavguides are 72.1mm high by 9mm wide. The
! walls between the waveguides are 2mm wide. The depths of the passive
! waveguides is 3.6 cm for the two on the left of the row, and 5.9cm for the
! two on the right, as seen from the plasma towards the launcher. The fixed
! phasing of a multijunction is such that the phases as seen from the plasma
! towards the launcher is from left to right 0 -90 -180 90(-270) .
!
! The vertical position of the middle of the rows of waveguides, with
! respect to the midplane of the machine (z=0) i as follows:
!   row 1   407mm  row 2  333mm  row 3  259mm  row 4  185mm
!   row 5   111mm  row 6   37mm  row 7  -37mm  row 8  -111mm
!   row 9  -185mm  row10 -259mm  row11 -333mm  row12  -407mm
! The 90 deg phase difference between the wave travelling towards the
! plasma in adjacent waveguides only holds when the match to the plasma is
! good. Otherwise the properties of the multijunction has to be taken into
! account. Assuming 90 degrees is a good approximation in most relevant
! situations though, as we normally need to operate with a good match.
!
! The laucher is shaped at the plasma interface to fit the plasma, this
! means that the waveguides are cut at an agle different to 90 deg at the
! plasma interface. What influence this has on the launched spectra I don't
! know - we havent computed it, but I dont think it alters the spectra much as
! long as we have a good match.
!
! When we talk about a phasing of the launcher we refer to the phase
! difference between the leftmost waveguides in adjacent multijunctions. Most
! of our experiments have been carried out at 0 deg phasing which gives a
! n_\parallel peak at 1.84.  --- Morten Lennholm;  ml@jet.uk
 
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     Argument list first:
      CHARACTER*8  coupler
      INTEGER nrays, nslices, turnnegs
      REAL*8                                                             &  
     &        fghz,                                                      &  
     &        phaseDeg, nparalls(nrays), amplitds(nrays),                &  
     &                  nparaMin
!     Argument list ends.
!
!     These following are pretty much the original variables.
      INTEGER    DIM,      NGPBX,      NPDIM
      PARAMETER (DIM = 40, NGPBX = 32, NPDIM = 601)
      INTEGER    NGTdeV,   NGTors, NGSslow, NGJET, NGTFTR
      PARAMETER (NGTdeV = 32, NGTors = 34, NGSslow = 16)
      PARAMETER (NGJET  = 32, NGTFTR = 32)
 
!     Three "user" couplers are put here as place holders.
!     The original entry for all was identical to TFTR.
!     The intent is that these entries can be modified
!     by the user.
!     Key variables are NGUSR1,   NGUSR2,   NGUSR3
!     Key names are       USR1SPEC, USR2SPEC, USR3SPEC
!
      INTEGER    NGUSR1,      NGUSR2,      NGUSR3
      PARAMETER (NGUSR1 = 32, NGUSR2 = 32, NGUSR3 =32)
!     CHARACTER*20 filename
      CHARACTER*75 pline
      COMPLEX*16                                                         &  
     &        e0(DIM),      r(DIM),    bt(DIM), et(DIM), bb(DIM),        &  
     &        aa(DIM,DIM), ai(DIM,DIM),                                  &  
     &        c0, c1, c2, c3,                                            &  
     &        c5, c6, c7, c8, c9,                                        &  
     &        xp, rho, enz, anx, anzi, denom
      INTEGER nplot, nguide, nd, nphase, nedgd, ngradn, nph, ngn
      INTEGER nloop, ned, nl, nmod, ngset, n, n2
      INTEGER idim, isw, iepsw, ifast, i, iskip, istop, idel, ix1, ix2
      INTEGER l, k, ll, ix, iy, ngpr1, ngpr2, npmin, npmax
      INTEGER lsmx, ndif
      INTEGER NumRay, NumBeg
      INTEGER ii
      INTEGER G
      REAL*8                                                             &  
     &        a1, a2, PI, delt, delphi, f0, epsr, x0, gradn, gradn0,     &  
     &        a, anmax, anplt, anz, anzac, anx3, air, bir, aip, bip,     &  
     &        al0, w0,                                                   &  
     &        xprobe, edged, pmax, xmin, xmax,                           &  
     &        edged0, dencr, pave, zmin, zmax,                           &  
     &        twidth, wdc, w, beta, beta1, ratio,                        &  
     &        alp, alph3, ph0, ph1, ph2, del, c10
      REAL*8                                                             &  
     &        dnz, pplas, pint, area, rr, ri, etr, eti, btr, bti,        &  
     &        arg, teff, pdave, epsrp, sign,                             &  
     &        sum, smx, t1, danplt, asum1, asum2, asum3, fracc, t2, t
      REAL*8                                                             &  
     &        az(10,DIM),  bz(10,DIM),  b(DIM),                          &  
     &        p0(DIM),  phi(DIM), phi0(DIM),  ei(DIM),                   &  
     &        rm(DIM),  rp(DIM),  p0r(DIM),                              &  
     &        etph(DIM),sp(NPDIM), s(NPDIM),    sn(NPDIM),               &  
     &        etm(DIM), btm(DIM), btph(DIM),  phold(DIM)
      REAL*8                                                             &  
     &        RelPha, MinNpar, TotPwr, MaxPwr
      DATA       a1,        a2,        PI, idim, nplot /                 &  
     &     0.355028_R8, 0.2588194_R8, 3.1415926_R8,   40,   301 /
      DATA     nd, delphi,    f0,  epsr,  x0,   gradn  /                 &  
     &       4000,    00._R8, 4600._R8,    1._R8,  0._R8,  1.0E12_R8/
      DATA      a, anmax, anplt, xprobe, edged                /          &  
     &        5.8_R8,   30._R8,   10._R8,     0._R8, 1.5E12_R8/
 
      DATA ( az(1,i), i = 1, NGPBX  )                               /    &  
     &     00.00_R8, 00.71_R8,  1.42_R8,  2.13_R8,  3.15_R8,  3.86_R8,   &  
     &     4.57_R8,  5.28_R8,                                            &  
     &      6.30_R8,  7.01_R8,  7.72_R8,  8.43_R8,  9.45_R8, 10.16_R8,   &  
     &    10.87_R8, 11.58_R8,                                            &  
     &     12.60_R8, 13.31_R8, 14.02_R8, 14.73_R8, 15.57_R8, 16.46_R8,   &  
     &    17.17_R8, 17.88_R8,                                            &  
     &     18.90_R8, 19.61_R8, 20.32_R8, 21.03_R8, 22.05_R8, 22.76_R8,   &  
     &    23.47_R8, 24.18_R8/
      DATA ( bz(1,i), i = 1, NGPBX  )                               /    &  
     &     00.50_R8,  1.21_R8,  1.92_R8,  2.63_R8,  3.65_R8,  4.36_R8,   &  
     &     5.07_R8,  5.78_R8,                                            &  
     &      6.80_R8,  7.51_R8,  8.22_R8,  8.93_R8,  9.95_R8, 10.66_R8,   &  
     &    11.37_R8, 12.08_R8,                                            &  
     &     13.10_R8, 13.81_R8, 14.52_R8, 15.23_R8, 16.25_R8, 19.96_R8,   &  
     &    17.67_R8, 18.38_R8,                                            &  
     &     19.40_R8, 20.11_R8, 20.82_R8, 21.53_R8, 22.55_R8, 23.26_R8,   &  
     &    23.97_R8, 24.68_R8/
      DATA ( az(2,i), i = 1, NGPBX  )                               /    &  
     &      0.00_R8,  0.61_R8,  1.22_R8,  1.83_R8,  2.44_R8,  3.05_R8,   &  
     &     3.66_R8,  4.27_R8,                                            &  
     &      4.88_R8,  5.49_R8,  6.10_R8,  6.71_R8,  7.32_R8,  7.92_R8,   &  
     &     8.53_R8,  9.14_R8,                                            &  
     &      9.75_R8, 10.36_R8, 10.97_R8, 11.58_R8, 12.19_R8, 12.80_R8,   &  
     &    13.41_R8, 14.02_R8,                                            &  
     &     14.63_R8, 15.24_R8, 15.85_R8, 16.46_R8, 17.07_R8, 17.68_R8,   &  
     &    18.29_R8, 18.90_R8/
      DATA ( bz(2,i), i = 1, NGPBX  )                               /    &  
     &      0.51_R8,  1.12_R8,  1.73_R8,  2.34_R8,  2.95_R8,  3.56_R8,   &  
     &     4.17_R8,  4.78_R8,                                            &  
     &      5.38_R8,  5.99_R8,  6.60_R8,  7.21_R8,  7.82_R8,  8.43_R8,   &  
     &     9.04_R8,  9.65_R8,                                            &  
     &     10.26_R8, 10.87_R8, 11.48_R8, 12.09_R8, 12.70_R8, 13.31_R8,   &  
     &    13.92_R8, 14.53_R8,                                            &  
     &     15.14_R8, 15.75_R8, 16.36_R8, 16.97_R8, 17.58_R8, 18.19_R8,   &  
     &    18.80_R8, 19.41_R8/
 
      DATA ( az(3,i), i = 1, NGTdeV  )                               /   &  
     &     00.00_R8, 00.70_R8,  1.40_R8,  2.10_R8,  2.95_R8,  3.65_R8,   &  
     &     4.35_R8,  5.05_R8,                                            &  
     &      5.90_R8,  6.60_R8,  7.30_R8,  8.00_R8,  8.85_R8,  9.55_R8,   &  
     &    10.25_R8, 10.95_R8,                                            &  
     &     11.80_R8, 12.50_R8, 13.20_R8, 13.90_R8, 14.75_R8, 15.45_R8,   &  
     &    16.15_R8, 16.85_R8,                                            &  
     &     17.70_R8, 18.40_R8, 19.10_R8, 19.80_R8, 20.65_R8, 21.35_R8,   &  
     &    22.05_R8, 22.75_R8/
      DATA ( bz(3,i), i = 1, NGTdeV  )                               /   &  
     &     00.55_R8,  1.25_R8,  1.95_R8,  2.65_R8,  3.50_R8,  4.20_R8,   &  
     &     4.90_R8,  5.60_R8,                                            &  
     &      6.45_R8,  7.15_R8,  7.85_R8,  8.55_R8,  9.40_R8, 10.10_R8,   &  
     &    10.80_R8, 11.50_R8,                                            &  
     &     12.35_R8, 13.05_R8, 13.75_R8, 14.45_R8, 15.30_R8, 16.00_R8,   &  
     &    16.70_R8, 17.40_R8,                                            &  
     &     18.25_R8, 18.95_R8, 19.65_R8, 20.35_R8, 21.20_R8, 21.90_R8,   &  
     &    22.60_R8, 23.30_R8/
!     DATA ( az(4,i), i = 1, NGTors  )                               /
!    ^      0.00,  0.95,  1.90,  2.85,  3.80,  4.85,  5.80,  6.75,
!    ^      7.70,  8.75,  9.70, 10.65, 11.60, 12.65, 13.60, 14.55,
!    ^     15.50, 16.55, 17.50, 18.45, 19.40, 20.45, 21.40, 22.35,
!    ^     23.30, 24.35, 25.30, 26.25, 27.20, 28.25, 29.20, 30.15,
!    ^     31.10, 32.05                                           /
!     DATA ( bz(4,i), i = 1, NGTors  )                               /
!    ^      0.85,  1.80,  2.75,  3.70,  4.65,  5.70,  6.65,  7.60,
!    ^      8.55,  9.60, 10.55, 11.50, 12.45, 13.50, 14.45, 15.40,
!    ^     16.35, 17.40, 18.35, 19.30, 20.25, 21.30, 22.25, 23.20,
!    ^     24.15, 25.20, 26.15, 27.10, 28.05, 29.10, 30.05, 31.00,
!     ^     31.95, 32.90                                           /
      DATA ( az(4,i), i = 1, NGTors  )                               /   &  
     &      0.00_R8,  1.05_R8,  2.10_R8,  3.15_R8,  4.20_R8,  5.60_R8,   &  
     &     6.65_R8,  7.70_R8,                                            &  
     &      8.75_R8, 10.15_R8, 11.20_R8, 12.25_R8, 13.30_R8, 14.70_R8,   &  
     &    15.75_R8, 16.80_R8,                                            &  
     &     17.85_R8, 19.25_R8, 20.30_R8, 21.35_R8, 22.40_R8, 23.80_R8,   &  
     &    24.85_R8, 25.90_R8,                                            &  
     &     26.95_R8, 28.35_R8, 29.40_R8, 30.45_R8, 31.50_R8, 32.90_R8,   &  
     &    33.95_R8, 35.00_R8,                                            &  
     &     36.05_R8, 37.10_R8/
 
      DATA ( bz(4,i), i = 1, NGTors  )                               /   &  
     &      0.85_R8,  1.90_R8,  2.95_R8,  4.00_R8,  5.05_R8,  6.45_R8,   &  
     &     7.50_R8,  8.55_R8,                                            &  
     &      9.60_R8, 11.00_R8, 12.05_R8, 13.10_R8, 14.15_R8, 15.55_R8,   &  
     &    16.60_R8, 17.65_R8,                                            &  
     &     18.70_R8, 20.10_R8, 21.15_R8, 22.20_R8, 23.25_R8, 24.65_R8,   &  
     &    25.70_R8, 26.75_R8,                                            &  
     &     27.80_R8, 29.20_R8, 30.25_R8, 31.30_R8, 32.35_R8, 33.75_R8,   &  
     &    34.80_R8, 35.85_R8,                                            &  
     &     36.90_R8, 37.95_R8/
 
      DATA ( az(5,i), i = 1, NGSslow  )                               /  &  
     &      0.00_R8,  0.30_R8,  0.60_R8,  0.90_R8,  1.20_R8,  1.50_R8,   &  
     &     1.80_R8,  2.10_R8,                                            &  
     &      2.40_R8,  2.70_R8,  3.00_R8,  3.30_R8,  3.60_R8,  3.90_R8,   &  
     &     4.20_R8,  4.50_R8/
      DATA ( bz(5,i), i = 1, NGSslow  )                               /  &  
     &      0.25_R8,  0.55_R8,  0.85_R8,  1.15_R8,  1.45_R8,  1.75_R8,   &  
     &     2.05_R8,  2.35_R8,                                            &  
     &      2.65_R8,  2.95_R8,  3.25_R8,  3.55_R8,  3.85_R8,  4.15_R8,   &  
     &     4.45_R8,  4.75_R8/
 
      DATA ( az(6,i), i = 1, NGJET  )                               /    &  
     &      0.00_R8,  1.10_R8,  2.20_R8,  3.30_R8,  4.40_R8,  5.50_R8,   &  
     &     6.60_R8,  7.70_R8,                                            &  
     &      8.80_R8,  9.90_R8, 11.00_R8, 12.10_R8, 13.20_R8, 14.30_R8,   &  
     &    15.40_R8, 16.50_R8,                                            &  
     &     17.60_R8, 18.70_R8, 19.80_R8, 20.90_R8, 22.00_R8, 23.10_R8,   &  
     &    24.20_R8, 25.30_R8,                                            &  
     &     26.40_R8, 27.50_R8, 28.60_R8, 29.70_R8, 30.80_R8, 31.90_R8,   &  
     &    33.00_R8, 34.10_R8/
 
      DATA ( bz(6,i), i = 1, NGJET  )                               /    &  
     &      0.90_R8,  2.00_R8,  3.10_R8,  4.20_R8,  5.30_R8,  6.40_R8,   &  
     &     7.50_R8,  8.60_R8,                                            &  
     &      9.70_R8, 10.80_R8, 11.90_R8, 13.00_R8, 14.10_R8, 15.20_R8,   &  
     &    16.30_R8, 17.40_R8,                                            &  
     &     18.50_R8, 19.60_R8, 20.70_R8, 21.80_R8, 22.90_R8, 24.00_R8,   &  
     &    25.10_R8, 26.20_R8,                                            &  
     &     27.30_R8, 28.40_R8, 29.50_R8, 30.60_R8, 31.70_R8, 32.80_R8,   &  
     &    33.10_R8, 34.20_R8/
 
      DATA ( az(7,i), i = 1, NGTFTR )                               /    &  
     &      0.00_R8,  0.70_R8,  1.40_R8,  2.10_R8,  2.80_R8,  3.50_R8,   &  
     &     4.20_R8,  4.90_R8,                                            &  
     &      5.60_R8,  6.30_R8,  7.00_R8,  7.70_R8,  8.40_R8,  9.10_R8,   &  
     &     9.80_R8, 10.50_R8,                                            &  
     &     11.20_R8, 11.90_R8, 12.60_R8, 13.30_R8, 14.00_R8, 14.70_R8,   &  
     &    15.40_R8, 16.10_R8,                                            &  
     &     16.80_R8, 17.50_R8, 18.20_R8, 18.90_R8, 19.60_R8, 20.30_R8,   &  
     &    21.00_R8, 21.70_R8/
 
      DATA ( bz(7,i), i = 1, NGTFTR )                               /    &  
     &      0.55_R8,  1.25_R8,  1.95_R8,  2.65_R8,  3.35_R8,  4.05_R8,   &  
     &     4.75_R8,  5.45_R8,                                            &  
     &      6.15_R8,  6.85_R8,  7.55_R8,  8.25_R8,  8.95_R8,  9.65_R8,   &  
     &    10.35_R8, 11.05_R8,                                            &  
     &     11.75_R8, 12.45_R8, 13.15_R8, 13.85_R8, 14.55_R8, 15.25_R8,   &  
     &    15.95_R8, 16.65_R8,                                            &  
     &     17.35_R8, 18.05_R8, 18.75_R8, 19.45_R8, 20.15_R8, 20.85_R8,   &  
     &    21.55_R8, 22.25_R8/
 
      DATA ( az(8,i), i = 1, NGUSR1 )                               /    &  
     &      0.00_R8,  0.70_R8,  1.40_R8,  2.10_R8,  2.80_R8,  3.50_R8,   &  
     &     4.20_R8,  4.90_R8,                                            &  
     &      5.60_R8,  6.30_R8,  7.00_R8,  7.70_R8,  8.40_R8,  9.10_R8,   &  
     &     9.80_R8, 10.50_R8,                                            &  
     &     11.20_R8, 11.90_R8, 12.60_R8, 13.30_R8, 14.00_R8, 14.70_R8,   &  
     &    15.40_R8, 16.10_R8,                                            &  
     &     16.80_R8, 17.50_R8, 18.20_R8, 18.90_R8, 19.60_R8, 20.30_R8,   &  
     &    21.00_R8, 21.70_R8/
 
      DATA ( bz(8,i), i = 1, NGUSR1 )                               /    &  
     &      0.55_R8,  1.25_R8,  1.95_R8,  2.65_R8,  3.35_R8,  4.05_R8,   &  
     &     4.75_R8,  5.45_R8,                                            &  
     &      6.15_R8,  6.85_R8,  7.55_R8,  8.25_R8,  8.95_R8,  9.65_R8,   &  
     &    10.35_R8, 11.05_R8,                                            &  
     &     11.75_R8, 12.45_R8, 13.15_R8, 13.85_R8, 14.55_R8, 15.25_R8,   &  
     &    15.95_R8, 16.65_R8,                                            &  
     &     17.35_R8, 18.05_R8, 18.75_R8, 19.45_R8, 20.15_R8, 20.85_R8,   &  
     &    21.55_R8, 22.25_R8/
 
      DATA ( az(9,i), i = 1, NGUSR2 )                               /    &  
     &      0.00_R8,  0.70_R8,  1.40_R8,  2.10_R8,  2.80_R8,  3.50_R8,   &  
     &     4.20_R8,  4.90_R8,                                            &  
     &      5.60_R8,  6.30_R8,  7.00_R8,  7.70_R8,  8.40_R8,  9.10_R8,   &  
     &     9.80_R8, 10.50_R8,                                            &  
     &     11.20_R8, 11.90_R8, 12.60_R8, 13.30_R8, 14.00_R8, 14.70_R8,   &  
     &    15.40_R8, 16.10_R8,                                            &  
     &     16.80_R8, 17.50_R8, 18.20_R8, 18.90_R8, 19.60_R8, 20.30_R8,   &  
     &    21.00_R8, 21.70_R8/
 
      DATA ( bz(9,i), i = 1, NGUSR2 )                               /    &  
     &      0.55_R8,  1.25_R8,  1.95_R8,  2.65_R8,  3.35_R8,  4.05_R8,   &  
     &     4.75_R8,  5.45_R8,                                            &  
     &      6.15_R8,  6.85_R8,  7.55_R8,  8.25_R8,  8.95_R8,  9.65_R8,   &  
     &    10.35_R8, 11.05_R8,                                            &  
     &     11.75_R8, 12.45_R8, 13.15_R8, 13.85_R8, 14.55_R8, 15.25_R8,   &  
     &    15.95_R8, 16.65_R8,                                            &  
     &     17.35_R8, 18.05_R8, 18.75_R8, 19.45_R8, 20.15_R8, 20.85_R8,   &  
     &    21.55_R8, 22.25_R8/
 
      DATA ( az(10,i), i = 1, NGUSR3 )                               /   &  
     &      0.00_R8,  0.70_R8,  1.40_R8,  2.10_R8,  2.80_R8,  3.50_R8,   &  
     &     4.20_R8,  4.90_R8,                                            &  
     &      5.60_R8,  6.30_R8,  7.00_R8,  7.70_R8,  8.40_R8,  9.10_R8,   &  
     &     9.80_R8, 10.50_R8,                                            &  
     &     11.20_R8, 11.90_R8, 12.60_R8, 13.30_R8, 14.00_R8, 14.70_R8,   &  
     &    15.40_R8, 16.10_R8,                                            &  
     &     16.80_R8, 17.50_R8, 18.20_R8, 18.90_R8, 19.60_R8, 20.30_R8,   &  
     &    21.00_R8, 21.70_R8/
 
      DATA ( bz(10,i), i = 1, NGUSR2 )                               /   &  
     &      0.55_R8,  1.25_R8,  1.95_R8,  2.65_R8,  3.35_R8,  4.05_R8,   &  
     &     4.75_R8,  5.45_R8,                                            &  
     &      6.15_R8,  6.85_R8,  7.55_R8,  8.25_R8,  8.95_R8,  9.65_R8,   &  
     &    10.35_R8, 11.05_R8,                                            &  
     &     11.75_R8, 12.45_R8, 13.15_R8, 13.85_R8, 14.55_R8, 15.25_R8,   &  
     &    15.95_R8, 16.65_R8,                                            &  
     &     17.35_R8, 18.05_R8, 18.75_R8, 19.45_R8, 20.15_R8, 20.85_R8,   &  
     &    21.55_R8, 22.25_R8/
!
      DATA pmax, xmin, xmax, anzac /                                     &  
     &      1.0_R8, -10._R8, +10._R8,  1.20_R8/
      DATA RelPha, MinNpar / 120._R8, 1.50_R8/
!
      INTEGER Curious, TRUE, FALSE
      DATA    Curious, TRUE, FALSE /                                     &  
     &              0,    1,     0 /
      REAL*8    RE81, RE82, RE8360
      REAL*8 AREAL

      RE8360 = 360._R8
!
!     Copy argument list into the local variables:
!
 1    continue
      f0 = fghz * 1000._R8
      nplot  = nslices
      NumRay = nrays
      MinNpar= abs(nparaMin)
      RE81   = phaseDeg
      RelPha = mod(RE81,    RE8360)
!     RelPha = amod(phaseDeg,360.)
      do 2 i = 1, DIM
        p0(i)=1.00_R8
        RE81    = AREAL(i)*RelPha
        phi0(i) = mod( RE81          , RE8360)
!       phi0(i) = amod(float(i)*RelPha, 360.)
 2    continue
 
      if      ( coupler .eq. 'PBXMFAST' )  then
        G = 1
        nguide = NGPBX
      else if ( coupler .eq. 'PBXMSLOW' )  then
        G = 2
        nguide = NGPBX
      else if ( coupler .eq. 'TOKDEVAR' )  then
        G = 3
        nguide = NGTdeV
      else if ( coupler .eq. 'TORSUPRA' )  then
        G = 4
        nguide = NGTors
        p0(1) =     -1.00_R8
        p0(NGTors)= -1.00_R8
 
        phi0(1) =    0.00_R8
        phi0(NGTors)=0.00_R8
        phi0(2) =    0.00_R8
        phi0(3) =   90.00_R8
        phi0(4) =  180.00_R8
        phi0(5) =  270.00_R8
        do i = 6, NGTors-1
          RE81    = phi0(i-4)+RelPha
          phi0(i) = mod( RE81,              RE8360)
!         phi0(i) = amod((phi0(i-4)+RelPha), 360.)
        enddo
      else if ( coupler .eq. 'SLOWSLOW' ) then
        G = 5
        nguide = NGSslow
      else if ( coupler .eq. 'JET_LHCD' ) then
        G = 6
        nguide = NGJET
        phi0(1) =    0.00_R8
        phi0(2) =   90.00_R8
        phi0(3) =  180.00_R8
        phi0(4) =  270.00_R8
        do i = 5, NGTors
          RE81    = phi0(i-4)+RelPha
          phi0(i) = mod( RE81,              RE8360)
!         phi0(i) = amod((phi0(i-4)+RelPha), 360.)
        enddo
      else if ( coupler .eq. 'TFTRLHCD' ) then
        G = 7
        nguide = NGTFTR
      else if ( coupler .eq. 'USRSPEC1') then
        G = 8
        nguide = NGUSR1
      else if ( coupler .eq. 'USRSPEC2') then
        G = 9
        nguide = NGUSR2
      else if ( coupler .eq. 'USRSPEC3') then
        G = 10
        nguide = NGUSR3
      else
        G = 1
        nguide = NGPBX
      endif
 
      if (nplot .gt. NPDIM) nplot = NPDIM
 
!
!         read in data
!      type 2000
! 2000 format(' enter data file: ')
!      accept 1000,filename
! 1000 format(a10)
!      open(unit=21,access='sequential',file=filename,
!     ^		status='old' )
!      read(21,2010) nguide,nd,
!     ^              delphi,f0,epsr,x0,  gradn,
!     ^              a,anmax,anplt,xprobe,
!     ^              edged
! 2010 format(/, 2i5,4g5.0,e10.4,4g5.0,g10.0 )
!		input max. for spectrum plot
!      write(6,2015)
! 2015  format(' Enter max spectral amplitude, npar(min), npar(max) ',
!     1       /' nz(accessible) to plot(<=0 for autoscaling):' )
!      read(5,*) pmax,xmin,xmax,anzac
!
!      if(abs(xmax).gt..01) go to 10
!      if(anzac.lt.1.0) anzac=1.0
!      xmin=-anplt
!      xmax=anplt
! 10   continue
!         check if phase iteration requested
      isw=0
      if(nd.gt.5001) nd=5001
      if(nguide.lt.0) isw=2
      if(nguide.lt.0) nguide=-nguide
      if(nguide.gt.idim) go to 999
!         iterate phase
      nphase=1
      if(delphi .gt. 14._R8) nphase = (180._R8/delphi) + 1.1_R8
!         include waveguide impedance?(iepsw=1)
      iepsw=1
      if(epsr.lt.0.0_R8) iepsw=0
      epsr=abs(epsr)
      if(epsr.lt.1.0_R8) epsr=1.0_R8
      dencr=f0**2/0.009_R8**2
!         iterate on edge density
      nedgd=1
      if(edged.lt.0.0_R8) nedgd=7
      edged0=abs(edged)
      if(nedgd.gt.1.and.(edged0/dencr.gt.111._R8.or.edged0/dencr         &  
     &     .lt.1._R8))                                                   &  
     & edged0=dencr
!         iterate density gradient
      ngradn=1
      if(gradn.lt.0.0_R8) ngradn=4
      gradn0=abs(gradn)
	iskip=0
	if(nphase.gt.1.or.ngradn.gt.1.or.nedgd.gt.1)iskip=1
!
!         read in data for each waveguide
      pave=0.0_R8
      ifast=1
      zmin=9999._R8
      zmax=-9999._R8
!      do 20 i=1,nguide
!      read(21,2020,END=18,err=999) p0(i),phi0(i),az(i),bz(i)
! 2020 format(6g10.0)
!	go to 19
!c	assign values based on wg no. 1
! 18   p0(i)=p0(i-1)
!      phi0(i)=phi0(i-1)+phi0(1)
!      az(i)=bz(i-1)+az(1)
!      bz(i)=bz(i-1)+bz(1)
! 19   continue
!c         check for passive guides
!      if(p0(i).lt.0.0) isw=1
!      if(p0(i).gt.0.0)pave=pave+p0(i)
!      if(az(i).lt.zmin) zmin=az(i)
!      if(bz(i).gt.zmax) zmax=bz(i)
!c         check for symmetry of the waveguides(ifast=1)
!      if(i.gt.1.and.abs((az(i)-az(i-1))-(bz(i)-bz(i-1))).gt.0.01)ifast=0
! 20   continue
!      pave=pave/nguide
!      twidth=abs(zmax-zmin)
      pave = 1._R8
      twidth = abs ( bz(G,nguide) - az(G,1) )
!
!         calculate constants
      nloop=1
      if(isw.gt.0) nloop=5
      c0=(0._R8,0._R8)
      c1=(0._R8,1._R8)
      w=2._R8*PI*f0*1.E6_R8
      wdc=w/3.E10_R8
      al0=3.E4_R8/f0
      beta=sqrt(epsr-(al0/2._R8/a)**2)
      beta1=beta
      if(iepsw.eq.0) beta1=sqrt(epsr)
      xp=c1*wdc*xprobe*beta
!
!         begin main loop
!
!         iterate on phase
      do 800 nph=1,nphase
      if(delphi.lt.14._R8) go to 26
      do 25 i=1,nguide
   25 phi0(i)=(nph-1._R8)*(i-1._R8)*delphi
   26 continue
!
!         iterate density gradient
      do 800 ngn=1,ngradn
      gradn=gradn0*10._R8**(ngn-1)
!
!         iterate edge density
      do 800 ned=1,nedgd
      edged=edged0*3.1622777_R8**(ned-1._R8)
      ratio=edged/dencr-1._R8
      if(ratio.lt.-0.999_R8) edged=dencr
      if(ratio.lt.-0.999_R8)ratio=0.0_R8
!          iterate for passive guides or phase feedback of density gradient
      do 788 nl=1,nloop
!         incident electric field in i'th guide in v/cm
      istop=1
      if(nl.eq.1) istop=0
      do 30 i=1,nguide
      b(i)=abs(az(G,i)-bz(G,i))
      if(nl.ne.1.and.p0(i).lt.0.0_R8.and.p0(i).le.p0r(i)) alp=0.9_R8
      if(nl.ne.1.and.p0(i).lt.0.0_R8.and.p0(i).gt.p0r(i)) alp=1.5_R8
      if(nl.ne.1.and.p0(i).lt.0.0_R8)p0(i)=alp*p0r(i)+(1._R8-alp)*p0(i)
      if(nl.eq.1.and.p0(i).lt.0.0_R8) p0(i)=-pave/10._R8
      if(nl.ne.1.and.p0(i).lt.0.0_R8.and.abs(rm(i)-1._R8)                &  
     &                                       .gt.0.01_R8)istop=0
      arg=4.E3_R8*abs(p0(i))*120._R8*PI/a/b(i)/beta
      phold(i)=phi(i)
      phi(i)=phi0(i)
      if(nl.eq.1.or.i.eq.1.or.isw.ne.2) go to 29
!     if(nl.eq.1.or.i.eq.1.or.isw.ne.2) go to 30
      ph0=phi0(i)-phi0(i-1)
      if(phi0(i).lt.phi0(i-1)) ph0=ph0+360._R8
      ph1=btph(i)-btph(i-1)
      if(btph(i).lt.btph(i-1)) ph1=ph1+360._R8
      ph2=phold(i)-phold(i-1)
      if(phold(i).lt.phold(i-1)) ph2=ph2+360._R8
      del=ph0-ph1
      idel=del/180._R8
      del=del-idel*180._R8
      if(abs(del).gt.1.0_R8)istop=0
      phi(i)=phi(i-1)+ph2+del*0.5_R8
      nmod=phi(i)/360._R8
      phi(i)=phi(i)-nmod*360._R8
   29 RE81 = PI*phi(i)/180._R8
   30 e0(i)=sqrt(arg)*exp(RE81*c1)
!  30 e0(i)=sqrt(arg)*exp(PI*c1*phi(i)/180.)
      if(istop.eq.1) go to 799
      alph3=(gradn/(f0/.009_R8)**2/wdc)**0.333333_R8
      RE81 = PI/3._R8
      c3=exp(RE81*c1)*a2/a1*alph3
!     c3=cexp(PI*c1/3.)*a2/a1*alph3
      RE81 = PI/2._R8
      c2=exp(RE81*c1)*a2/a1*alph3
!     c2=cexp(PI*c1/2.)*a2/a1*alph3
      c9=2._R8*c1*x0*wdc
!
!         zero matricies
      do 51 l=1,nguide
      r(l)=c0
!
        do 50 k=1,nguide
        ai(l,k)=c0
        aa(l,k)=c0
 50     continue
 51   continue
!
!         perform nz integration from -anmax to anmax by simpson's rule
      dnz=anmax/nd
      anz=0.0_R8
!         if guides are symmetrical,
!         then don't need to compute all the matrix elements
      ngset=nguide
      if(ifast.eq.1)ngset=1
      do 100 n=1,nd,2
      do 100 n2=1,2
      anzi=anz*c1*wdc
!         abs(nz) > 1
      if(abs(anz).le.1.0_R8) go to 60
      anx=c1*sqrt(anz*anz-1._R8)
      anx3=(anz**2-1._R8)**(0.1666667_R8)
      if(anx3.lt.0.01_R8) anx3=0.01_R8
      c5=c3/anx3
      if(ratio.le.0.1_R8) go to 70
!         density step present
      w0=-anx3**2*ratio/alph3**2
      call airy(air,bir,aip,bip,w0)
      c5=alph3*(aip+c1*bip)/(air+c1*bir)/anx3
      go to 70
!         abs(nz) < 1
   60 anx=sqrt(1._R8-anz**2)
      anx3=(1._R8-anz**2)**(0.166667_R8)
      if(anx3.lt.0.01_R8) anx3=0.01_R8
      c5=c2/anx3
      if(ratio.le.0.1_R8) go to 70
!         density step present
      w0=anx3**2*ratio/alph3**2
      call airy(air,bir,aip,bip,w0)
      if(abs(air).lt.0.001_R8) air=0.001_R8
      c5=-c1*alph3*aip/air/anx3
   70 denom=anx*anz*anz
      if(abs(denom).lt.0.001_R8) denom=0.001_R8
      rho=exp(c9*anx)*(1._R8-c5)/(1._R8+c5)
      c6=n2*(1._R8-rho)/(1._R8+rho)/denom
      do 80 l=1,ngset
      RE81 = bz(G,l)
      RE82 = az(G,l)
      c7=exp(anzi*RE81)   -exp(anzi*RE82)
!     c7=cexp(anzi*bz(G,l))-cexp(anzi*az(G,l))
      do 80 k=1,nguide
      RE81 = bz(G,k)
      RE82 = az(G,k)
      c8=c7*(exp(-anzi*RE81)    -exp(-anzi*RE82  ))
!     c8=c7*(cexp(-anzi*bz(G,k))-cexp(-anzi*az(G,k)))
   80 ai(l,k)=ai(l,k)+c6*(c8+conjg(c8))
  100 anz=anz+dnz
!
!         guides are symmetrical - use symmetry to compute matrix elements ai
      if(ifast.eq.0) go to 200
      do 150 l=2,nguide
      do 130 k=l,nguide
  130 ai(l,k)=ai(l-1,k-1)
      do 150 k=1,l-1
  150 ai(l,k)=ai(k,l)
!         multiply 'ai' by constant, set up matrix equation
  200 c10=dnz/3._R8/PI/wdc/beta1
      do 300 l=1,nguide
      bb(l)=e0(l)
      aa(l,l)=e0(l)
      do 300 k=1,nguide
      ai(l,k)=ai(l,k)/b(l)*c10
      bb(l)=bb(l)-e0(k)*ai(l,k)
  300 aa(l,k)=aa(l,k)+e0(k)*ai(l,k)
!       print 8000,((l,k,aa(l,k),bb(l),ai(l,k),k=1,nguide),l=1,nguide)
!  8000 format(2i5,6e12.4)
!
!         invert matrix 'aa'  -  put in array 'ai'
      call invert(idim,nguide,aa,ai)
!
!         multiply aa(inverse) * bb  =  r
      do 400 l=1,nguide
      do 400 k=1,nguide
!       print 8000,l,k,ai(l,k),aa(l,k)
  400 r(l)=r(l)+ai(l,k)*bb(k)
!
!         compute parameters in each guide
      pplas=0.0_R8
      pint=0.0_R8
      area=0.0_R8
      do 500 k=1,nguide
      if(p0(k).gt.0.0_R8)area=area+a*b(k)
      rm(k)=abs(r(k))
      rr=REAL(r(k))
      ri=AIMAG(r(k))
      rp(k)=atan2(ri,rr)*180._R8/PI
      if(p0(k).gt.0.0_R8)pint=pint+p0(k)
      p0r(k)=p0(k)*rm(k)**2
      if(p0(k).gt.0.0_R8)pplas=pplas+p0(k)-p0r(k)
      et(k)=e0(k)*(exp(xp)+r(k)*exp(-xp))
      ei(k)=abs(e0(k))
      etr=REAL(et(k))
      eti=AIMAG(et(k))
      etm(k)=abs(et(k))
      etph(k)=atan2(eti,etr)*180._R8/PI
      bt(k)=-e0(k)*(exp(xp)-r(k)*exp(-xp))
      btm(k)=abs(bt(k))
      btr=REAL(bt(k))
      bti=AIMAG(bt(k))
  500 btph(k)=atan2(bti,btr)*180._R8/PI
      teff=100._R8*pplas/PInt
      pdave=pint/area
!
!         print general information for run
      epsrp=epsr
      if(iepsw.eq.0)epsrp=-epsr
!	loop=1+(nguide+15)/16
!	do 777 lll=1,loop
 
!	ix=1
!	iy=765
!                                       Take out the write statements, normally
                                        if (Curious .eq. TRUE) then
	write(nTSCscrn,2100) pint,pplas
	write(nTSCscrn,2101) f0,nguide
	write(nTSCscrn,2102) gradn,x0
	write(nTSCscrn,2103) edged,twidth
	write(nTSCscrn,2104) pdave,teff
	write(nTSCscrn,2105) a,epsrp
	write(nTSCscrn,2106) dnz,anmax
	write(nTSCscrn,2107) xprobe,nl
 2100 format(' incident rf power(kW) =',f9.0,3x,                         &  
     & ' rf power to plasma(kw)=' ,f9.0)
 2101 format(' frequency(MHz)        =' ,f9.0,3x,                        &  
     & ' no. of waveguides     =' ,i9)
 2102 format(' density gradient(cm-4)=' ,1pe9.2,3x,                      &  
     & ' plasma-guide gap(cm)  =' ,0pf9.3)
 2103 format(' edge density(cm-3)    =' ,1pe9.2,3x,                      &  
     & ' total grill width(cm) =' ,0pf9.3)
 2104 format(' ave.pow.den.(kW/cm**2)=' ,f9.3,3x,                        &  
     & ' % efficiency          =' ,f9.0)
 2105 format(' waveguide height(cm)  =' ,f9.3,3x,                        &  
     & ' dielectric constant   =' ,f9.3)
 2106 format(' nz integration step   =' ,f9.3,3x,                        &  
     & ' maximum nz            =' ,f9.3)
 2107 format(' rf probe position(cm) =' ,f9.3,3x,                        &  
     & ' iteration             =' ,i9)
!
                                        endif
!
!         print results for each guide
!	ngpr1= 1+16*(lll-1)
!	ngpr2=16+16*(lll-1)
        ngpr1 = 1
        ngpr2 = nguide
	if(ngpr2.gt.nguide) ngpr2=nguide
	if(ngpr1.gt.nguide) go to 543
!                                       Take out more wri
!
                                        if (Curious .eq. TRUE) then
        write (nTSCscrn,2109)
        write (nTSCscrn,2110)
 2109    format(1x)
 2110    format(' guide  p0 p0ref mag(r) ph(r)',                         &  
     &          '   e0  phe0  etot  phet  btot  phbt b(cm)')
        do 520 i=ngpr1,ngpr2
        write(nTSCscrn,2120)                                             &  
     &                i,p0(i),p0r(i),rm(i),rp(i),ei(i),phi(i),           &  
     &                etm(i),etph(i),btm(i),btph(i),b(i)
 2120     format(i4,1x,f6.0,f5.0,f7.3,2f6.0,f5.0,4f6.0,f6.2)
 520    continue
                                        endif
!
 543    continue
!	if(lll.ne.loop) go to 776
!
!         plot power vs. nz
      smx=0.0_R8
      sum=0.0_R8
      t1=a*beta/beta1*1.E-4_R8/24._R8/wdc/PI**2
      c2=c1*t1*alph3
      c9=c1*x0*wdc
	danplt=anplt*2._R8/nplot
	npmin=1
	npmax=1
!               integrate 1/nz and 1/nz**2 over spectrum
	asum1=0.0_R8
	asum2=0.0_R8
	asum3=0.0_R8
	fracc=0.0_R8
 
      do 600 i=1,nplot
      anz=-anplt+(i-1)*anplt/(nplot*0.5_R8-1._R8)
      sn(i)=anz
        if(anz.lt.xmin+danplt)npmin=i
        if(anz.lt.xmax)npmax=i
      c8=c0
      if(abs(anz).lt.1.1_R8) go to 580
      anx=sqrt(anz**2-1._R8)*c1
      t2=1.0_R8/anz**2
      anzi=anz*c1*wdc
      c6=c0
      c7=c0
 
      do 550 k=1,nguide
      RE81 = az(G,k)
      RE82 = bz(G,k)
      c5=exp(-anzi*RE81)   -exp(-anzi*RE82)
!     c5=cexp(-anzi*az(G,k))-cexp(-anzi*bz(G,k))
      c7=c7+anx*beta1*e0(k)*(1._R8-r(k))*c5
!     c7=c7+anx*beta1*e0(k)*(1.-r(k))*c5
  550 c6=c6+e0(k)*(1._R8+r(k))*c5
 
      enz=(c6+c7)/4._R8/c1
      rho=(c6-c7)/(c6+c7)
      c8=enz*(exp(c9*anx)+rho*exp(-c9*anx))
      if(ratio.le.0.1_R8) go to 560
      w0=-(anz**2-1._R8)**0.333333_R8*ratio/alph3**2
      call airy(air,bir,aip,bip,w0)
      c5=c2*(aip-c1*bip)/(air-c1*bir)
      go to 570
  560 RE81 = PI/3._R8
      c5=c2*a2/a1*exp(-RE81*c1)
! 560 c5=c2*a2/a1*cexp(-PI*c1/3.)
  570 c8=c8*conjg(c8)*t2*c5/(anz**2-1._R8)**0.666667_R8
  580 sp(i)=REAL(c8)
      if(sp(i).gt.smx) smx=sp(i)
      sum=sum+sp(i)
      s(i)=sum
        sign=anz/abs(anz)
        delt=sign*sp(i)/anz**2
        asum1=asum1+delt
      if(abs(anz).lt.anzac) go to 600
        asum2=asum2+delt
        asum3=asum3+delt
        fracc=fracc+sp(i)
600   continue
!
!	ave. nz over all power nz>1
                asum1=abs(asum1/sum)**0.5_R8
!	ave. nz over accessible power |nz|>nza
                asum2=abs(asum2/fracc)**0.5_R8
!	ave. nz over accessible power, normalize to total power
                asum3=abs(asum3/sum)**0.5_R8
                fracc=fracc/sum
!
!         plot nz spectrum in the plasma
      RE81= smx
      lsmx=log10(RE81)
!     lsmx=log10(smx)
      t=10._R8**lsmx
      ll=smx/t+1
      smx=ll*t
	if(pmax.gt.1.0_R8) smx=pmax
!
!
!                call initia
!
!       call ezsets(120,920,50,480,xmin,xmax,0.,smx,1)
	ndif=npmax-npmin+1
!       call ezcurv(sn(npmin),sp(npmin),ndif)
 
      do 700 i=npmin,npmax
        s(i)=s(i)*smx/sum
 700  continue
!       call ezpnts(sn(npmin),s(npmin),ndif)
	ix1=xmax-xmin+.01_R8
        ix1=ix1/2
	ix2=2
	if(ix1.lt.5)ix2=10
!       call ezaxes(ix1,ix2,5,2)
!       call ezwrit(520,8,'nz$',1,0)
!       call ezwrit(5,225,'power(kw/nz)$',1,1)
!
		ix=650
		iy=485
		write(pline,2221) asum1
!                call EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2221 format('<|nz|**-2>t=',f6.3)
		write(pline,2222) asum2
!                call EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2222 format('<|nz|**-2>a=',f6.3)
		write(pline,2223) asum3
!                call EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2223 format('<|nz|**-2>at=',f6.3)
		write(pline,2224) anzac
!                call EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2224 format('  nz(acc)=',f7.3)
		write(pline,2225) fracc
!                call EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2225 format('fract(acc)=',f6.3)                                            
!
 776  if(nloop.eq.1.and.iskip.eq.0) then
!         call tinput(idum)
        endif
! 777    continue
  788 continue
 
!         make hard copy
 799    continue
!       call finish
!
!     Plot the spectrum as a serites of bars
!     call ezsets(120,920,100,600,xmin,xmax,0.,smx,1)
!     call ezbars(sn(npmin),sp(npmin),ndif, 'y')
!     call ezaxes(ix1,ix2,5,2)
!     call ezwrit(520,25,'nz$',1,0)
!     call ezwrit(25,350,'rel power$',1,1)
!     write(pline,'(''phase, nrays, nslic: '',f8.4,1x,i3,1x,i3,''$'')')
!    ^             phaseDeg, nrays, nslic
!     -----------------------------------------------------------------|
!     call ezwrit(120,700, pline, 0, 0)
 
!     call finish
 
  800 continue
!
!     Find the NumRay largest peaks
!     But first eliminate the low nparallels
      do 820 i= npmin, npmax
        if (abs(sn(i)) .le. MinNpar) sp(i) = 0.00_R8
 820  continue
!     NumRay from argument list; NumBeg labels the first of NumRays in stack
      NumBeg = npmax - (NumRay-1)
!
!     Put the highest sp's at top of the array
      call PikSr2NR ( ndif  , sp(npmin) , sn(npmin)  )
!
!     Take the NumRay highest sp's and put them in order
!     accorting to the sn's
      call PikSr2NR ( NumRay, sn(NumBeg), sp(NumBeg) )
!
!     Normalize so the NumRay highest sp's add to 1.00,
!     and find out what is the largest value
      TotPwr = 0.00_R8
      MaxPwr = 0.00_R8
!xxxx do 830 i = npmin + ndif - NumRay, npmin + ndif !Terpstra found this wrong
      do 830 i = NumBeg, npmax
        TotPwr = TotPwr + sp(i)
 830  continue
        ii = 1
!xxxx do 840 i = npmin + ndif - NumRay, npmin + ndif ! This one too!
      do 840 i = NumBeg, npmax
        sp(i) = sp(i)/TotPwr
        nparalls(ii) = sn(i)
        amplitds(ii) = sp(i)
        if(sp(i) .gt. MaxPwr) MaxPwr = sp(i)
        ii = ii+1
 840  continue
 
      if ( abs(RelPha - 180._R8) .gt. 10._R8.and.                        &  
     &     abs(RelPha -   0._R8) .gt. 10._R8.and.                        &  
     &     turnnegs             .eq. 1           ) then
        do 900 i = 1, NumRay
          nparalls(i) = abs(nparalls(i))
 900    continue
      endif
!
!
      return
  999 call LSCstop(' error in the Brambilla code')
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
