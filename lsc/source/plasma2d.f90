!
!     ------------------------------------------------------------------
!
!     Plasma2d                              ---------------------------|
      SUBROUTINE plasma2d(r,z,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
!
!     In which input cylindrical coords r,z (in m)
!     yield outputs psi, the poloidal flux (in webers per radian),
!     Br,Bz (in T) and RBphi (in T-m), and omc in GHz frequency,
!     Te (in keV),
!     \omega_{pe}^2 = pe2 in GHz^2,
!     \omega_{pi}^2 = pi2 in GHz^2,
!     \alpha_{ion}  = 3 \omega_{pi}^2/\omega^2 v_{Ti}^2/\omega^2
!     = aio (in inverse GHz^4)
!     \alpha_{elect}= 3/4 \omega_{pe}^2/\omega_{ce}^2
!     v_{Te}^2/\omega^2
!     = ael (in inverse GHz^4) [must divide by f_{ce}^4]
!
!
!     Fields are related to the flux by
!     Br = -(1/R) d psi/d Z
!     Bz =  (1/R) d psi/d R
!
!     A 'Bean' Equilibrium of TSC is used; up-down symmetric.
!     The grid data file goes from the midplane up.
!
!     Ignat       TSC     TSC
!     z-numbers   symm    full
!
!     .   .   .   .   .   .   .   .   . NumZ      nzp     nzp
!     .   .   .   .   .   .   .   .   .           nz      nz
!     .   .   .   .   .   .   .   .   .           nz-1    nz-1
!     .   .   .   .   .   .   .   .   .
!     .   .   .   .   .   .   .   .   .           3
!     mid plane _   _   _   _   _   _   _   _ NumZh     2     (nz+1)/2+1
!     .   .   .   .   .   .   .   .   .           3
!     .   .   .   .   .   .   .   .   .
!     .   .   .   .   .   .   .   .   . 3                 4
!     .   .   .   .   .   .   .   .   . 2         nz      3
!     .   .   .   .   .   .   .   .   . 1         nzp     2
!
!     1   2   3   4                 NumR      Ignat numbering
!     2   3   4                  nx  nxp      TSC   numbering
!
!
!     Flux array numbering: xsv
! psimin                          psilim
!    |                              |
!     .   .   .   .   .   .   .   .   .
!     1   2   3   4            NpsiM1 NpsiJ   Ignat numbering
!     2   3   4   5            npsit  npsip   TSC   numbering
!
!omment: From 000431 jardin steve@nersc (ppl) on 06/24/91 at 19:11:01
!	David:
!       TSC writes out for you the arrays
!       XSV(l),l=2,npsip       poloidal flux
!       anecc(l),l=2,npsip     electron density
!       tekev(l),l=2,npsip        "      temperature
!etc
!       These are all defined at the same locations.  Location 2
!       is near the center (but not exactly at the center) and
!       location npsip is just outside the plasma boundary.
!
!       TSC wants to read in the arrays
!
!       POWTSC(l),l=2,npsit
!       CURRTSC(l),l=2,npsit
!
!       so note that these quantities are defined at the same locations
!       as the ones written out, but that there is one less location read
!       than what is written.  This is because we don't read back in the
!       value of power and current density at location npsip just outside
!       the plasma.
!
!       note that      npsip = npsit+1
!                      npsim = npsit-1
!
!       depending on how you read your arrays in, you might store
!       and write out your stuff from   1,npsim   instead of 2,npsit.
!       If so , just keep in mind that your quantities are shifted
!       in the index 1 from TSC.
!                                             Steve (25jun91)
!
!     Identities:   NumZ = nz(full)
!     NumZh= nz(symm)
!     NumR = nx
!     NpsiJ = npsit = npsip-1
!     NumZh=(NumZ+1)/2
!     NumZ = 2*NumZh-1
!     Description of TSC variables:
!     (note--everything in MKS if not specified)
!
!     achrg:  charge of eqch ion species (hydrogen is 1)
!     amass:  mass of each ion species in grams (H mass is 1.6e-24 GRAM)
!     anecc:  electron density in particles/cc      --> eNeVec
!     anicc:       ion density in particles/cc      --> eNiVec
!     apl:       plasma current in amperes
!     bgzero:    vacuum field strength at rgzero
!     bgzero * rgzero == gzero in other routines
!     Btesl in common /RayWrk/
!     dVlVec: volume of a \psi-surface as passed from TSC: vptemp*(d\psi)
!     iVlVec: integral of volume from center to halfway to next indexed surface
!     gary:   toroidal field function R*Bt (MKS) -- G of psi
!     gpary:  derivative of gary wrt xsv --         G-prime of psi
!     iplim:  limiter switch   pos-plasma rests on limiter   neg-diverted
!     isym:   symmetry option,  0-no symmetry    1- up/down symmetry
!     kcycle: TSC time cycle number (starts at 0)
!     equhd:   80 character description of run from TSC title card
!     npsit:
!     npsitm: number of flux surfaces for surface averaged quantities
!     (note npsitm = npsit-1)
!     |guard point                            |guard point
!     .   .   .   .   .   .   .   .   .   .   .
!     1   2   3         nx-3    nx-1  nx nxp  nxp+1
!     ^                               ^
!     2                              npsit
!     read                             read
!
!     nspc:   total number of ion species
!     nx:     number of cartesian mesh points in x direction (for PsiGrd)
!     nz:     number of cartesian mesh points in z direction (for PsiGrd)
!     (called nxpl and nzpl in TSC write statements)
!     PARAMETERS:
!     PWORDS:  P       R   number of comment words read in
!     PIMP:     A     E    number of impurities  == nspc - 1
!     PPSI       R   T     number of psi values for fns of psi
!     PNX:        A E      number of x points in PsiGrd grid
!     PNZ:         M       number of z points in PsiGrd grid
!
!     pary:   pressure (MKS units) --               p of psi
!     ppary:  derivative of pressure wrt xsv --     p-prime of psi
!     psep,xsep,zsep:  flux value,x and z coordinates of separatrix (2)
!     PsiGrd:       poloidal flux per radian
!     psimin:\  poloidal flux per radian at magnetic axis and P/V boundary,
!     psilim:/  defined in common /RayWrk/
!     rgzero: nominal major radius of machine (= xplas)
!     Rmaj in common /RayWrk/
!     rho:    sqrt of the normalized *toroidal* flux
!     tekev:  electron temp in kev                  --> TeeVec
!     tikev:  ion temp in kev
!     times:  TSC problem time in sec
!     vptemp: derivative of volume wrt xsv
!     vptemp: note, this is centered derivative as of Sep93!!
!     vptemp: vol of shell = vptemp(i) * (xsv(i+1) - xsv(i-1))/2.
!     vptemp: name chosen in TSC to distinguish it as a derivative wrt poloidal
!             flux.  TSC mostly works on toroidal flux.
!     voltlp: toroidal loop voltage in volts
!     xary(2):leftmost boundary of cartesian mesh
!     xary(nxp): rightmost boundary of cartesian mesh
!     xmag,zmag:      (x,z) coordinates of magnetic axis
!     xsep:    see psep
!     xsv:     poloidal flux per radian
!     zary(2):   bottom boundary of cartesian mesh
!     zary(nzp):    top boundary of cartesian mesh
!     zsep:    see psep
!
!     PsiVec is a vector containing the psi values where quantities
!     are given for pressure and diamagnetism
!     ___Vec is explained in the following:
!     RBpVec is a vector containing R*Bphi (phi=toroidal): m-T
!     pe2Vec is a vector containing omega_{pe}^2 in GHz^2 frequency
!     pi2Vec is a vector containing omega_{pi}^2 in GHz^2 frequency
!     TeeVec is a vector of electron temperature in keV
!     AioVec is a vector of thermal correction from ions in GHz^4 frequency
!
      USE MKSetc
      USE params
      USE PIetc
      USE plcmx
      USE ProfBody
      USE RayWrk
      USE TSCgrap
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER RETURN_EVERYTHING, FORCE_RECALC, NO_FORCE, OUT_OF_RANGE
      PARAMETER(RETURN_EVERYTHING = 0, FORCE_RECALC = 1, NO_FORCE = 0,   &  
     &     OUT_OF_RANGE = -1)
      INTEGER iread, isave, jsave
      DATA isave, jsave                                                  &  
     &     / OUT_OF_RANGE, OUT_OF_RANGE /
!     initalize out of range
      REAL*8    psderiv(0:2, 0:2), psimat(0:3, 0:3)
      COMMON / plxxr / psderiv, psimat
!     COMMON / plxxi / isave, jsave ! commented out june17 92 ibm barfed
!     enforce static allocation
      REAL*8                                                             &  
     &     HUGE
      REAL*8                                                             &  
     &     r,z, psi, Br, Bz, RBphi, omc, Tee, pe2, pi2, aio, ael
      REAL*8                                                             &  
     &     Rold,   Zold, psiold,  Brold,  Bzold, omcold
!     DATA  achrg(1),   achrg(2),   achrg(3)        ! These presented here
!    ^     / 1.,         8.,         6./            ! for information purposes.
!           [1] H      [16] O      [12] C           ! Data not actually used
!     DATA  amass(1),   amass(2),   amass(3)        ! in the program.
!    ^     /1.6726e-24,  26.76e-24,  20.07e-24/     !
 
      DATA  HUGE,   iread                                                &  
     &     /1.E+30_R8,    0 /
 
!
!     ***************************************************************
      if (iread .eq. 1 ) go to 100
      iread = 1
      Rold  = 0._R8
      Zold  = 0._R8
      Psiold = 0._R8
      Brold  = 0._R8
      Bzold  = 0._R8
      omcold = 0._R8
!
!
 
 100  Continue
      If ( r .eq. Rold .and. z .eq. Zold ) then
         Psi    = psiold
         Br     = Brold
         Bz     = Bzold
         omc = omcold
      else
         Rold = r
         Zold = z
         call grapLSC(RETURN_EVERYTHING, r, z, psi, psderiv, NO_FORCE,   &  
     &        rmin, dr, zmin, dz, isave, jsave, nx, NumZ, PNX, PsiGrd,   &  
     &        psimat)
         Br = - psderiv(0, 1) / r
         Bz = psderiv(1, 0) / r
      endif
      call plasma1d (psi, psiold, RBphi,Tee,pe2,pi2,aio,ael)
      omc = OmcFac * SQRT( Br**2 + Bz**2 + (RBphi/r)**2 )
      Ael = Ael / omc**4
      Te = Tee
!
      Brold  = Br
      Bzold  = Bz
      omcold = omc
      return
      END
!
!
!     plasma2d ends                         ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
