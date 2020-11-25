C     RayBins                               ---------------------------|
C                                                                      |
C                                                                      |
c     nzones  INTEGER
c             number of zones (input)
c     nrays   INTEGER
c             number of rays (input)
c     npsiJ   INTEGER
c             number of spatial grid points (input from TSC);
c             .le. PPSI which sets the max array size of TSC.
c             The final  J  is for Jardin.  Equivalent to his   npsit  .
c     npsi    INTEGER
c             number of psi bins as used in Dql calculation;
c             can be more or less than npsii from the TSC grid.
c             .le. NPSIDIM
c             NOT to be confused with: NpsiJ, npsit
c     izone   INTEGER
c             zone index; low --> starting of ray; high --> ending
c     iray    INTEGER
c             ray index; increases from low starting n_{\parallel}
c     iLastRy INTEGER
c             index of the last ray traced; used when tracing rays gradually.
c     ipsi    INTEGER
c             psi index; increases from low psi at center to psimax at edge
c     ivind   INTEGER RAYARRY
c             table of index into Vpar along each ray
c     izind   INTEGER RAYARRAY                                           ejv
c             izind(izone, iray) maps ray intersection number
c             to spatial index
c             As a ray progresses from the outside in, izone increases
c             by one each time a psi zone is crossed.  Typically, izind( , ):
c
c                    RAY#1   RAY#2   RAY#3   RAY#4   RAY#5   RAY#6
c             ZONE#1   61      61      61      61      61      61
c             ZONE#2   60      60      60      60      60      60
c             ZONE#3   59      59      59      59      61      59
c             ZONE#4   58      58      58      60      --      58
c             ZONE#5   57      57      59      61              57
c             ZONE#6   56      58      60      --              56
c             ZONE#7   55      59      61                      57
c
c     iznew
c     izold   are not arrays, but are related to izind( , ).  As the ray
c             progresses, iznew is computed and compared to izold. If there
c             is a change, then a zone has been passed, and information
c             is stored.
c     power   REAL    RAYARRAY
c             the power on each ray as a function of zone number
c             Indexing: (izone, iray)
c             psi, etc at the zones intersection
c             = psi(izind(izone, iray))
c             Units:  Watts
c     dlnPdsK REAL    RAYARRAY
c             d ln(P)/ds for unit value of { \p f_e/\p v_{\parallel} }
C     dlnPdsX REAL    RAYARRAY
C             d ln(P)/ds if maXwellian, or at maXimum
C                                          without quasilinear burnthrough.
C             = -2 (\p D/\p Epar) / dDdkABS * Im{Epar}
C                where Epar is K_{zz} is K_{33}
C             Units: m^{-1}
C     dlnPds  REAL    RAYARRAY
C             d ln(P)/ds after quasilinear effects,
C             always larger than dlnPdsX by the thermal slope to ql slope
C             Units: m^{-1}
c     ezsq    REAL    RAYARRAY
c             (izone, iray)
c             ratio of square of parallel electric field
c             to power flowing along ray
c             Units: MKS -- (volts/meter)^2 / watts
c     npar    REAL    RAYARRAY
c             (izone, iray)
c             c_light k_parallel / omega (parallel index)
c     dnpar   REAL
c             scalar n_parallel with for each ray
c             used in Dql computation; computed from input
c     vpar    REAL    RAYARRAY
c             (izone, iray)
c             parallel phase velocity.  Used in q-l calculation.
c             Units: v_parallel / c_light
c     spectrumREAL
c             function which returns normalized power given n_parallel.
c             integral (spectrum(n) * dn) = 1.
c     spwidth REAL
c             scalar n_parallel spectral width (input variable).
c             Used by spectrum
c     spcentr REAL
c             central n_parallel value (input).
c             Used by spectrum
c     npargr  SUBROUTINE
c             function which generates initial nparallel grid
c             of size nray extending from nparmin to nparmax
c     nparmin REAL
c             initial minimum nparallel actually assigned to a ray
c     nparmax REAL
c             initial maximum nparallel actually assigned to a ray
c     rFudgDmp REAL
c             (izone, iray)    --- real to Fudge Damping
c             factor less than 1 showing the average power in a zone.
c             This factor is created in the ray picture, and is
c             communicated to the QL picture. I
c             If d is the damping decrement, ie, exp(-d), then
c             rFudgDmp = < exp ( -d ) >
c             rFudgDmp = (1 - exp(-d) ) / d  =~ 1 - d/2 + d^2/6 for d small
c
c     P(iz+1,ir) = P(iz,ir) exp ( -d )
c     Dql = P(iz,ir) \eta  v /( dvsym dVol )  * rFudgDmp
c     and, hopefully, Pql == Pray
c     note, \eta is the polarization term .... dD/dEpar / (omega dD/dw ) with
c     constants added.
c
c     Here follow quantities of interest for diagnostics and understanding
c     but not required for the quasilinear calculation:
c     RofRay (NZONDIM) radius location of ray vs zone number
c     ZofRay (NZONDIM) z (height)                " "
c     PofRay (NZONDIM) Phi (toroidal angle)      " "
c     NperRy (NZONDIM) n_{\perp}                 " "
c     NparRy (NZONDIM) n_{\parallel}             " "
c     rtPsRy (NZONDIM) {\psi}^0.5                " "
c     NeofRy (NZONDIM) n_e                    vs zone number
c     BthRay (NZONDIM) B_\theta               vs zone number
c     BphRay (NZONDIM) B_\phi                 vs zone number
c     PowrRy (NZONDIM) exp(-2\int k_i \cdot dr)  " "  (linear damping)
c     TimeRy (NZONDIM) time given as \omega t    " "
c     DistRy (NZONDIM) distance in meters        " "
c     DetrRy (NZONDIM) determinate D over largest term (relative error) ""
c     Diagnostic quantities end.
      INTEGER
     ^        ipsi, iray, iLastRy, iznew, izold,
     ^        izone,
     ^        lnewray,
     ^        npsi, nrays,
     ^        nzones, nrl, nru, nzl, nzu, nslices,
     ^        izind(NZONDIM, NRAYDIM), ivind(NZONDIM, NRAYDIM),
     ^        izcind(NZONDIM, NRAYDIM)
 
      REAL*8
     ^        dnpar, dtdV,
     ^        ezsq   (NZONDIM, NRAYDIM),
     ^        npar   (NZONDIM, NRAYDIM),
     ^        power  (NZONDIM, NRAYDIM),
     ^        dlnPds (NZONDIM, NRAYDIM),
     ^        dlnPdsK(NZONDIM, NRAYDIM),
     ^        dlnPdsX(NZONDIM, NRAYDIM), rFudgDmp(NZONDIM,NRAYDIM)
 
      REAL*8
     ^        RofRay (NZONDIM), ZofRay (NZONDIM), PofRay (NZONDIM),
     ^        NparRy (NZONDIM), NperRy (NZONDIM), rtPsRy (NZONDIM),
     ^        PowrRy (NZONDIM),
     ^        TimeRy (NZONDIM), DistRy (NZONDIM), DetrRy (NZONDIM),
     ^        NeofRy (NZONDIM), BthRay (NZONDIM), BphRay (NZONDIM)
 
      COMMON /rayIbin/
     ^        ipsi, iray, iLastRy, izind, izcind, ivind,
     ^        iznew, izold, izone,
     ^        lnewray,
     ^        npsi, nrays,
     ^        nzones, nrl, nru, nzl, nzu, nslices
      COMMON /rayRbin /
     ^        dnpar           , dtdV            ,
     ^        ezsq            ,
     ^        npar            ,
     ^        power           ,
     ^        dlnPds          , dlnPdsK         , dlnPdsX         ,
     ^        rFudgDmp        ,
     ^        RofRay          , ZofRay          , PofRay          ,
     ^        NparRy          , NperRy          , rtPsRy          ,
     ^        PowrRy          ,
     ^        TimeRy          , DistRy          , DetrRy          ,
     ^        NeofRy          , BthRay          , BphRay
c
c     rayIbin   common block for integer output
c     rayRbin   common block for real output arrays
c
c                                                                      |
c     RayBins ends                          ---------------------------|
 
 
 
 
 
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
