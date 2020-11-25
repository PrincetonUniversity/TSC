!#include "f77_dcomplx.h"
!     E2byPr                                ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE E2byPr
      USE dielec
      USE Doflags
      USE FeBins
      USE MKSetc
      USE params
      USE PIetc
      USE PlPr
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE TSCgrap
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL DispRela, plasma2d, ivtabl
      INTEGER  ivtabl
!     E2byPr    Finds for all zones the ratio E_{z}^2/P;
!     damping constants  dlnPdsX  dlnPdsK;
!     n_{\parallel}; and polarization.
!     Written by D. W. Ignat, April 1991.
!     The time and path-position of entering and leaving
!     a zone is obtained by linear interpolation.  This is to allow
!     for best accuracy in obtaining  dt  and  ds .
!     This version does not anticipate crossing more than one zone
!     boundary in a step, but if this happens a warning is issued.
!
!     The approach is to average the quantity needed for each zone by
!     accumulating a total and dividing by the number of entries.  This
!     alone is sufficient for n_{\parallel} and polarization.
!     For E_z^2/P ... ezsq ... one needs to multiply by the time spent
!     in the zone, and divide by the volume of the zone.
!     For dlnPdsK  dlnPdsX we have chosen to multiply by   ds  at the end.
!     These variables are really, then,
!     dP/P -- Kernel to accept df/dv later; and dP/P -- assuming MaXwellian.
!
!     Zone boundaries are between the psi points, as follows.
!
!     +     .     +     .     +     .     +     .     +     .     +
! psimin                                                        psimax
!     1           2           3                     npsi-1        npsi
!  zone 1---><--zone 2--><--zone 3-->                        <--zone npsi
!
!     Ray power at the beginning is indexed 1.  As ray crosses into
!     the next zone, the value is indexed 2.  Power deposited in outer most
!     zone is therefore P_1 - P_2
!
!     sNew, sOld                New and old path length from start
!     tNew, tOld                and time duration from start.
!     RzindNew, RzindOld        The real and integer values of the zone
!     IzindNew, IzindOld        index.  These are related to the izind array.
!                               zind is small at the center and large at edge.
!                               izone=1 implies edge value of zind for
!                               the start of the ray.
!     Rzind* = (Psi - PsiMinx) / (PsiMaxx - PsiMinx) * (Npsi - 1) + 1.5
!     Rzind* = (Psi - PsiMinx) /  DelPsi                          + 1.5
!     Izind* = ifix(Rzind*)
!
!          Rzind*
!            ^
!     npsi   +
!            |    .                                  .
!     npsi-1 +
!            |         .    double jump              \
!     npsi-2 +            /                            double jump
!            |                                  .
!     .      +
!            |              .              .
!     .      +
!            |                   .    .
!     .      +
!            |
!     4      +
!   izind=3  |
!     3      +
!   izind=2  |
!     2      .____.____.____.____.____.____.____.____.____.
!                 0    1    2    3    4    5    6    7    8 > s (or) t
!
!
!     sSlope, tSlope    (New - Old)/(RzindNew - RzindOld)
!     sLeave(nzones), sEnter(nzones) are redundant: sLeave(i)=sEnter(i+1)
!     tLeave(nzones), tEnter(nzones) are redundant: tLeave(i)=tEnter(i+1)
!     accum(nRayQt)  an accumulator for averaging
!                    when the zone number is constant
!     RayQt(nRayQt)  the fresh parameter to be averaged
!     NinAc(nRayQt)  the number of entries now in the accumulator
!
!     RayQts are:
!     1 for ezsq    (epsz)^{-1} * 2.
! old 2 for dlnPdsK dD/dK33 /abs(dD/dk) * Im K33
! old 3 fzor dlnPdsX but the _K means kernel for adding df/dv later
! old               and the _X means MaXwellian assumed.
!     2 for dlnPdsK dD/dK33 /(w/2)(dD/dw) / (dV/dwt)f * Im K33
!     3 for dlnPdsX but the _K means kernel for adding df/dv later
!                   and the _X means MaXwellian assumed.
!     4 npar
!     5 epolX
!     6 epolY
!
!     CLIGHT    speed of light, in 10^8 m/s
!     ee        electron epsilon == 1 + \omega_{pe}^2/\omega^2
!     eps0      (\mu_o c^2)^{-1} farads per meter; epsilon sub zero
!     epsz      epsilon sub z; converts E_z^2 to energy density
!     ex        E_x / E_z ; x-polarization
!     ey      i E_y / E_z ; y-polarization
!     PI        3.1415926
!     cEparIK   constant converts df_e/dv to Im{K_{33}}
!               ImEpar == - PI \omega_{pe}^2/k_{\parallel}^2 df/dv
!                      == - cEparIK df/dv / kpar2
!
!
!     [ K + n n - n^2 I ]  \cdot E = 0
!       =   - -       =          -
!
!     so if E = [ ex , i ey , 1 ] E_z , then
!     ex = - (Kzz - nperp^2)/(nperp npar)
!     ey = ex Kxy / (Kyy - n^2)
!
!
!     (1/P) (dP/dt) = - 2 Im(K_par) dD/dK_par / (dD/d\omega)
!     (1/P) (dP/ds) = - 2 Im(K_par) dD/dK_par / abs(dD/dk)
!           (dP/dV) = - 2 Im(K_par) dD/dK_par / (dD/d\omega) \times
!                         (eps_z  /2) E_z^2
!
!     P =  U  \vec{v_g} \cdot \vec{A} where A is area of ray
!     P = <U>  (dV/dt)  on flux surface average
!
!     U = (1/2) (eps_0 /2) E^* \cdot [K + d/d\omega(\omega K)] \cdot E
!                          -          =                    =         -
!          ^        ^                 ^           ^
!          |                      magnetic     electric and particle
!      avg amplitudes
!
!       = (E_z^2/2) (eps_0 /2) [ ex, -i ey, 1] \cdot
!   -                      -
!   |  2 ee    i Kxy    0  |   ex
!   | -i Kxy   2 ee     0  | i ey
!   |    0       0      2  |   1
!   -                      -
!
!     where   ee = (1 + \omega_{pe}^2 / \omega_{ce}^2)
!
!     U = E_z^2  (eps_0 /2) [ 1 + ee (ex^2 + ey^2) - Kxy ex ey ]
!       = E_z^2  (eps_z /2)
!
!     E_z^2 /P  = [(eps_z /2) dV/dt]^{-1}        == ezsq
!     (note that dt is normalized with \omega....to get MKS put back.)
!
!     dVol(NPSIDIM) is an array such that dVol(j) is the volume
!     centered on psi_j, etc.  This is used in calculating
!     dV/dt.  Once a ray enters a zone, it is assumed to be spread out
!     over all the volume of that zone.
!
!     Manifestly evident consistency with the QL form is consistent
!     with
!     U = E_z^2 (eps_L /2)
!     where
!                eps_L/eps0 = \omega/2 dD/d\omega / (dD/dEpar)
!
!
!     For some special studies, the components along field,
!     perp to field and along gradient, perp to field and along
!     flux surface are important.  Call the unit vectors \
!     e_11, e_rr, e_tt respectively.
!     Related unit vectors are toroidal, along gradient, poloidal,
!     e_ph, e_rr, e_pl
!
!     The code is written in unit vectors
!     e_ph, e_R , e_Z, usually referred to as
!     (e_R , e_Z, e_ph) , (y(1), y(2), y(3)) , (y(4), y(5), y(6))
!
!     e_rr =       Bz/Bp e_R     -    Br/Bp e_Z +     0 e_ph
!     e_tt = Br/Bp Bph/B e_R  + Bz/Bp Bph/B e_Z + -Bp/B e_ph
!     e_11 =       Br /B e_R  +        Bz/B e_Z + Bph/B e_ph
!
!     e_pl =       Br/Bp e_R  +       Bz/Bp e_z +     0 e_ph
!
!     B    =  sqrt ( Br^2 + Bz^2 + Bph^2 )
!     Bp   =  sqrt ( Br^2 + Bz^2 )
      CHARACTER*70 ErrMsg
      INTEGER nRayQt, NinAc
      PARAMETER (nRayQt=6)
      REAL*8    RayQt(nRayQt), accum(nRayQt)
      INTEGER ifirst, i, jzn, jry,                                       &  
     &        IzindNew, IzindOld, IzindJmp
      REAL*8                                                             &  
     &        RzindNew, RzindOld, RzindCrs,                              &  
     &        sNew    , sOld    , tNew    , tOld    ,                    &  
     &        sSlope  ,           tSlope  ,                              &  
     &        sEnter(NZONDIM+1) , sLeave(NZONDIM)   ,                    &  
     &        tEnter(NZONDIM+1) , tLeave(NZONDIM)
      EQUIVALENCE (tEnter(2), tLeave(1))
      EQUIVALENCE (sEnter(2), sLeave(1))
      DATA    ifirst   /                                                 &  
     &          1      /
      REAL*8                                                             &  
     &        Btot, Bpol, Bphi,  dDdEpar, dDdEparOld,                    &  
     &        ee, eps0, ex, ey, woc,                                     &  
     &        Kpar, Kper, Kpar2,Kper2,                                   &  
     &        psie2,qpar,veow2
      REAL*8                                                             &  
     &        psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      REAL*8                                                             &  
     &        DispRela, det, dum
      REAL*8                                                             &  
     &        NparDead, NperDead, MaxwDead, MxdPdz, SMALL
      DATA    NparDead, NperDead, MaxwDead, MxdPdz, SMALL   /            &  
     &            9.8_R8,     147._R8,  1.0E-12_R8, 0.900_R8,            &  
     &            1.0E-30_R8/
!    ^            9.8 ,     147.,  1.0e-12, 0.001,  1.0e-30 / ! from 96
!     MxdPdz = 0.9 means the most Maxwellian Power Lost in a zone
!              is 10%....
 
      REAL*8    RE41, RE42
      REAL*8                                                             &  
     &        ZERO, ONE
      DATA    ZERO, ONE/                                                 &  
     &         0.0_R8, 1.0_R8/
      REAL*8 AREAL
!
      if (ifirst .eq. 1) then
!                                       Initialize constants; note omega fixed
        ifirst = 0
        cEparIK = 4._R8* PI * PI * ECOULB**2/ELECMS * 1.0E-08_R8
        eps0    = 1._R8/(4._R8*PI*1.0E-7_R8*CLIGHT*CLIGHT*1.0E+16_R8)
        omega   = fghz*2._R8*PI*1.0E+9_R8
        woc     = sqrt(woc2)
      endif
!
!     Begin initializations associated with starting a new ray --------|
!                                                                      |
!                                       lnewray is a flag set in DoRay
      if (lnewray .gt. 0) then
          lnewray = 0
          call plasma2d                                                  &  
     &            (y(1),y(2), psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
          tOld   = y(7)
          sOld   = y(8)
          Bpol  = sqrt (Br*Br + Bz*Bz)
          Bphi  = RBphi/y(1)
          Btot  = sqrt ( Bpol*Bpol + Bphi*Bphi)
          Kpar  =( y(4)*Br + y(5)*Bz            + y(6)/y(1)*Bphi)/Btot
!
          Kpar2 = Kpar**2
          Kper2 = Y(4)**2 + Y(5)**2 + (Y(6)/Y(1))**2   -  Kpar2
          Kper  = sqrt ( abs(Kper2) )
 
          izone  = 1
!         RzindOld = (Psi - PsiMinx)/DelPsi + 1.0 zones bounded by psi--archaic
!         RzindOld = (Psi - PsiMinx)/DelPsi + 1.5 zones centerd on psi
          RzindOld = (Psi - PsiMinx)/DelPsi + 1.5_R8
          RE41 = RzindOld
          IzindOld = int(RE41)
!         IzindOld = ifix(RzindOld)
 
          izind(izone,iray) = IzindOld
          sEnter(izone) = sOld
          tEnter(izone) = tOld
!                                       These are not needed, but are
!                                       plotted sometimes.
!         PowrRy(izone)   = 0.          !! prior to apr93
          PowrRy(izone)   = 1._R8
          RofRay(izone)   = y(1)
          ZofRay(izone)   = y(2)
          PofRay(izone)   = y(3)
          NparRy(izone)   = Kpar/woc
          NperRy(izone)   = Kper/woc
          rtPsRy(izone)   = sqrt( (psi-psimin)/(psilim-psimin) )
          TimeRy(izone)   = y(7)
          DistRy(izone)   = y(8)
          NeofRy(izone)   = pe2/Pe2Fac * 1.0E+14_R8
          BthRay(izone)   = sqrt(Br**2+Bz**2)
          BphRay(izone)   = RBphi/RofRay(izone)
 
          det = DispRela ( y(1) , y(2) , y(4) , y(5) , y(6)  )
          d1 = abs(d1)
          d2 = abs(d2)
          d4 = abs(d4)
          det = det / max( d1, d2,d4 )
          DetrRy(izone)   = det
 
!                                       The quantities to be collected
!                                       (RayQt) are averaged over the
!                                       zone by summing into accum,
!                                       dividing by NinAc.  This clears.
          do 50 i=1,nRayQt
            accum(i) = 0._R8
            RayQt(i) = 0._R8
 50       continue
          NinAc = 0
          dtdV  = 0._R8
          return
 
      endif
!                                                                      |
!     End   initializations associated with starting a new ray --------|
!
!                                       Find location and various parameters.
!                                       This is the NORMAL BEGINNING point, in
!                                       that all initializations have been made
!                                       and the ray is progressing through
!                                       zones.
!
      call plasma2d (y(1),y(2), psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      Bpol  = sqrt (Br*Br + Bz*Bz)
      Bphi  = RBphi/y(1)
      Btot  = sqrt ( Bpol*Bpol + Bphi*Bphi)
      Kpar  = ( y(4)*Br + y(5)*Bz            + y(6)/y(1)*Bphi)   / Btot
!
      Kpar2 = Kpar**2
      Kper2 = Y(4)**2 + Y(5)**2 + (Y(6)/Y(1))**2   -  Kpar2
      Kper  = sqrt ( abs(Kper2) )
      Qpar = Kpar2 - woc2*Eper
      dDdEparOld = (Qpar*(Qpar + Kper2) - woc4*Exy**2)
      dDdEpar    = (Qpar + kper2)*kpar2*kper2 / ( kper2 - woc2*Epar )
 
      EparI = 0._R8
      veow2 = 0.0445E-04_R8*tee/fghz**2
      psie2 = 2._R8*veow2*kpar2
!                                       Trying to avoid overflows here.
      if ( psie2 .gt. 0.02_R8) then
        psie2 = 1._R8/psie2
        EparI = 2._R8*RTPI*epsq*psie2*sqrt(psie2)*exp(-psie2)
      endif
!
!     dDdEparOld = (Qpar*(Qpar + Kper2) - woc4*Exy**2)  is the term
!     multiplying the Epar or K_{33} or K_{zz} in D.
!     This is not manifestly positive, but it is believed to be positive.
!     Now called Old.  We have had trouble with this evaluating wrong. 24apr92
!
!     Using the dispersion relation, dDdEpar can be rewritten
!     dDdEpar = (Qpar + kper2) kpar2 kper2 / ( kper2 - woc2 Epar )
!
!     dDdkABS is partial of D wrt vector k, abs value
!     The - is because of damping; the 2. is because power is field ^2
!
!     EparIK = cEparIK / Kpar2
      ex    = - (woc2*Epar - Kper2)/(Kper*Kpar)
!     ey    = ex*Exy / ( Eper- (Kper2+Kpar2)/woc2 ) equivalent to following:
      ey    =(ex*(Eper - Kpar2/woc2) + Kper*Kpar/woc2 ) / Exy
      ee    = 1._R8+ epsq/ecyc2
      epsz  = eps0*( 1._R8+ ee*(ex*ex+ey*ey) - Exy*ex*ey )
!     epQL  =        0.5 * wdDdw / dDdEparOld ! changed form 25Apr92
      epQL  =        0.5_R8* wdDdw / dDdEpar
      epsL  = eps0*  epQL
 
!     RayQt(1)            = 2./epsz
      RayQt(1)            = 2._R8/epsL
!     RayQt(2)            = + 2.  * dDdEpar  / dDdkABS * EparIK
!     RayQt(3)            = - 2.  * dDdEpar  / dDdkABS * EparI
!     RayQt(2)            = + EparIK / epQL jul7 92
      RayQt(2)            = + 1.00_R8/ epQL
      RayQt(3)            = - EparI  / epQL
      RayQt(4)            = Kpar/woc
      RayQt(5)            = ex
      RayQt(6)            = ey
 
      tNew   = y(7)
      sNew   = y(8)
!
      RofRay(izone)   = y(1)
      ZofRay(izone)   = y(2)
      PofRay(izone)   = y(3)
      NparRy(izone)   = Kpar/woc
      NperRy(izone)   = Kper/woc
      rtPsRy(izone)   = sqrt( (psi-psimin)/(psilim-psimin) )
      TimeRy(izone)   = y(7)
      DistRy(izone)   = y(8)
      NeofRy(izone)   = pe2/Pe2Fac * 1.0E+14_R8
      BthRay(izone)   = Bpol
      BphRay(izone)   = Bphi
      det = DispRela ( y(1) , y(2) , y(4) , y(5) , y(6)  )
      d1 = abs(d1)
      d2 = abs(d2)
      d4 = abs(d4)
      det = det / max( d1, d2,d4 )
      DetrRy(izone)   = det
!
!
      RzindNew = (Psi - PsiMinx)/DelPsi + 1.5_R8
      RE41     =  RzindNew
      IzindNew = int(RE41)
!     IzindNew = ifix(RzindNew)
!                                       START major IF/ELSE/ENDIF branch.
!                                       If we are in the same zone as
!                                       before, put RayQt in accumulator.
      if ( IzindNew .eq. IzindOld ) then
        do 100 i=1,nRayQt
          accum(i)=accum(i)+RayQt(i)
 100    continue
          NinAc = NinAc + 1
!                                       We are in a new zone.
!                                       Divide the accumulator by # entries
!                                       but take care if NinAc = 0.
!                                       Interpolate to find crossings.
!                                       The ELSE of the major IF/ELSE/ENDIF branch
      else
!
        do 101 i=1,nRayQt
          if (NinAc .eq. 0) then
            accum(i)=RayQt(i)
          else
            accum(i)=accum(i)/AREAL(NinAc)
          endif
 101    continue
        NinAc = 0
 
        sSlope   = ( sNew - sOld ) / ( RzindNew - RzindOld )
        tSlope   = ( tNew - tOld ) / ( RzindNew - RzindOld )
        RE41     = RzindNew
        IzindNew = int(RE41)
!       IzindNew = ifix(RzindNew)
!                                       The crossing point is different
!                                       depending on whether the ray is going
!                                       outward:   IzindNew > IzindOld
!                                       or inward: IzindNew < IzindOld
            if (IzindNew .gt. IzindOld) then
                RzindCrs = AREAL(IzindOld) + 1._R8
            else
                RzindCrs = AREAL(IzindOld)
            endif
        sLeave(izone) = sOld + sSlope * (RzindCrs - RzindOld)
        tLeave(izone) = tOld + tSlope * (RzindCrs - RzindOld)
!
!                                       EQUIVALENCE accomplishes this:
!                                       sEnter(izone+1) = sLeave(izone)
!                                       tEnter(izone+1) = tLeave(izone)
!
!                                       Compute the desired parameters
!                                       for the zone we just left.
!                                       npar according to the velocity table 15Jan93
 
        npar(izone,iray)    = accum(4)
        ivind(izone,iray) = ivtabl(npar(izone, iray))
        npar(izone,iray)    = 1._R8/vpar(ivind(izone,iray))
        izind(izone,iray) = IzindOld
 
        ezsq(izone,iray)    = accum(1) *                                 &  
     &                        (tLeave(izone)-tEnter(izone)) / omega /    &  
     &                         dVol(IzindOld)
        dlnPdsK(izone,iray) = accum(2) * cEparik /                       &  
     &                        (npar(izone,iray) * woc)**2 *              &  
     &                        (tLeave(izone)-tEnter(izone))
        dlnPdsX(izone,iray) = accum(3) *                                 &  
     &                        (tLeave(izone)-tEnter(izone))
 
 
        if (ezsq(izone,iray) .lt. 0._R8) then
          write(ErrMsg,'('' Ez2<0!zn ind iry epsL epsZ:''                &
     &    , i4,i4,i3,1pe10.2,1x,1pe10.2)')                               &
     &    izone, IzindNew, iray, epsL, epsZ
          call LSCwarn( ErrMsg)
          ezsq(izone,iray)   = 0.0_R8
          dlnPdsK(izone,iray)= 0.0_R8
          dlnPdsX(izone,iray)= 0.0_R8
        endif
        dtdV                = (tLeave(izone)-tEnter(izone)) /            &  
     &                        ( abs(dVol(izindOld)) + SMALL )
!                                       Dont forget to save the integers.
 
 
!                                       We are done with the information
!                                       accumulated for the old zone, so
!                                       set the accumulator to the value
!                                       of the parameters we just calc-
!                                       ulated.  Reset counter.
        do 102 i=1,nRayQt
          accum(i)=RayQt(i)
 102    continue
          NinAc = 1
 
!                                       If we have more zones to go, then
!                                       INCREMENT THE ZONE COUNTER, and
!                                       fill the arrays with the same
!                                       values we just found, in case
!                                       we run out of time steps (eg).
        if (izone .lt. nzones) then
          izone = izone + 1
          RofRay(izone)        = RofRay(izone-1)
          ZofRay(izone)        = ZofRay(izone-1)
          PofRay(izone)        = PofRay(izone-1)
          NparRy(izone)        = NparRy(izone-1)
          NperRy(izone)        = NperRy(izone-1)
          rtPsRy(izone)        = rtPsRy(izone-1)
          NeofRy(izone)        = NeofRy(izone-1)
          BthRay(izone)        = BthRay(izone-1)
          BphRay(izone)        = BphRay(izone-1)
 
          ezsq   (izone,iray)  = ezsq   (izone-1,iray)
          dlnPdsK(izone,iray)  = dlnPdsK(izone-1,iray)
          dlnPdsX(izone,iray)  = dlnPdsX(izone-1,iray)
          npar   (izone,iray)  = npar   (izone-1,iray)
          ivind  (izone,iray)  = ivind  (izone-1,iray)
          dum = 1._R8+dlnPdsX(izone-1,iray)
          dum = max(MxdPdZ,dum)
          dum = min(dum,ONE)
 
          PowrRy (izone)       = PowrRy(izone-1)*dum
          izind  (izone,iray)  = IzindNew
        else
!                                       But, if we used up all the zones,
!                                       quit as fast as you can.
          izone = nzones
          Lstop = 1
        endif
!
!     .                                 If Npar is so large that damping would
!     .                                 be total, stop the calculation now.
        if ( abs(npar(izone,iray)) .gt. NparDead ) then
          Lstop = 1
        endif
!     .                                 If Nper is really large, stop now.
!     .                                 No matter what.
          if ( abs(nperRy(izone))    .gt. NperDead ) then
            Lstop = 1
          endif
!     .
!     .                                 If all power is gone in Maxwellian
!     .                                 stop the calculation now.
          if ( PowrRy(izone) .lt. MaxwDead ) then
            Lstop = 1
            if (PrFlg(RAYWR) .ge. TRUE) then
              write(nLSCcomm,'('' Maxwellian power=0; stop early'')')
            endif
          endif
      endif
!                                       END   major IF/ELSE/ENDIF branch.
 
!
!                                       See if we jumped more than
!                                       one zone.  Report if so.
      IzindJmp = IzindNew - IzindOld
      if (IzindJmp*IzindJmp .gt. 1) then
         if (PrFlg(RAYWR) .ge. TRUE ) then
            write(nLSCcomm,'('' Jumped more than one zone;'',            &  
     &           '' izone, IzindNew, jump = '',                          &  
     &             i4, 1x, i4, 1x, i4)')izone, IzindNew, IzindJmp
!                                       If we jumped more than one zone
!                                       and the ray is well advanced, then
!                                       there is probably trouble coming.
         endif
        if(izone .ge. nzones/10) then
         call LSCwarn( ' More than 1 zone jumped and izone >> 1 ')
        endif
      endif
!                                       Set the  _Old parameters to the
!                                       existing _New ones.
      sOld  = sNew
      tOld  = tNew
      IzindOld = IzindNew
      RzindOld = RzindNew
!
      return
      END
!
!                                                                      |
!     E2byPr ends                           ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok
