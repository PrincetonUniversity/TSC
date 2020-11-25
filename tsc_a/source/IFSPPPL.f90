      subroutine IFSPPPL( znine, zncne, znbne, zrlt, zrln,               &  
     &                    zq, zshat, zeps, zkappa, omegaexb,             &  
     &                    ne19, tekev, tikev, rmajor, btesla,            &  
     &                    switches, grhoi, gvti, gnu, gtau,              &  
     &                    chii, chie,                                    &  
     &                    zchiicyc, zchii1, zchii2, zchie1, zchie2,      &  
     &                    zrlt1, zrlt2, ierr )
!
!
!
!     switches(j) settings:
!
!        switches(1) = 0  No output
!                    = 1  Print page of inputs, computational steps and output
!        switches(2) = 0  Use first set of dimensioned input
!                         (ne19, tekev, tikev, btesla)
!                    = 1  Use second set of dimensioned/dimensionless input
!                         (grhoi, gvti, gnu, gtau)
!                         Note that, with either approach, rmajor is required.
!        switches(3) = 0  Use 1995 IFS/PPPL model
!                    = 1  Use 1994 IFS/PPPL model (IAEA paper)
!        switches(4) = 0  Use the value of gnu as given (if applicable)
!                    = 1  Multiply gnu by 0.84; this allows direct comparison
!                         with the fortran code provided by B. Dorland.
!                         There is a factor of 0.84 difference between the
!                         definition of collisionality in the 1995 published
!                         paper and the coded implementation provided by
!                         B. Dorland
!        switches(5) = 0  Do not relax restrictions on the value of znu
!                    = 1  Allow the value of znu to be larger than 10.0
!
!        switches(30), switches(31), switches(32)
!                       Reserved as output channel numbers
!
!        switches(j) should be defaulted to zero.
!
!     Dimensionless input parameters:
!        znine:     Ratio of thermal deuterium density to electron density
!        zncne:     Ratio of thermal carbon density to electron density
!        znbne:     Ratio of beam deuterium density to electron density
!        zrlt:      Ratio of major radius to ion temperature scale length
!        zrln:      Ratio of major radius to electron density scale length
!        zq:        Plasma safety factor
!        zshat:     Magnetic shear (defined as r/q dq/dr)
!        zeps:      Local inverse aspect ratio
!        zkappa:    Local plasma elongation
!        gtau:      Ratio of Ti/Te (used only when switches(2)=1)
!        gnu:       Local collisionality (used only when switches(2)=1)
!
!     Dimensioned input parameters:
!        ne19:      Electron density, in 10^19 meters^-3
!        tekev:     Electron temperature, in keV
!        tikev:     Deuterium ion temperature, in keV
!        rmajor:    Major radius, in meters
!        btesla:    Toroidal magnetic field, in Tesla
!        grhoi:     Local ion gyroradius, in meters (used when switches(2)=1)
!                   Input to replace calculated value for vthermi/omegaci
!        gvti:      Ion thermal speed, in m/s (used only when switches(2)=1),
!                      defined as sqrt(Ti/mi) -- nonstandard definition
!        omegaexb:  The E x B shearing rate, in inverse seconds
!
!     Dimensioned outputs:
!        chii:      Ion thermal diffusivity, in m^2/sec
!        chie:      Electron thermal diffusivity, in m^2/sec
!
!     Dimensionless outputs:
!        zchiicyc:  \chi_i, normalized by L_n/(\rho_i^2 v_{ti})
!        zchii1:    Hydrogenic \chi_i, normalized by R/(\rho_i^2 v_{ti})
!        zchii2:    Carbon-driven \chi_i, normalized by R/(\rho_i^2 v_{ti})
!        zchie1:    Hydrogenic \chi_e, normalized by R/(\rho_i^2 v_{ti})
!        zchie2:    Carbon-driven \chi_e, normalized by R/(\rho_i^2 v_{ti})
!        zrlt1:     Normalized threshold gradient for hydrogenic mode
!        zrlt2:     Normalized threshold gradient for carbon-driven mode
!        ierr:      Error code = 0 if all is OK
!                              = 1 if inputs are outside the ranges prescribed
!                                    in Kotschenreuther, et al,
!                                    Physics of Plasmas, v2, p2381 (1995)
!
!
!     Declare variables
!     NAMING CONVENTION: Dimensionless variables begin with a 'z'
!                        Dorland-style input variables begin with a 'g'
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER switches(32), ierr
      REAL*8 znine, zncne, znbne, zrlt, zrln, zq, zshat, zeps,           &  
     &       ne19, tekev, tikev, rmajor, rhoi, vthermi, omegaexb,        &  
     &       ztau, znu, zzeffstr, zrlnstr, ztaub, zkappa,                &  
     &       zchii, zchii1, zchii2, zscriptg1, zscriptg2, zx,            &  
     &       zscriptz, zf, zg, zh, zscriptd, zscripte,                   &  
     &       zrlt1, zrlt2, zchie, zchie1, zchie2, zscriptj,              &  
     &       btesla, chie, chii, omegaci, zcorrect, zkapfact,            &  
     &       grhoi, gvti, gtau, gnu, zscripty, gamma,                    &  
     &       zc1, zc2, zc3, zscriptg, zchiicyc, zscripts
      CHARACTER*35 calcflag
!
!
!     Set some initial values
!
      ierr = 0
      omegaexb = abs(omegaexb)
      calcflag = '                                   '
      if ( switches(2) .eq. 1)                                           &  
     &     calcflag = '* Calcd from gnu,gvti,gtau,grhoi * '
!
!
!     Check some of the inputs for validity
!
      if (zq .lt. 0.7_R8) then
         zq = 0.7_R8
         ierr = 1
      end if
      if (zq .gt. 8.0_R8) then
         zq = 8.0_R8
         ierr = 1
      end if
      if (zshat .lt. 0.5_R8) then
         zshat = 0.5_R8
         ierr = 1
      end if
      if (zshat .gt. 2.0_R8) then
         zshat = 2.0_R8
         ierr = 1
      end if
      if (zeps .lt. 0.1_R8) then
         zeps = 0.1_R8
         ierr = 1
      end if
      if (zeps .gt. 0.3_R8) then
         zeps = 0.3_R8
         ierr = 1
      end if
      if (zrln .lt. 1.0E-6_R8) then
         zrln = 1.0E-6_R8
         ierr = 1
      end if
      if (zrln .gt. 6.0_R8) then
         zrln = 6.0_R8
         ierr = 1
      end if
!
!
!     Initial calculations (check tau, nu, zeff for validity)
!
      if ( switches(2) .eq. 0 ) then
         vthermi = sqrt(tikev*1.6022E-16_R8/3.3436E-27_R8)
         omegaci = 1.6022E-19_R8* btesla / 3.3436E-27_R8
         rhoi = vthermi/omegaci
         if ( (tikev .lt. 1.0E-4_R8) .or. (tekev .lt. 1.0E-4_R8) ) then
            znu = 10.0_R8
            ierr = 1
         else
            znu = 2.1_R8*rmajor*ne19/(tekev**1.5_R8* tikev**0.5_R8)
            if ((znu .lt. 0.5_R8) .and. (switches(3) .eq. 0)) then
               znu = 0.5_R8
               ierr = 1
            end if
            if ( (znu .gt. 10.0_R8) .and. (switches(5) .eq. 0) ) then
               znu = 10.0_R8
               ierr = 1
            end if
         end if
         if ( tekev .lt. 1.0E-4_R8) then
            ztau = 4.0_R8
            ierr = 1
         else
            ztau = tikev/tekev
            if (ztau .lt. 0.5_R8) then
               ztau = 0.5_R8
               ierr = 1
            end if
            if (ztau .gt. 4.0_R8) then
               ztau = 4.0_R8
               ierr = 1
            end if
         end if
         gvti = vthermi
         grhoi = rhoi
         gtau = ztau
         gnu = znu
      end if
!
      if ( switches(2) .eq. 1 ) then
         vthermi = gvti
         rhoi = grhoi
         if (rhoi .lt. 1.0E-10_R8) then
            rhoi = 1.0E-10_R8
            ierr = 1
         end if
         znu = gnu
         if (switches(4) .eq. 1) znu = 0.84_R8*gnu
         if ((znu .lt. 0.5_R8) .and. (switches(3) .eq. 0)) then
            znu = 0.5_R8
            ierr = 1
         end if
         if ( (znu .gt. 10.0_R8) .and. (switches(5) .eq. 0) ) then
            znu = 10.0_R8
            ierr = 1
         end if
         ztau = gtau
         if (ztau .lt. 0.5_R8) then
            ztau = 0.5_R8
            ierr = 1
         end if
         if (ztau .gt. 4.0_R8) then
            ztau = 4.0_R8
            ierr = 1
         end if
         tikev = 3.3436E-27_R8*gvti*gvti/(1.6022E-16_R8)
         tekev = tikev/ztau
         omegaci = gvti/grhoi
         btesla = omegaci*3.3436E-27_R8/1.6022E-19_R8
         ne19 = znu*tekev**1.5_R8*tikev**0.5_R8/(2.1_R8*rmajor)
      end if
!
      ztaub = ztau/(1.0_R8- znbne)
      zzeffstr = 1.0_R8+ 30.0_R8*zncne/(znine + 6.0_R8*zncne)
      if (zzeffstr .lt. 1.0_R8) then
         zzeffstr = 1.0_R8
         ierr = 1
      end if
      if (zzeffstr .gt. 4.0_R8) then
         zzeffstr = 4.0_R8
         ierr = 1
      end if
      zrlnstr = min(6.0_R8,zrln)
!
!
!     Are we using the 1994 IFS/PPPL model?
!     If so, go to the start of that section
!
      if ( switches(3) .eq. 1 ) go to 2000
!
!     The 1995 IFS/PPPL model
!
!     The threshold temperature gradient for the hydrogenic mode
!
      zf = 1._R8- ( 0.942_R8*zzeffstr**0.516_R8/zshat**0.671_R8)*        &  
     &          ( 2.95_R8*zeps**1.257_R8/znu**0.235_R8- 0.2126_R8)
      zg = ( 0.671_R8+ 0.570_R8*zshat - 0.189_R8*zrlnstr )**2 +          &  
     & 0.392_R8+                                                         &  
     &          0.335_R8*zrlnstr - 0.779_R8*zshat + 0.210_R8*zshat*      &  
     & zshat
      zh = 2.46_R8* (zzeffstr/2.0_R8)**0.7_R8* ztaub**0.52_R8*           &  
     &          (1.0_R8+ 2.78_R8/(zq*zq))**0.26_R8
!
      zrlt1 = zf*zg*zh
!
 
!
!     The threshold temperature gradient for the carbon mode
!
      zscriptd = max( 1.0_R8, 3.0_R8-0.6667_R8*zrlnstr )
      zscripte = 1.0_R8+ 6.0_R8*max( 0.0_R8, 2.9_R8-zzeffstr )
!
      zrlt2 = 0.75_R8*(1.0_R8+ztaub)*(1.0_R8+zshat)*zscriptd*zscripte
!
 
!
!     Calculate zkapfact (=1.0 when kappa=1.0 -- circular)
!
      zkapfact = 1.0_R8/ ( 1._R8+ ((zkappa-1._R8)*zq/3.6_R8)**2 )
!
!     Calculate the shear flow stabilization factor
!
      gamma = (vthermi/rmajor)*(0.25_R8/(ztau*(1._R8+0.5_R8*zshat*zshat)  &  
     & ))*                                                               &  
     &     ((1._R8+3._R8*max(0.0_R8,zeps-0.16667_R8)) / (1._R8+          &  
     & max(0._R8,zq-3._R8)/15._R8))*                                     &  
     &     (zrlt-zrlt1)
      zscripts = 0.0_R8
      if (gamma .gt. omegaexb) zscripts = 1.0_R8- omegaexb/gamma
!
!     Calculate zchii1
!
      zx = zrlt - zrlt1
      zscriptg1 = zx
      if ( zx .lt. 0.0_R8) zscriptg1 = 0.0_R8
      if ( zx .gt. 1.0_R8) zscriptg1 = sqrt(zx)
!
      zscriptz = min( 1.0_R8, (3.0_R8/zzeffstr)**1.8_R8)
!
      zchii1 = ( 11.8_R8*zq**1.13_R8/ztaub**1.07_R8) / (1.0_R8+ zshat**  &  
     & 0.84_R8) *                                                        &  
     &        (1.0_R8+ 6.72_R8*zeps/(zq**0.96_R8*znu**0.26_R8)) *        &  
     &        zscriptz * zscriptg1 * zkapfact * zscripts
!
!     Calculate zchii2
!
      zx = zrlt - zrlt2
      zscriptg2 = zx
      if ( zx .lt. 0.0_R8) zscriptg2 = 0.0_R8
      if ( zx .gt. 1.0_R8) zscriptg2 = sqrt(zx)
!
      zchii2 = (7.88_R8)/( ztaub**0.8_R8* (1.0_R8+ zshat) )*             &  
     &        max( 0.25_R8, zzeffstr-3.0_R8) * zscriptg2 *               &  
     &        zkapfact * zscripts
!
!     Now calculate chii
!
      zchii = max(zchii1,zchii2)
      chii = (rhoi*rhoi*vthermi/rmajor)*zchii
!
!     Get chii normalized in Cyclone units
!
      zchiicyc = zchii/zrln
!
!
!     Calculate chie correction factor
!
      zcorrect = (7.0_R8- zzeffstr)/6.0_R8
!
!     Calculate zscriptj
!
      zscriptj = max ( 2.0_R8, (1.0_R8+ 0.3333_R8*zrlnstr) )
!
!     Calculate zchie1 and zchie2
!
      zchie1 = zchii1 * 0.72_R8* max( 0.1667_R8, zeps ) * znu**0.14_R8*  &  
     &        (zq/zshat)**0.3_R8* ztau**0.4_R8*                          &  
     &        ( 1._R8+ 0.333_R8* zrlnstr)
!
      zchie2 = zchii2 * 0.263_R8* ztau * znu**0.22_R8* zscriptj
!
!     Finally, calculate chie
!
      zchie = max(zchie1,zchie2) * zcorrect
      chie = (rhoi*rhoi*vthermi/rmajor)*zchie
!
!
!     (Optional) Print out page of calculation results
!
      if ( switches(1) .eq. 1 ) then
 
         if ( switches(30) .gt. 0) then
            write (switches(30),1000) chii, ierr, chie,                  &  
     &           rhoi, vthermi,                                          &  
     &           omegaci, calcflag,                                      &  
     &           tikev, calcflag,                                        &  
     &           tekev, calcflag,                                        &  
     &           ne19, calcflag,                                         &  
     &           btesla, calcflag,                                       &  
     &           omegaexb,                                               &  
     &           zchiicyc,                                               &  
     &           zchii1, zchii2,                                         &  
     &           zchie1, zchie2, zrlt, zrlt1, zrlt2,                     &  
     &           zkapfact, zcorrect,                                     &  
     &           gamma, zscripts,                                        &  
     &           zscriptj, zscriptz,                                     &  
     &           zscriptg1, zscriptg2,                                   &  
     &           zscriptd, zscripte, zf, zg, zh,                         &  
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
         if ( switches(31) .gt. 0) then
            write (switches(31),1000) chii, ierr, chie,                  &  
     &           rhoi, vthermi,                                          &  
     &           omegaci, calcflag,                                      &  
     &           tikev, calcflag,                                        &  
     &           tekev, calcflag,                                        &  
     &           ne19, calcflag,                                         &  
     &           btesla, calcflag,                                       &  
     &           omegaexb,                                               &  
     &           zchiicyc,                                               &  
     &           zchii1, zchii2,                                         &  
     &           zchie1, zchie2, zrlt, zrlt1, zrlt2,                     &  
     &           zkapfact, zcorrect,                                     &  
     &           gamma, zscripts,                                        &  
     &           zscriptj, zscriptz,                                     &  
     &           zscriptg1, zscriptg2,                                   &  
     &           zscriptd, zscripte, zf, zg, zh,                         &  
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
         if ( switches(32) .gt. 0) then
            write (switches(32),1000) chii, ierr, chie,                  &  
     &           rhoi, vthermi,                                          &  
     &           omegaci, calcflag,                                      &  
     &           tikev, calcflag,                                        &  
     &           tekev, calcflag,                                        &  
     &           ne19, calcflag,                                         &  
     &           btesla, calcflag,                                       &  
     &           omegaexb,                                               &  
     &           zchiicyc,                                               &  
     &           zchii1, zchii2,                                         &  
     &           zchie1, zchie2, zrlt, zrlt1, zrlt2,                     &  
     &           zkapfact, zcorrect,                                     &  
     &           gamma, zscripts,                                        &  
     &           zscriptj, zscriptz,                                     &  
     &           zscriptg1, zscriptg2,                                   &  
     &           zscriptd, zscripte, zf, zg, zh,                         &  
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
 1000    format (/,' Results from the 1995 IFS/PPPL Model',/,            &  
     &       /,' Dimensional:',                                          &  
     &       /,' chi-i =',1pe21.14,' (m^2/s)','   error code =',i2,      &  
     &       /,' chi-e =',1pe21.14,' (m^2/s)',/,                         &  
     &       /,' rho-i =',1pe21.14,' meters',                            &  
     &       /,' v-therm-i =',1pe21.14,' m/s',                           &  
     &       /,' omega-ci =',1pe21.14,' rad/s  ',a35,                    &  
     &       /,' T-i =',1pe21.14,' keV         ',a35,                    &  
     &       /,' T-e =',1pe21.14,' keV         ',a35,                    &  
     &       /,' n-e =',1pe21.14,' 10^19 m^-3  ',a35,                    &  
     &       /,' B-tor =',1pe21.14,' Tesla     ',a35,                    &  
     &       /,' omega-ExB =',1pe21.14,' rad/s  ',/,                     &  
     &       /,' Dimensionless:',                                        &  
     &       /,' chi-i =',1pe21.14,' in Cyclone units',                  &  
     &       /,' zchi-i1 =',1pe21.14,' zchi-i2 =',1pe21.14,              &  
     &       /,' zchi-e1 =',1pe21.14,' zchi-e2 =',1pe21.14,/,            &  
     &       /,' R/L_T =',1pe21.14,                                      &  
     &       /,' R/L_Tcrit1 =',1pe21.14,' R/L_Tcrit2 =',1pe21.14,/,      &  
     &       /,'  zkapfact =',1pe21.14,'  zcorrect =',1pe21.14,          &  
     &       /,' gamma s^-1=',1pe21.14,'  zscripts =',1pe21.14,          &  
     &       /,'  zscriptj =',1pe21.14,'  zscriptz =',1pe21.14,          &  
     &       /,' zscriptg1 =',1pe21.14,' zscriptg2 =',1pe21.14,          &  
     &       /,'  zscriptd =',1pe21.14,'  zscripte =',1pe21.14,/,        &  
     &       /,'       zf =',1pe21.14,'      zg =',1pe21.14,             &  
     &       /,'       zh =',1pe21.14,'    ztau =',1pe21.14,             &  
     &       /,'    ztaub =',1pe21.14,'     znu =',1pe21.14,             &  
     &       /,' zzeffstr =',1pe21.14,' zrlnstr =',1pe21.14,             &  
     &       /)
 
      end if
!
      go to 3000
!
!
!     The 1994 IFS/PPPL model
!
 2000 continue
!
      zc1 = 1.26_R8+ abs(zzeffstr-3._R8)*                                &  
     &             (-0.27_R8+0.075_R8*zzeffstr-0.044_R8*zzeffstr*        &  
     & zzeffstr)
      zc1 = max(0.57_R8,zc1)
!
      zc2 = 2.0_R8*(3.5_R8-zzeffstr) + 0.1_R8/ztau
      if ( zzeffstr .lt. 3.0_R8) zc2 = 1.0_R8+ 0.1_R8/ztau
      if ( zzeffstr .gt. 3.5_R8) zc2 = (0.1_R8+ 0.2_R8*(zzeffstr-3.5_R8)  &  
     & )/ztau
!
      zc3 = 0.61_R8- 0.27_R8/( 1.0_R8+ exp( 8.0_R8*(3.3_R8-zzeffstr) ) )    
!
!
!     Calculate zscripty
!
      zscripty = abs(0.1976_R8-0.4550_R8*zshat+0.1616_R8*zrln)**         &  
     & 0.769_R8                                                          &  
     &           + 0.7813_R8+ 0.2762_R8*zshat + 0.3967_R8*zshat*zshat
!
!     Calculate the critical temperature gradient
!
      zrlt1 = 2.778_R8* zc1 * ztaub**zc3 * sqrt(0.5_R8+1._R8/zq) *       &  
     &         abs(1._R8- 0.85_R8*zeps/zshat**0.25_R8) * zscripty
!
!
!     Calculate zscriptg and chii
!
      zx = zrlt - zrlt1
      zscriptg = zx
      if ( zx .lt. 0.0_R8) zscriptg = 0.0_R8
      if ( zx .gt. 1.0_R8) zscriptg = sqrt(zx)
!
      zchii1 = 14._R8*zq*zc2/(ztaub*(2._R8+zshat))                       &  
     &         * (1._R8+ zeps/3._R8) * zscriptg
      chii = zchii1 * rhoi*rhoi*vthermi/rmajor
!
!     Get chii normalized in Cyclone units
!
      zchiicyc = zchii1/zrln
!
!     Zero out zchii2, zchie2 and zrlt2 (they are not used in the 1994 model)
!
      zchii2 = 0.0_R8
      zchie2 = 0.0_R8
      zrlt2 = 0.0_R8
!
!
!     Calculate chie and zchie
!
      zchie1 = zchii1 * 0.27_R8* ztau**0.7_R8
      chie = chii * 0.27_R8* ztau**0.7_R8
!
!
!     (Optional) Print out page of calculation results
!
      if ( switches(1) .eq. 1 ) then
 
         if ( switches(30) .gt. 0) then
            write (switches(30),5000) chii, ierr, chie,                  &  
     &           rhoi, vthermi,                                          &  
     &           omegaci, calcflag,                                      &  
     &           tikev, calcflag,                                        &  
     &           tekev, calcflag,                                        &  
     &           ne19, calcflag,                                         &  
     &           btesla, calcflag,                                       &  
     &           zchiicyc, zrlt, zrlt1,                                  &  
     &           zscriptg, zc1, zc2, zc3,                                &  
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
         if ( switches(31) .gt. 0) then
            write (switches(31),5000) chii, ierr, chie,                  &  
     &           rhoi, vthermi,                                          &  
     &           omegaci, calcflag,                                      &  
     &           tikev, calcflag,                                        &  
     &           tekev, calcflag,                                        &  
     &           ne19, calcflag,                                         &  
     &           btesla, calcflag,                                       &  
     &           zchiicyc, zrlt, zrlt1,                                  &  
     &           zscriptg, zc1, zc2, zc3,                                &  
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
         if ( switches(32) .gt. 0) then
            write (switches(32),5000) chii, ierr, chie,                  &  
     &           rhoi, vthermi,                                          &  
     &           omegaci, calcflag,                                      &  
     &           tikev, calcflag,                                        &  
     &           tekev, calcflag,                                        &  
     &           ne19, calcflag,                                         &  
     &           btesla, calcflag,                                       &  
     &           zchiicyc, zrlt, zrlt1,                                  &  
     &           zscriptg, zc1, zc2, zc3,                                &  
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
 5000    format (/,' Results from the 1994 IFS/PPPL Model',/,            &  
     &       /,' Dimensional:',                                          &  
     &       /,' chi-i =',1pe21.14,' (m^2/s)','   error code=',i2,       &  
     &       /,' chi-e =',1pe21.14,' (m^2/s)',/,                         &  
     &       /,' rho-i =',1pe21.14,' meters',                            &  
     &       /,' v-therm-i =',1pe21.14,' m/s',                           &  
     &       /,' omega-ci =',1pe21.14,' rad/s  ',a35,                    &  
     &       /,' T-i =',1pe21.14,' keV         ',a35,                    &  
     &       /,' T-e =',1pe21.14,' keV         ',a35,                    &  
     &       /,' n-e =',1pe21.14,' 10^19 m^-3  ',a35,                    &  
     &       /,' B-tor =',1pe21.14,' Tesla     ',a35,/,                  &  
     &       /,' Dimensionless:',                                        &  
     &       /,' chi-i =',1pe21.14,' in Cyclone units',                  &  
     &       /,' R/L_T =',1pe21.14,'  R/L_Tcrit1 =',1pe21.14,/,          &  
     &       /,'  zscriptg =',1pe21.14,/,                                &  
     &       /,' zc1 =',1pe21.14,' zc2 =',1pe21.14,                      &  
     &       /,' zc3 =',1pe21.14,/,                                      &  
     &       /,'    ztau =',1pe21.14,'    ztaub =',1pe21.14,             &  
     &       /,'     znu =',1pe21.14,' zzeffstr =',1pe21.14,             &  
     &       /,' zrlnstr =',1pe21.14,                                    &  
     &       /)
 
      end if
!
!
 3000 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
