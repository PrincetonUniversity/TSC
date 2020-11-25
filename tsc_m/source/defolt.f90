      subroutine defolt
!
!***********************************************************************
!                                                                      *
!.....sets default values for all input variables                      *
!                                                                      *
!***********************************************************************
!
      USE CLINAM
      USE RUNAWAY
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!...........................................................
!
!....... scalar integer
!
!...........................................................
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER l,j,ll, i
!============
      ineg   = 0
      ires   = 0
      isurf  = 0
      igone  = 0
      irayt = 0
      iplot = 0
      nskipsfi = 20
      iwall  = 0
      nx     = 46
      nz     = 16
      ndiv   = 1
      nc0    = 1
      nthe   = 100
      istart = 1
      ntpts = 2
      ncnt = 0
      itrmod = 2
      nplate = 0
      ngrmax = 0
      nmult = 0
!
!.....time programming points
      tpros(1) = 0._R8
      tpros(2) = 1.E6_R8
!...........................................................
!
!....... scalar floating point
!
!...........................................................
      alx    = 3.0_R8
      alz    = 1.0_R8
      dtmin  = .001_R8
      dtmax  = 1.0_R8
      dtfac  = 0.5_R8
      xlim   = 0._R8
      zlim   = 100._R8
      th     = 0.5_R8
      ffac = 1.0_R8
      reboun = 0._R8
      erate  = 0._R8
      small  = 1.E-8_R8
      amu    = 0._R8
      eta    = 0._R8
      zeff   = 1.0_R8
      amgas = 1.0_R8
      zgas = 1.0_R8
      tfmult = 16.0_R8
!
!.....dominant impurity
      zimp = 6.0_R8
!
!.....neutral beam form factor
      abeam = 0.25_R8
      dbeam = 0.1_R8
      nebeam = 1.0_R8
!....current drive parameters
      ebeamkev = 80._R8
      ambeam = 1._R8
      fracpar = 0._R8
      ibootst = 0
      iicrh = 0
      ifwcd = 0
      iecrh = 0
      ieccd = 0
      idefnpe = 0
!
!.....density exponents
      alphar = 0.5_R8
      betar = 2.0_R8
!
!
!.....vacuum temperature
      tevv = 1.0_R8
!
!.....reference density
      udsd = 1.0E19_R8
!
!.....switches for alpha calculation and ballooning calculation
      ialpha = 0
      ibalsw = 0
!.....resistivity enhanced for q<qsaw
      qsaw = 1.0_R8
!...........................................................
!
!....... acoef array
!
!............................................................
      do 100 l=1,4999+pncoil
  100 acoef(l) = 0._R8
!
!......ratio between initial electron/ion pressure
      acoef(2) = 0.5_R8
!......relaxation factors for pressure and g function(isurf=1 only)
      acoef(4) = 100._R8
      acoef(5) = 0.0_R8
!......parallel diffusion multiplier (isurf=0 only)
      acoef(6) = 10._R8
!.....mix between dufort frankel and for time cent space in flux diffusion
      acoef(7) = .85_R8
!......implicit parameter for surface averaged variables time advancement
      acoef(8) = 1._R8
!......numerical viscosity coefficient
      acoef(9) = 40._R8
!
!
!......ratio of incompressible to compressible viscosity
      acoef(10) = 0.50_R8
!.....constant in def of rebounn
      acoef(11) = 0.5_R8
!
!......initialization current for lrswtch .gt.0
      acoef(12) = 1._R8
!
!.....boundary flux switch
      acoef(13) = 4._R8
!
!.....max voltage for loop voltage plot
      acoef(14) = 10._R8
!.....min and max for loop voltage
      acoef(15) = -100._R8
      acoef(16) = +100._R8
!
!......switch for cubic time point interpolation
      acoef(17) = 0._R8
!.....multiplier of psidot boundary term
      acoef(18) = 0.0_R8
!.....switch for ucor
      acoef(20) = 0._R8
!.....maximum amach and ekin allowable before error switch set
      acoef(21) = 1._R8
      acoef(22) = 100._R8
!.....parameters used in equilibrium functions
      acoef(25) = 1.0_R8
      acoef(26) = 1.0_R8
      acoef(27) = 1.0_R8
!.....convergence criterion for initial equilibrium
      acoef(23) = 1.E-7_R8
      acoef(24) = 1.E-6_R8
!.....bypass initial filiment growth rate calculation if non-zero
      acoef(28) = 1._R8
!.....stop calculation if time exceeds this value
      acoef(29) = 1000._R8
!......anomolous transport coefficients chie,chii,d
      acoef(35) = 1._R8
      acoef(37) = 1._R8
      acoef(39) = .2_R8
!.....resgap coeffiecient of coil gap resistance
      acoef(41) = 0.5_R8
!.....xprof for profile plot
      acoef(42) = 0.0_R8
      acoef(43) = 0.0_R8
!.....no of zones to search over for x-point
      acoef(45) = 2._R8
!.....max for taue plot(sec)
      acoef(46) = 2._R8
!.....max power for burn control simulation
      acoef(47) = 1.E12_R8
!
!.....no. of plasma flux contours
      acoef(48) = 20
!
!......no of vacuum flux contours
      acoef(49) = 20
!
!.....relaxation factors for initial equilibrium
      acoef(50) = 0.5_R8
      acoef(51) = 0.5_R8
!.....coefficent in plasma current burn control
      acoef(54) = 0._R8
!
!.....reflectivity of walls for cyclotron radiation calculation
      acoef(55) = 0.9_R8
!.....beginning and ending times for hyper
      acoef(57) = 1.E6_R8
      acoef(58) = 0.00_R8
!
!.....convergence criteria for hyper iteration
      acoef(59) = 1.E-6_R8
      acoef(60) = 4000._R8
!
!.....if non-zero, zmag time history plotted even for isym=1
      acoef(61) = 0._R8
!
!.....ratio of toroidal to compressible viscocity
      acoef(62) = 1._R8
      acoef(63) = .667_R8
!
!......coefficients needed for hyperresistivity term
      acoef(64)=0.1_R8
      acoef(65)=0.3_R8
      acoef(66)=1.0_R8
      acoef(67)=1.2_R8
      acoef(68)=1.0_R8
      acoef(69) = .010_R8
!
!
!.....relaxation coef for resistivity when lrswtch.ne.0
      acoef(70) = .0001_R8
!
!
!.....max temp for resistivity calculation
      acoef(71) = 1.E6_R8
!.....fraction of ion neoclassical for electrons
      acoef(73) = 1._R8
!......special limiter adjustment switch
      acoef(74) = 0._R8
!.....npert no of cycles resistivity is enhanced to let perturbation in
      acoef(75) = 0.0_R8
!.....drag terms
      acoef(90) = 0.2_R8
      acoef(91) = 0.2_R8
      acoef(92) = 0.2_R8
!
!.....confinement time for helium ash
      acoef(93) = 1._R8
!
!.....disruption time and qsaw after disruption
      acoef(95) = 1.E6_R8
      acoef(96) = 2.0_R8
!.....halo ratio
      acoef(97) = 0.0_R8
!.....halo temperature
      acoef(98) = 1.0_R8
!.....multiplier of bootstrap current
      acoef(103) = 1.0_R8
!.....multiplier of ajavlh
      acoef(106) = 1._R8
!
!.....anom aux heating transport coeff
      acoef(107) = 1._R8
      acoef(108) = 1.E6_R8
!.....concentration of deuterium for alpha heating calculation
      acoef(113) = 0.49_R8
      acoef(114) = 0._R8
!
!  resistivity flattening inside qsaw radius
      acoef(120)= 0.5_R8
!......coefficients for anomolous transport model 2
      acoef(121) = .08_R8
      acoef(122) = 0.42_R8
      acoef(123) = .5_R8
      acoef(124) = 2.0_R8
      acoef(125) = 0.0_R8
      acoef(126) = 2.0_R8
      acoef(127) = 0._R8
      acoef(128) = 0.5_R8
      acoef(129) = 2.0_R8
      acoef(700) = 50._R8
      acoef(701) = 10._R8
!
!.....these coefficients needed for ilhcd=2
      nslhrt = acoef(700)
      nslhpc = acoef(701)
!
!.....these coefficients needed for iffac=1
      acoef(801) = .005_R8
      acoef(802) = 0.9_R8
      acoef(803) = 1.1_R8
      acoef(804) = 1000._R8
      acoef(805) = 1.0_R8
      acoef(806) = 0.5_R8
!.....multiplies ablation rate
      acoef(809) = 1._R8
!.....multiplies eta
      acoef(810) = 1.0_R8
!.....to limit plasma by q-value
      acoef(811) = 0._R8
!.....these coefficients needed for subrotuine iterate
      acoef(821) = 1.75_R8
      acoef(822) = 1.4_R8
      acoef(823) = 1.E-8_R8
      acoef(824) = 4000._R8
      acoef(831) = 1.85_R8
      acoef(832) = 1.4_R8
      acoef(833) = 1.E-8_R8
      acoef(834) = 2000._R8
      acoef(841) = 1.62_R8
      acoef(842) = 1.38_R8
      acoef(843) = 1.E-8_R8
      acoef(844) = 4000._R8
!
      acoef(852:861) = 0.0_R8
      acoef(851) = 0.0_R8
!.....multiplier of Bremsstrahlung radiation
      acoef(877) = 1._R8
!
!     edge electron temp for transport calculation if .gt.0
      acoef(880)=0._R8
!     edge density fraction for transport calculation
      acoef(881)=0.1_R8
!     edge ratio of total pressure to electron pressure
      acoef(882)=2._R8
!...heat conduction multiplier
      acoef(890) = 1._R8
!...heat conduction denominator
      acoef(891) = 100._R8
!
!   coefficient for Coppi model  
	acoef(893) = 0.0_R8
!......coefficient of gamma in heat flux equation
      acoef(897) = 3._R8/2._R8
!
!
!.....coefficients for h-mode transport
      acoef(3003) = 0._R8
      acoef(3011) = .95_R8
!
!.....relaxation in time and space for glf23 and mmm95
      acoef(3010) = 0.1_R8
      acoef(3009) = .01_R8
!
!.....settings and relaxation for global confinement constraint on chi
      acoef(3007) = 0._R8
      acoef(3008) = 1._R8
!
!.....coefficients that appear in zeffa
      acoef(3012) = 0._R8
      acoef(3013) = 0._R8
      acoef(3014) = 0._R8
!     
      acoef(3015) = 0.0_R8
!
!......for the hyper term in advsfa
      acoef(3101) = 0.       !  to turn on hyper, set to 0.008 to 0.020
      acoef(3102) = 5000.     ! pedestal temp in eV when hyper is used
!
!     acoef used in  the CDBM model
      acoef(3200) = 1.0      ! icdbm = select the model 
      acoef(3201) = 12.0     ! multiplicative factor for thermal diffusivity 
      acoef(3202) = 1.0      ! parameter in ExB shear term
      acoef(3203) = 1.0      ! tor. flux outside which Coppi-Tang is used  
      acoef(3204) = 0.75     ! tor. flux inside which CDBM is applied
      acoef(3205) = 1.0      ! if >0 then chiped is replaced by 
!				acoef(3205)*q(0)/qmin - ised for itrmod=17
      acoef(3206) = 0.0      ! if >0 then write output quantities from CDBM
      acoef(3207) = 0.0      ! if >0 then the match between CDBM and Coppi-Tang 
!				is done with TANH function instead of linear

      acoef(3208) = 0.0      ! if =1 then uses pedestal values from EPED1  
                             !    =2 uses Sugihara - PPCF 45 L55-L62 (2003)
                             !    =3 uses Maget - NF 53 093011 (2013)  
      acoef(3209) = 5.0      ! pedestal width (cm) for acoef(3208)=3.0      
      acoef(4948) = 0.       ! switch to use thermal diffusivity profile 
!
!     acoef(4950)=time to run transp
!     acoef(4951)=integration interval in transp
      acoef(4950) = 100000000.0
      acoef(4951) = 0.0
!     grid for plasma state
      acoef(4953) = 0.0
!     zimp
      acoef(4954) = 0.0
!     acoef(4955)=hydrogen fraction used for outputing to the plasma state
!     note: acoef(113)= deuterium fraction in DT mix
      acoef(4955) = 0.0
      acoef(4956) = 0.0
!     newton steps for glf
      acoef(4957) = 1.0
      acoef(4958) =-1.0
!     switches for transp
      acoef(4959) = 0.0      !swim
      acoef(4960) = 0.0      !transp
      acoef(4961) = 0.0      !write ps
      acoef(4962) = 0.0      !use tsc beam
      acoef(4963) = 0.0      !use tsc fi
      acoef(4964) = 0.0      !use tsc fw
      acoef(4965) = 0.0      !use tsc lh
      acoef(4966) = 0.0      !use tsc ech
      acoef(4967) = 0.0      !use external beam data
      acoef(4968) = 0.0      !use external lh data
      acoef(4969) = 0.0      !use external fw data
      acoef(4970) = 0.0      !use external ne data
!     fw heating partition
      acoef(4971) = 1.0        ! savefw
      acoef(4972) = 0.0      !use external rotation data
!     feedback on stored energy in auxheat
      acoef(4973) = 0.0        ! beam
      acoef(4974) = 0.0        ! rf
!     internal source models
      acoef(4975) = 0.0        ! nb
      acoef(4976) = 0.0        ! rf
!     density profile
      acoef(4977) = 0.0  
!     
      acoef(4978) = 0.0        !use radiation data      
!
      acoef(4979) = 0.0        !use external EC data
!     some glf parameters
      acoef(4981) = 0.0        ! itport
      acoef(4982) = 1.0
      acoef(4983) = 1.0
      acoef(4984) = 0.0
      acoef(4985) = 0.0
      acoef(4986) = 0.0        ! angrotp
      acoef(4987) = 0.0        ! alpha
      acoef(4988) = 0.0        ! irotstab
      
      acoef(4990) = 0.0		! if set to 1, uses total power from user files (rf, beam) to normalize profiles 
      				! instead of useing the input power from the inputa file 

      do 200 j=1,ppsi
  200 etafac(j) = 1._R8
!
      do ll=1,ptpts
      dzerv(ll) = 1._R8
      vloopv(ll) = 0._R8
      xmagz(ll) = 0._R8
      zmagz(ll) = 0._R8
      rzerv(ll) = 0._R8
      azerv(ll) = 0._R8
      ezerv(ll) = 0._R8
      zmagz(ll) = 0._R8
      rnorm(ll) = 0._R8
      gzerov(ll) = 0._R8
      tevv0(ll) = 0._R8
      ffac0(ll) = 0._R8
      zeffv(ll) = 0._R8
      alpharv(ll) = 0._R8
      betarv(ll) = 0._R8
      frcparv(ll) = 0._R8
      thalov(ll) = 0._R8
      whalov(ll) = 0._R8
      heactv(ll) = 0._R8
      tpedv(ll) = 0._R8
      qaddv(ll) = 0._R8
      fhmodeiv(ll) = 0._R8
      pwidthcv(ll) = 0._R8
      chipedv(ll) = 0._R8
      do i=1,8
        fraciv(i,ll) = 0._R8
      enddo
      enddo
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
