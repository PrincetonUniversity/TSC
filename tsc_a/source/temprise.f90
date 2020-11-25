      subroutine temprise
!......3.80 auxdef
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.....     TEMPIRSE.FOR                                 02/06/87
!.....
!.....     R.D.Pillsbury, Jr.
!.....     Plasma Fusion Center
!.....     MIT  (617)-253-8605
!.....
!.....     Interface routine with TSC to take coil group currents
!.....     and calculate the adiabatic temperature rise in each
!.....     "coil". Quotes are used since TSC may use multiple
!.....     filaments to model a coil. However, the easiest approach
!.....     is to treat each filament as an independent coil with
!.....     input builds, copper and steel fractions, and initial
!.....     temperatures for each.
!.....
!.....     This routine is called with parameters:
!.....
!.....     ntime     - time counter - assumes at ntime=1 , time=0
!.....     ncoil     - actual number of coils
!.....     dt        - time step (s) - i.e. = t_present - t_last_call
!.....     dxcoil(n) - radial build of nth coil (m)
!.....     dzcoil(n) - radial build of nth coil (m)
!.....     fcu(n)    - copper fraction of nth coil
!.....     fss(n)    - Inconel fraction of nth coil
!.....     ccoil(n)  - current per turn of nth coil at the present time (A)
!.....     aturnsc(n)- number of turns in nth coil
!.....     tempc(n)   - temperature of nth coil (K)
!.....                 on input at previous time
!.....                 on output at present time
!.....
!.....     This coding minimizes the necessity of matching common blocks
!.....     and cliches.
!.....
!.....     Routine preforms an adiabatic temperature rise calculation
!.....     for uniform current density coils.  Individual coils may have
!.....     different copper fractions. Routine preforms a mixture of the
!.....     heat capacities of copper and steel for internally re-inforced
!.....     coils. NOTE: This routine assumes any steel is 100% effective
!.....     and, hence, the conduction of heat into the steel is perfect.
!.....     To approximate diffusion effects into the steel, a reduced fss
!.....     can be used.
!.....
!.....     Routine calculates an incremental temperature rise with thermal
!.....     properties evaluated at the previous time step according to:
!.....
!.....  dT = fcu*(rho(T)*jcu(t)**2)/(fcu*gamcu*Cpcu(T)+fss*gamss*Cpss(T))*dt
!.....
!.....     where
!.....
!.....     temp   - temperature at previous time step (K).
!.....     fcu    - conductor fraction (vol of conductor)/(total vol)
!.....     fss    - steel fraction (vol of steel)/(total vol)
!.....     rho(T) - electrical resisitivity of conductor at temp. (Ohm-m)
!.....     jcu(t) - conductor current density at present time, t. (A/m**2)
!.....     gamcu  - mass density of conductor (kg/m**3).
!.....     gamss  - mass density of steel (kg/m**3).
!.....     Cpcu(T)- specific heat of conductor at temp (J/kg/K)
!.....     Cpss(T)- specific heat of steel at temp (J/kg/K)
!.....     dt     - time step (s).
!.....
!.....
!.....  MODIFICATIONS:
!.....   05/05//89  - Modified to allow different conducting material in
!.....                different coils. In order to have no additional
!.....                input fields, the parameter fcu is used.  If
!.....                0.000 < fcu < 0.999   - OFHC Copper
!.....                1.000 < fcu < 1.999   - AL25 (Glidcop)
!.....                2.000 < fcu < 2.999   - Berylium Copper
!.....
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,ig,iabs,itype,igg
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 gamss,dtimes,xi,ffcu,rho,gamcu,cpcu,cpss,areacu,xjcu
      REAL*8 heatcap,deltat
!============
      data gamss/7800._R8/
!.....
!.....     NB - the copper area for each coil should be checked at
!.....     input to insure a non-zero area -- otherwise jcu blows up.
!.....
!
      dtimes=times-tlast
!
!
!ccccccinput initial temperatures, cu and ss fractions, radial+vert
!ccccccbuilds of coils
!
!
!
      if(ifrst(11).le.1) then
        ifrst(11)=2
        numpf = 0
        do 7 n=1,ncoil-nwire
          ig = iabs(igroupc(n))
          if(ig.gt.numpf) numpf = ig
          if(dxcoil(n).le.0) ineg=30
          if(dzcoil(n).le.0) ineg=30
          if(fcu(n) .le. 0) ineg=30
    7   continue
      else
        do 100 n=1,ncoil-nwire
          xi    = ccoil(n)*udsi
          itype = int(fcu(n))+1
          ffcu  = fcu(n)-itype+1
          call getrho(itype,tempc(n),rho,gamcu)
          call cucpvst(tempc(n),cpcu)
          call sscpvst(tempc(n),cpss)
          areacu     = ffcu*dxcoil(n)*dzcoil(n)
          xjcu       = xi/areacu
          heatcap    = ffcu*gamcu*cpcu+fss(n)*gamss*cpss
          deltat     = ffcu*(rho*xjcu**2/heatcap)*dtimes
          tempc(n)   = tempc(n)+deltat
          rscoils(n) = tpi*xcoil(n)*rho/areacu
          rscoil(n)  = rscoils(n)*usdr
  100   continue
      end if
      do 11 igg=1,numpf
        do 12 n=1,ncoil-nwire
          nn = n
          ig = iabs(igroupc(n))
          if(ig.eq.igg) go to 13
   12   continue
        templ(igg) = 0._R8
        go to 11
   13   templ(igg) = tempc(nn)
   11 continue
      tlast=times
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
