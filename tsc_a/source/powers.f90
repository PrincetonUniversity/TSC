      subroutine powers(formx,formne,formte,formj)
      USE CB_ALL
!
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 formne,formte,formj,formx
      REAL*8 coulom,zj0,pscale,palfa,gsyn,rtlam,zphi
!============
      dimension formx(*),formne(*),formte(*),                            &  
     &   formj(*)
!============      
!     REAL*8 kappa
!     common/cb1/rminor,rmajor,kappa,area,volm,curent,btor,qcyl
!     common/cb2/cboltz,pi
!     common/cb3/alpj,alpt,alpn,tkv0,tkvav,dene0,deneav,deneln
!     common/cb4/ndens,ntemp,nform
!     common/cb5/zeff,zimp,aibar,fdeu,fdt,rnzne,rnine,resnc,refl
!     common/cb6/ctau
!     common/cb7/wtot,fpal,fpoh,fpbr,fpsync
!     common/cb8/zpin,zpaux
!     common/cb9/tscoh,tscal,tscbr,tscsyn
!     common/cb10/coh,cal,cbr,csyn,cwtot
!..plagiarised version of Darren Stotler's subrtoutine plspow
!..computes ohmic heating, alpha heating, radiation losses and
!..total plasma thermal energy for given density and temperature.
!
!
!..dene0     = central density (in m**-3)
!..deneav    = vol averaged density (in m**-3)
!..gn1       = <n>/n(0)
!..tkv0      = central temperature (in kev)
!..tkvav     = vol averaged temperature (in kev)
!..alpt      = T(0)/<T> - 1
!..alpn      = n(0)/<n> - 1
!..alpj      = J(0)/<J> - 1
!..coulom    = coulomb log in mks-kev
!..zj0       = central current density
!..curent    = total plasma current (in amps)
!..area      = cross sectional area of plasma (in m**2)
!..gj1       = <J>/J0
!..fpoh      = ohmic heating power (in watts)
!..resnc     = neoclassical resistivity correction
!..volm      = volume of plasma (in m**3)
!..pscale    = scale factor for ohmic heating power
!..cboltz    = boltzmann constant (in J/kev)
!..fdeu      = D density/ D-T density
!..fdt       = D-T density/ electron density
!..nform    = number of radial zones for profile form arrays
!..formte    = T(r)/T(0) array
!..formx     = r/a       array
!..formne    = ne(r)/ne(0) array
!..fpbr      = bremsstrahlung radiation (in Watts)
!..fpsync    = synchrotron radiation losses (in watts)
!..rminor    = minor radius of plasma (in m)
!..rmajor    = major radius of plasma (in m)
!..btor      = toroidal magnetic field (in T)
      coulom=37.8_R8-log(sqrt(deneav)/tkvav)
!..central current density
      zj0=curent/area*(alpj+1)
!..Ohmic heating power (W) (see Uckan, fusion Technology 14,299,(1988))
      fpoh=1.65E-9_R8*zeff*coulom*resnc*tkv0**(-1.5_R8)*zj0**2*volm/     &  
     &   (1.+2.*alpj-1.5*alpt)
      fpoh=coh*fpoh
!..scale factor for alpha power
      pscale=3.5E3_R8*cboltz*fdeu*(1-fdeu)*(fdt*dene0)**2*volm
!..alpha power (W)
      fpal=palfa(pscale,formte,formx,formne)
      fpal=cal*fpal
      fpbr=1.0E6_R8*0.0168_R8*(dene0/1.E20_R8)**2                        &  
     &   *sqrt(tkv0/10.)*zeff*volm/(1.+2.*alpn+0.5*alpt)
      fpbr=cbr*fpbr
!..synchrotron radiation (in W) using Trubnikov formula from
!..Houlberg CIT report P-870615-ppl-01; replace n and T with volume
!..averages.
      gsyn=5.198E-3_R8*tkvav**1.5_R8*sqrt(1._R8+22.61_R8*rminor          &  
     &   /(rmajor*sqrt(tkvav)))
      rtlam=7.78E-9_R8*sqrt(deneav*rminor/btor)
      zphi=(gsyn/rtlam)*sqrt(1-refl)
      fpsync=6.214E-17_R8*deneav*tkvav*btor**2*zphi*volm
      fpsync=csyn*fpsync
!..total thermal energy
      wtot=1.5_R8*deneav*(1._R8+rnine)*cboltz*tkvav*volm
      wtot=cwtot*wtot
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
