      subroutine popcon(a1,a2,a3,a4,a5,a6,a7p,a8p,a9,a10,a11,a12,a13,    &  
     &   a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24)

      USE CB_ALL
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ncontr,kcycle,nforms,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 a2,a3,a4,a5,a6,a7p,a8p,a9,a10,a11,a12,a13,a1
      REAL*8 formx,formne,formte,formj,cpaux
      REAL*8 tkvary,denary,conval,betalm,denlim,denmin,denmax,temin
      REAL*8 temax,delta,avdens,tsct0,avtemp,a14,a15,a16,a17,a18
      REAL*8 times,a19,a20,a21,a22,a23,tscwtot,a24,dxform,dline
      REAL*8 fpaux
      character*80 s100
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
!
      dimension formx(101),formne(101),formte(101),formj(101)
      dimension cpaux(50,50),tkvary(50),denary(50),conval(11)
      dimension betalm(50),denlim(50)
!
      data denmin/0.0_R8/,denmax/10.0_R8/,temin/0.0_R8/,temax/20.0_R8/
      data conval/-50._R8,-40._R8,-30._R8,-20._R8,-10._R8,0._R8,10._R8,  &  
     & 20._R8,30._R8,40._R8,50._R8/

!.......disabled 1/18/10
      return
      pi = 3.1415926535_R8
      cboltz = 1.6021E-16_R8
      ndens = 50
      ntemp = 50
      nform = 101
      ncontr = 11
      resnc = 2.5_R8
!
      ctau=a1
      alpn=a2
      alpt=a3
      alpj=1.5_R8*alpt
      curent=a5
      btor=a6
      rmajor=a7p
      rminor=a8p
      kappa=a9
      delta=a10
      avdens=a11
      tsct0=a12
      avtemp=tsct0*(1+alpn)/(1+alpn+alpt)
      zeff=a13
      aibar=a14
      zimp=a15
      fdeu=a16
      refl=a17
      kcycle=a18
      times=a19
      tscoh=a20
      tscal=a21
      tscbr=a22
      tscsyn=a23
      tscwtot=a24
!
      area=kappa*pi*rminor**2*(1-delta**2/8._R8)
      volm=2*pi*rmajor*area
      qcyl=2.5E6_R8*rminor**2*btor*(1._R8+kappa**2*                      &  
     &   (1.+2.*delta**2-1.2*abs(delta)**3))/(curent*rmajor)
!..fdt = D-T density/electron density
      fdt=(zimp-zeff)/(zimp-1)
!..rnzne = impurity density/electron density
      rnzne=(zeff-fdt)/zimp**2
!..rnine = ion density/electron density
      rnine=fdt+rnzne
!
!..define "form" arrays for density,temperature ans current profiles
      call forms(formx,formne,formte,formj)
!..calculate line average of density form array
      dxform=1._R8/(nform-1)
      nforms=nform
      call simson(dline,formne,dxform,nforms)
!
!..calculate popcon "powers" for tsc operating point
      deneav=avdens
      dene0=deneav*(alpn+1)
      deneln=dene0*dline
      tkvav=avtemp
      tkv0=tkvav*(1+alpn+alpt)/(1+alpn)
      coh=1._R8
      cal=1._R8
      cbr=1._R8
      csyn=1._R8
      cwtot=1._R8
      call powers(formx,formne,formte,formj)
!..evaluate "fudge factors" for popcon "powers" to agree with TSC
      coh=tscoh/fpoh
      cal=tscal/fpal
      cbr=tscbr/fpbr
      csyn=tscsyn/fpsync
      cwtot=tscwtot/wtot
!
!..loop over average density and average temperature
      do 1 i=1,ntemp
      tkvav=temin+i*(temax-temin)/ntemp
      tkv0=tkvav*(1._R8+alpn+alpt)/(1._R8+alpn)
      tkvary(i)=tkvav
!..beta limit
      betalm(i)=3._R8*curent*1.E-6_R8/(rminor*btor)/                     &  
     &   (100.*2.*4.*pi*1.e-7*cboltz*(1.+rnine)*tkvav/btor**2)/1.e20
      denlim(i)=1.5_R8*btor/rmajor
      do 1 j=1,ndens
      deneav=(denmin+j*(denmax-denmin)/ndens)*1.E20_R8
      dene0=deneav*(alpn+1)
      deneln=dene0*dline
      denary(j)=deneav/1.E20_R8
!..compute ohmic heating power, alpha heating power, bremsstrahlung
!..radiation losses, and the total plasma thermal energy for each
!..density and temperature.
      call powers(formx,formne,formte,formj)
      call pauxc(fpaux)
      cpaux(i,j)=fpaux*1.E-6_R8
1     continue
!
      call maps(0._R8,20._R8,0._R8,10._R8,.15_R8,.80_R8,.15_R8,.80_R8)
      call dders(-1)
       call setold((temax-temin)*.35_R8,-(denmax-denmin)*.10_R8,1,0,2,0)    
       write(s100,1001)
      call gtext(s100,80,0)
 1001  format("<nT>/<n> (kev)")
       call setold(-(temax-temin)*.10_R8,(denmax-denmin)*.35_R8,1,0,2,1)    
       write(s100,1002)
      call gtext(s100,80,0)
 1002  format("<n> (*e20 m-3)")
      call setold(-(temax-temin)*.15_R8,(denmax-denmin)*1.23_R8,1,0,1,0)    
      write(s100,1003) curent/1.E6_R8,btor,rmajor,rminor,kappa,delta
      call gtext(s100,80,0)
1003  format("Ip(MA), Bt(T), R(m), a(m), k, d ",6(2x,f6.2))
      call setold(-(temax-temin)*.15_R8,(denmax-denmin)*1.18_R8,1,0,1,0)    
      write(s100,1004) times,alpt,alpn,alpj,ctau
      call gtext(s100,80,0)
1004  format("t(sec), alpt,  alpn , alpj, ctau      ",5(2x,f6.2))
      call setold(-(temax-temin)*.15_R8,(denmax-denmin)*1.13_R8,1,0,1,0)    
      write(s100,1005) ncontr,conval(1),conval(ncontr)
      call gtext(s100,80,0)                                                
1005  format(i3, "paux contour values equally spaced between ",f6.2,     &
     &"  and  "f6.2)
      call setold(-(temax-temin)*0.15_R8,(denmax-denmin)*1.08_R8,1,0,1,  &  
     & 0)
      write(s100,1006) coh,cal,cbr,csyn,cwtot
      call gtext(s100,80,0)
1006  format("coh,  cal,  cbr,  csyn,  cwtot",5(2x,f6.2))
      call rcontr(ncontr,conval,7,cpaux,ntemp,tkvary,1,ntemp,1,          &  
     &   denary,1,ndens,1)
!..plot beta limit
      call setpch(0,3,-1,-1)
      call tracec("beta",tkvary,betalm,ntemp,-1,-1,0._R8,0._R8)
      call tracec("murakami",tkvary,denlim,ntemp,-1,-1,0._R8,0._R8)
!..plot actual position
      call pointc("x",avtemp,avdens/1.E20_R8,1,                          &  
     &           -1,-1,0.0_R8,0.0_R8)
      write(88,2013) kcycle
2013  format("popcon plot, cycle =",i6)
      call frscj(6)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
