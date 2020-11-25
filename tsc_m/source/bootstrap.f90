      subroutine bootstrap(ajavbs, iboot, ierr, nout,                    &  
     &  zeff,rm2avx,epsx,rmajx,qx,ftrapx,                                &  
     &  pvale,ppvale,tempxe,tempxpe,                                     &  
     &  pvali,ppvali,tempxi,tempxpi)
!                                             6/27/94 scj
!
!     DESCRIPTION OF RETURNED PARAMETERS:
!     ajavbs   surface averaged parallel current due to bootstrap
!              divided by surface averaged B dot grad phi
!              (divide by mu0=4*pi*1.e-7 to get amperes/meter)
!     ierr     =0 for normal calculation
!              =1 for error exit
!
!     DESCRIPTION OF INPUT PARAMETERS:
!     iboot    =1 for collisionless bootstrap calculation
!              =2 to include collisonal corrections
!              =3 Sauter model
!     nout     logical unit for error output
!
!     zeff     effective charge
!     rm2avx   surface average of 1/R**2
!     epsx     inverse aspect ratio of this flux surface
!     rmajx    major radius of this flux surface (m)
!     qx       safety factor
!     ftrapx   trapped particle fraction
!     pvale    electron pressure times mu0
!     ppvale   derivative of electron pressure wrt psi  (PF/radian)
!     pvali    ion pressure times mu0
!     ppvali   derivative of ion pressure wrt psi       (PF/radian)
!     tempxe   electron temperature (kev)
!     tempxpe  derivative of electron pressure wrt psi  (PF/radian)
!     tempxi   ion temperature      (kev)
!     tempxpi  derivative of ion temperature wrt psi    (PF/radian)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iboot,ierr,nout,igoof
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zeff,rm2avx,epsx,rmajx,qx,ftrapx,pvale,ppvale,tempxe
      REAL*8 tempxpe,pvali,ppvali,tempxi,tempxpi,ajavbs,amu0,pi
      REAL*8 echarge,emass,pmass,xeps0,zave,xf,alfi,denom,anum1
      REAL*8 anum2,pepope,pipopi,tepote,tipoti,a13,b13,c13,a23,b23
      REAL*8 c23,xmi,xclog,xdense,xdensi,xvthe,xvthi,xtaue,xtaui
      REAL*8 xnuse,xnusi,xl31,xl32,col1,col2,alfhh,ft33,ft31,ft32ee
      REAL*8 ft32ei,ft34,alp0,alp,fl32ee,fl32ei,xl34,ratiop
!============
      data amu0/1.256637061E-6_R8/
      data pi/3.1415926535_R8/
      data echarge/1.6022E-19_R8/
      data emass/9.1095E-31_R8/
      data pmass/1.6726E-27_R8/
      data xeps0/8.854E-12_R8/
      data igoof/0/
!
      ierr = 0
      if(tempxe.le.0 .or. tempxi.le.0) go to 101
      if(ftrapx.le.0 .or. ftrapx.ge.1) go to 102
      if(pvale .le.0 .or. pvali .le.0) go to 103
      if(rm2avx .le. 0) go to 104
      if(zeff .lt. 1.0_R8) go to 105
!     if(iabs(iboot).ne.1 .or. iabs(iboot).ne.2
!    + .or. iabs(iboot) .ne. 4) go to 106
      if(iboot .lt. 1 .or. iboot .gt. 3) go to 106
      zave = (pvale/tempxe)/(pvali/tempxi)
      xf=ftrapx/(1.0_R8-ftrapx)
      alfi = 1.172_R8/(1.0_R8+0.462_R8*xf)
      denom = 1.414_R8*zeff + zeff**2 + xf*(0.754_R8+2.657_R8*zeff+      &  
     & 2._R8*zeff**2)                                                    &  
     &                          + xf**2*(0.348_R8+1.243_R8*zeff+         &  
     & zeff**2)
      anum1 = xf*(0.754_R8+ 2.21_R8*zeff + zeff**2                       &  
     &      + xf*(0.348_R8+ 1.24_R8*zeff + zeff**2) )
      anum2 = xf*(0.884_R8+ 2.07_R8*zeff)
      pepope=ppvale/pvale
      pipopi=ppvali/pvali
      tepote=tempxpe/tempxe
      tipoti=tempxpi/tempxi
      if(iboot .eq. 1) then
      ajavbs=-(pvale/rm2avx)*(anum1*(pepope+(pvali/pvale)                &  
     &          *(pipopi - alfi*tipoti) ) - anum2*tepote) / denom
      return
      endif
!
      if(zeff .le. 2.0_R8) then
        a13=1.02_R8-0.13_R8*(zeff-1.0_R8)
        b13=1.07_R8-0.45_R8*(zeff-1.0_R8)
        c13=1.07_R8-0.38_R8*(zeff-1.0_R8)
        a23=0.57_R8-0.05_R8*(zeff-1.0_R8)
        b23=0.61_R8-0.27_R8*(zeff-1.0_R8)
        c23=0.61_R8-0.23_R8*(zeff-1.0_R8)
      else
        a13=0.89_R8-0.05_R8*(zeff-2.0_R8)
        b13=0.62_R8-0.03_R8*(zeff-2.0_R8)
        c13=0.69_R8-0.09_R8*(zeff-2.0_R8)
        a23=0.52_R8-0.02_R8*(zeff-2.0_R8)
        b23=0.34_R8-0.005_R8*(zeff-2.0_R8)
        c23=0.38_R8-0.05_R8*(zeff-2.0_R8)
      endif
      xmi=1.0_R8
      xclog=17.0_R8
      xdense=pvale/(amu0*tempxe*echarge*1.E3_R8)
      xdensi=pvali/(amu0*tempxi*echarge*1.E3_R8)
      xvthe=sqrt(2.0_R8*echarge*1.E3_R8*tempxe/emass)
      xvthi=sqrt(2.0_R8*echarge*1.E3_R8*tempxi/(xmi*pmass))
      xtaue=(12.0_R8*pi**1.5_R8*xeps0**2*(emass)**2*xvthe**3)/           &  
     &(4.0_R8*xdense*zeff*(echarge)**4*xclog)
      xtaui=1.414_R8*(12.0_R8*pi**1.5_R8*xeps0**2*(xmi*pmass)**2*        &  
     &xvthi**3)/(4.0_R8*xdensi*(zave*echarge)**4*xclog)
      xnuse=1.414_R8*(rmajx*qx)/(xtaue*xvthe*epsx**1.5_R8)
      xnusi=1.414_R8*(rmajx*qx)/(xtaui*xvthi*epsx**1.5_R8)
      xl31=anum1/denom
      xl32=anum2/denom
      col1=(1.0_R8/(1.0_R8+a13*sqrt(xnuse)+b13*xnuse))*                  &  
     &(1.0_R8/(1.0_R8+c13*xnuse*epsx**1.5_R8))
      col2=(1.0_R8/(1.0_R8+a23*sqrt(xnuse)+b23*xnuse))*                  &  
     &(1.0_R8/(1.0_R8+c23*xnuse*epsx**1.5_R8))
      alfhh=((alfi-0.35_R8*sqrt(xnusi))/(1.0_R8+0.7_R8*sqrt(xnusi))      &  
     &-2.1_R8*xnusi**2*epsx**3)*(1.0_R8/(1.0_R8+xnusi**2*epsx**3))*      &  
     &(1.0_R8/(1.0_R8+xnuse**2*epsx**3))
      if(iboot .eq. 2) then
        ajavbs=-(pvale/rm2avx)*(xl31*col1*(pepope+(pvali/pvale)            &  
     &  *(pipopi-alfhh*tipoti))-(2.5_R8*xl31*col1-(2.5_R8*xl31-xl32)*col2)  &  
!    &                                                                   &  
     &  *tepote)
        return
      endif
!
!.....IBOOT = 3
!
      ft33=ftrapx/(1._R8+(0.55_R8-0.1_R8*ftrapx)*sqrt(xnuse)+0.45_R8*    &  
     & (((1._R8-ftrapx)*xnuse)/zeff**1.5_R8))
!
      ft31=ftrapx/(1._R8+(1._R8-0.1_R8*ftrapx)*sqrt(xnuse)+0.5_R8*       &  
     & (((1._R8-ftrapx)*xnuse)/zeff))
!
      ft32ee=ftrapx/(1._R8+0.26_R8*(1._R8-ftrapx)*sqrt(xnuse)+0.18_R8*   &  
     & (((1._R8-0.37_R8*ftrapx)*xnuse)/sqrt(zeff)))
!
      ft32ei=ftrapx/(1._R8+(1._R8+0.6_R8*ftrapx)*sqrt(xnuse)+0.85_R8*    &  
     & (((1._R8-0.37_R8*ftrapx)*xnuse)*(1._R8+zeff)))
!
      ft34=ftrapx/(1._R8+(1._R8-0.1_R8*ftrapx)*sqrt(xnuse)+0.51_R8*      &  
     & (((1._R8-0.5_R8*ftrapx)*xnuse)/zeff))
!
      alp0=-(1.17_R8*(1._R8-ftrapx))/(1._R8-0.22_R8*ftrapx-0.19_R8*      &  
     & ftrapx**2)
!
      alp=(((alp0+0.25_R8*(1._R8-ftrapx**2)*sqrt(xnusi))/                &  
     & (1._R8+0.5_R8*sqrt(xnusi)))                                       &  
!    & - 0.315_R8*xnusi**2*ftrapx**6)*(1._R8/(1._R8+0.15_R8*xnusi**2*    &  
!      See Erratum by Sauter, et al, Phys Plasmas Vol 9 p 5140 (2002)
     & + 0.315_R8*xnusi**2*ftrapx**6)*(1._R8/(1._R8+0.15_R8*xnusi**2*    &  
     & ftrapx**6))
!
      xl31=(1._R8+(1.4_R8/(1._R8+zeff)))*ft31-(1.9_R8/(1._R8+zeff))*     &  
     & ft31**2                                                           &  
     & +(0.3_R8/(1._R8+zeff))*ft31**3+(0.2_R8/(1._R8+zeff))*ft31**4
!
      fl32ee=(((0.05_R8+0.62_R8*zeff)/(zeff*(1._R8+0.44_R8*zeff)))*      &  
     & (ft32ee-ft32ee**4)                                                &  
     & +(1._R8/(1._R8+0.22_R8*zeff))*(ft32ee**2-ft32ee**4-1.2_R8*        &  
     & (ft32ee**3                                                        &  
     & -ft32ee**4))                                                      &  
     & +(1.2_R8/(1._R8+0.5_R8*zeff))*ft32ee**4)
!
      fl32ei=((-(0.56_R8+1.93_R8*zeff)/(zeff*(1._R8+0.44_R8*zeff)))*     &  
     & (ft32ei                                                           &  
     & -ft32ei**4)                                                       &  
     & +(4.95_R8/(1._R8+2.48_R8*zeff))*(ft32ei**2-ft32ei**4-0.55_R8*     &  
     & (ft32ei**3                                                        &  
     & -ft32ei**4))                                                      &  
     & -(1.2_R8/(1._R8+0.5_R8*zeff))*ft32ei**4)
!
      xl32=fl32ee+fl32ei
!
      xl34=(1._R8+(1.4_R8/(1._R8+zeff)))*ft34-(1.9_R8/(1._R8+zeff))*     &  
     & ft34**2                                                           &  
     & +(0.3_R8/(1._R8+zeff))*ft34**3+(0.2_R8/(1._R8+zeff))*ft34**4
!
      ratiop = (1.-(pvale/(pvale+pvali)))/(pvale/(pvale+pvali))

      ajavbs=-(pvale/rm2avx)*(xl31*(pepope+(pvali/pvale)*pipopi)+        &  
!     & xl32*tepote+xl34*alp*tipoti)
!     changed 02/04/2010 SCJ
     & xl32*tepote+xl34*alp*ratiop*tipoti)

      return
!
!.....error exit
  101 continue
      write(nout,1101) tempxe,tempxi
 1101 format(" Error in bootstrap, tempxe,tempxi=",1p2e12.4)
      go to 200
  102 continue
      write(nout,1102) ftrapx
 1102 format(" Error in bootstrap, ftrapx=",1pe12.4)
      go to 200
  103 continue
      write(nout,1103) pvale,pvali
 1103 format(" Error in bootstrap, pvale,pvali =",1p2e12.4)
      go to 200
  104 continue
      write(nout,1104) rm2avx
 1104 format(" Error in bootstrap, rm2avx=",1pe12.4)
      go to 200
  105 continue
      write(nout,1105) zeff
 1105 format(" Error in bootstrap, zeff=",1pe12.4)
      go to 200
  106 continue
      write(nout,1106) iboot
 1106 format(" Error in bootstrap, iboot=",i5)
!
  200 igoof=igoof+1
      if(igoof.gt.10) ierr=1
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
