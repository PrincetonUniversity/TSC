      subroutine control(AIP,TT,TVDE,ZPP,RPP,WVSPAIP,ZXP,SHAPEP,GAPINP,  &  
     &                  ZP00,zp,DFZP0,dfzp,TAY)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE DIAGCOM
      USE LOOP1
      USE PROBE1
      USE ARRAYS

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ktime,kpr1
      INTEGER int
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tt,tvde,zpp,rpp,wvspaip,zxp,shapep,gapinp,zp00,zp
      REAL*8 dfzp0,dfzp,tay,aip
      REAL*8 psf5atopsf8a,psf4atopsf8a,psf2atopsf8a,psf9btopsf8b
      REAL*8 psf5btopsf8b,psf4btopsf8b,psf3btopsf8b,dfz1b,dfz1,zcur
      REAL*8 f9a,zppv,vzp,dfz2,dfz,vertbo,dfzo,delta,vertao,dfr1
      REAL*8 dfr2,dfr,zb66,psf6ab,frat,prebip,f6amod,f6bmod,bias7p
      REAL*8 f7amod,f7bmod,sdiv,si,sz,zx1a,zx1,zx,dfzx,zb,dfgap1
      REAL*8 gapinpo,dfgap,btop,dftop1,gaptop,prtopa,prtop
!============
!     common/diagcom/verta,vertb,dfzpd
!     common / arrays / PF(18),PFE(18),FCOM(18),vchopper(18),FCOMSV(18)
!       COMMON /loop1/
!    *PSF2A,PSF4A,PSF5A,PSF6NA,
!    *PSF8A,PSF1B,PSF2B,PSF3B,PSF4B,PSF5B,PSF6NB,
!    *PSF8B,PSF9B,PSI7A,PSF6FA,PSI7B,PSF6FB
!     COMMON /probe1/
!    *AMPI9A ,AMPI67A ,AMPI6FA ,AMPI66M ,AMPI6FB ,
!    *AMPI67B ,AMPI7FB ,AMPI9B ,AMPI89B ,
!    *AMPI8B ,AMPI4B ,AMPI3B ,AMPI2B ,AMPI1B ,
!    *AMPI11M ,AMPI1A ,AMPI2A ,AMPI5A
!
      data ktime/0/
 
!.........
! Values for flux ratios for shot 83964 @ 1800 msec
!
      psf5atopsf8a = 1.001_R8
      psf4atopsf8a = 0.79956_R8
      psf2atopsf8a = 0.02136_R8
      psf9btopsf8b = 0.64697_R8
      psf5btopsf8b = 0.69885_R8
      psf4btopsf8b = 0.25940_R8
      psf3btopsf8b = 0.04883_R8
 
!
!
!  dfz CONTROL
      kpr1=0
!_____________________________________________________________
      dfz1b=0.532_R8*(13.69_R8*PSI7A)-0.537_R8*(13.5641_R8*PSI7B)        &  
     &  +1.506_R8*(13.5928_R8*PSF2A)-1.467_R8*(13.9567_R8*PSF2B)
!
      dfz1=2.446_R8*(5.91_R8*AMPI67A )-2.469_R8*(5.91_R8*AMPI67B )+      &  
     &  1.181_R8*(5.85_R8*AMPI2A )-1.178_R8*(5.87_R8*AMPI2B )+dfz1b
!
      zcur=2._R8*dfz1/(AIP/5.E5_R8)
      zp=10._R8*zcur
!
!
      f9a=0.300_R8
      zppv=ZPP-5._R8*PF(9)*4.3723E-1_R8
!---------------------
!
      vzp=(zp-ZP00)/TAY
      dfz2=0.2_R8*(2.5_R8*zppv)*(AIP/5.E5_R8)
      dfz=dfz1-dfz2
      if(ktime.ne.0) vertb = 0.5_R8*(vertb+vertbo)
      vertbo = vertb
      if(ktime.ne.0) dfz = 0.5_R8*(dfz+dfzo)
      dfzo = dfz
      delta=20._R8*dfz/(AIP/5.E5_R8)
!
!----------------------------------------------------------
      dfzp=20._R8*dfz
      dfzpd=(dfzp-DFZP0)/(TAY*1.E-3_R8)
!-------
      IF(dfzp.GT.14._R8)THEN
      dfzpd=0._R8
      dfzp=14._R8
      END IF
!-------
!
      verta=1._R8/20._R8*( dfzp+1.E-2_R8*dfzpd)
      if(ktime.ne.0) verta = 0.5_R8*(verta+vertao)
      vertao = verta
!
      vertb=1._R8/20._R8*( 4.E-3_R8*dfzpd )
      if(ktime.ne.0) vertb = 0.5_R8*(vertb+vertbo)
      vertbo = vertb
!----------------------------------------
      if(TT.ge.TVDE)then
      verta=0._R8
      vertb=0._R8
      end if
!-----------
!-------------------------------------------------
!  dfr CONTROL
!___________________________________________________________________
!
!
      dfr1=1._R8*(4.93_R8*AMPI11M )+0.842_R8*(5.85_R8*AMPI2A )           &  
     &  +0.820_R8*(6._R8*AMPI1B )
!
      dfr2=0.342_R8*(5.72_R8*AMPI6FA )+0.333_R8*(5.88_R8*AMPI6FB )       &  
     &  +1.063_R8*(5.94_R8*AMPI66M )
!
      dfr=2._R8*( -dfr2+0.1_R8*(1.034_R8*dfr1)*(2.5_R8*RPP) )
      delta=(19.343_R8/dfr1)*(13.5955_R8*PSF6NA+0.997_R8*13.5986_R8*     &  
     & PSF6NB                                                            &  
     &  -13.6_R8*PSF1B-dfr)
!__________________________________________________________
!
! f6amod, f6bmod, f7amod, f7bmod
!_____________________________________________________
!
      zb66=(2._R8/10._R8)*(2._R8*5.94_R8*AMPI66M )*(2.5_R8*zppv)
!
      psf6ab=13.5955_R8*PSF6NA+13.5986_R8*PSF6NB
!
!  *** WVSPAIP IS CONTROL PARAMETER
      frat=(2._R8/10._R8)*WVSPAIP*psf6ab
      prebip=0._R8
!----------------
!
      f6amod=2._R8*(13.5955_R8*PSF6NA-0.5_R8*13.6_R8*PSF1B               &  
     &  +vertb-0.15_R8*zb66-prebip)
!
      f6bmod=2._R8*(13.5986_R8*PSF6NB-0.5_R8*13.6_R8*PSF1B-vertb+        &  
     & 0.15_R8*zb66)
!
      bias7p=0._R8
!------------
      f7amod=1.5_R8*(13.69_R8*PSI7A)-2.25_R8*frat+verta-bias7p
      f7bmod=1.5_R8*(13.5641_R8*PSI7B)-2.25_R8*frat-bias7p
!....note: f7bmod zeroed out by M-Matrix....Jan12 1997
      f7bmod = 0._R8
!
!_______________________________________________________________
!
!   DFZXB TOM X-POINT HEIGHT
!________________________________________________________________
!  PROGRAMMED VERTICAL POSITION OF LOWER X-POINT
!
      sdiv=-6.279_R8*(5.84_R8*AMPI89B )+1.330_R8*(5.9_R8*AMPI8B )        &  
     &  +1.801_R8*(5.82_R8*AMPI9B )-                                     &  
     &  0.405_R8*(5.85_R8*AMPI4B )+(0.738_R8)*(AIP/5.E5_R8)
!
      si=1.241_R8*(AIP/5.E5_R8)+0.5_R8*sdiv
!
      sz=(0.0185_R8*zcur)*(sdiv)
!
      zx1a=0.380_R8*(AIP/5.E5_R8)+sz+1.786_R8*sdiv
!
      zx1=zx1a/si
      zx=100._R8*(0.5_R8*zx1-1.746_R8)
      int=int+1
!
!
      dfzx=2._R8*(-0.811_R8*(AIP/5.E5_R8)+0.2_R8*sz-0.08_R8*ZXP*si)
      delta=125._R8*dfzx/si
!
!
!   INSIDE GAP CONTROL
!________________________________________________________
      zb=zx1*(5.84_R8*AMPI89B )
!
!
      dfgap1=0.331_R8*(13.6411_R8*PSF8B)+0.333_R8*(13.5792_R8*PSF9B)+    &  
     & 0.604_R8*zb
!  INNER GAP PROGRAM
!  *** gapinpo CM =10*GAPINP V -10.9 CM-->
      gapinpo=10._R8*GAPINP-10.9_R8
      dfgap=dfgap1-0.172_R8*5.81_R8*AMPI1A *GAPINP
      delta=10._R8*dfgap/(0.172_R8*5.81_R8*AMPI1A )
!
!
!      TOP GAP
!____________________________________________________________
      btop=5.88_R8*AMPI9A +0.473_R8*(5.87_R8*AMPI5A )
      dftop1=0.107_R8*(13.6411_R8*PSF8B)+0.108_R8*(13.5792_R8*PSF9B)-    &  
     & 1.86_R8*zb
      gaptop=(dftop1-13.5817_R8*PSF8A)/(0.323_R8*btop)-18.326_R8*zcur+   &  
     & 34.836_R8
      if(abs(gaptop).gt.10000._R8)then
      end if
!
      prtopa=SHAPEP-1.833_R8*zcur+3.484_R8
      prtop=5._R8*(0.2_R8*dftop1+0.1_R8*(6.456_R8*btop)*prtopa)
      delta=(13.5817_R8*PSF8A-prtop)*(-3.098_R8)/btop
!   commands.........
!*********************************************************
!*********************************************************
      FCOM(1)=-200._R8*dfgap
      FCOM(2)=100*(vertb+13.5928_R8*PSF2A-psf2atopsf8a*13.5817_R8*PSF8A)    
      FCOM(3)=0._R8
      FCOM(4)=100*(psf4atopsf8a*13.5817_R8*PSF8A-13.8136_R8*PSF4A)
      FCOM(5)=100*(psf5atopsf8a*13.5817_R8*PSF8A-13.8235_R8*PSF5A)
      FCOM(6)=200*(f6amod-dfr)
      FCOM(7)=200*(f7amod)
      FCOM(8)=200*(prtop-13.5817_R8*PSF8A)
      FCOM(9)=0._R8
      FCOM(10)=100*(13.6_R8*PSF1B)
      FCOM(11)=100*(-vertb+13.9567_R8*PSF2B)
      FCOM(12)=100*(12.9165_R8*PSF3B-psf3btopsf8b*13.6411_R8*PSF8B)
      FCOM(13)=100*(psf4btopsf8b*13.6411_R8*PSF8B-13.8734_R8*PSF4B)
      FCOM(14)=100*(psf5btopsf8b*13.6411_R8*PSF8B-13.7418_R8*PSF5B)
      FCOM(15)=200*(f6bmod-dfr)
      FCOM(16)=200*(f7bmod-verta)
      FCOM(17)=100*(13.5792_R8*PSF9B-psf9btopsf8b*13.6411_R8*PSF8B)
      FCOM(18)=200._R8*(verta-dfzx)
!
!
!
!
!.....special printout
!     if (mod(ktime,10) .eq. 0)
!    1write(96,6666) ktime,fcom(6),fcom(7),fcom(15),fcom(16),
!    1     f6amod, f7amod,f6bmod,f7bmod,dfr,verta,vertb,
!    2     frat,zb66,dfzp ,dfzpd,
!    3     dfz1,dfz2,dfz1b,zppv,aip,
!    1     WVSPAIP,psf6ab,PSF6NA,PSF6NB
!6666 format("  special printout....cycle #", i5,/
!    1" fcom(6),(7),(15),(16) =       ",1p4e12.4,/,
!    2" f6amod,f7amod,f6bmod,f7bmod = ",1p4e12.4,/,
!    3" dfr,verta,vertb frat =        ",1p4e12.4,/,
!    4" zb66,dfzp,dfzpd,dfz1 =        ",1p4e12.4,/,
!    5" dfz2,dfz1b,zppv,aip  =        ",1p4e12.4,/,
!    6" WVSPAIP,psf6ab,PSF6NA,PSF6NB =",1p4e12.4)
!     ktime = ktime + 1
!     if(ktime.gt.100) stop
!
      RETURN
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
