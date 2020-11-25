      subroutine pauxc(fpaux)
      USE CB_ALL
!
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 fpaux
      REAL*8 tauna,ftau,taukac
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
!     common/cb9/tscoh,tscal,tscbr,tscsync
!     common/cb10/coh,cal,cbr,csyn,cwtot
!..solve power balance for paux
!
!..neoalcator
      tauna=7.E-22_R8*deneln*rminor*rmajor**2*qcyl
!..kaye-all complex
      ftau=ctau*.0521_R8*sqrt(aibar)*kappa**0.25_R8*(curent/1.E6_R8)**   &  
     & .85_R8                                                            &  
     &   *(deneln/1.e19)**0.1*btor**0.3*rminor**0.3                      &  
     &   *rmajor**0.85
      taukac=(ftau*(1.E6_R8/wtot)**.5_R8)**2
      fpaux=wtot/taukac-(fpal+fpoh-fpbr-fpsync)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
