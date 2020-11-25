      subroutine forms(formx,formne,formte,formj)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE CB_ALL

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 formne,formte,formj,formx
      REAL*8 dxform,rda
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
      dimension formx(*),formne(*),formte(*),                            &  
     &   formj(*)
!============      
!..set up arrays defining profile forms for r/a, ne/ne0, te/te0, and j/j0
      do 1 i=1,nform
      dxform=1._R8/(nform-1)
      rda=(i-1)*dxform+1.E-8_R8
      formx(i)=rda
      formne(i)=(1._R8+1.E-6_R8-rda**2)**alpn
      formte(i)=(1._R8+1.E-6_R8-rda**2)**alpt
      formj(i)=(1._R8+1.E-6_R8-rda**2)**alpj
1     continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
