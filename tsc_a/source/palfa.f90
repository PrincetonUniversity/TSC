      function palfa(pscale,formte,formx,formne)
!
      USE CB_ALL
      IMPLICIT NONE
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,nforms
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pscale,formte,formx,formne,palfa
      REAL*8 sgvn2,tval,sigdt,dx,palhv
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
      dimension formte(*),formx(*),formne(*)
      dimension sgvn2(1000)
!============      
      do 1 i=1,nform
      tval=tkv0*formte(i)
!
      if(abs(tval).lt.1.E-10_R8) then
      sgvn2(i)=0._R8
      goto 1
      end if
!
      sigdt=exp(-21.377692_R8/tval**.2935_R8-25.204054_R8                &  
     &   -7.1013427e-2*tval+1.9375451e-4*tval**2                         &  
     &   +4.9246592e-6*tval**3-3.9836572e-8*tval**4)*1.e-6
      sgvn2(i)=2._R8*formx(i)*formne(i)**2*sigdt
1     continue
!
!..perform radial integral
!
      dx=1._R8/(nform-1)
      nforms=nform
      call simson(palhv,sgvn2,dx,nforms)
      palfa=pscale*palhv
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
