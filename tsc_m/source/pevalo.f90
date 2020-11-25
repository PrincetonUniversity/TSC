      subroutine pevalo(psv,itype,pval,ppval,ig,jg)
!.....number 8.50
!
      USE CLINAM
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.....for itype=1  evaluate from analytic function
!.......  itype=2 evaluate from mu profile
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,ig,jg,n1,n2,lval,lv,lm
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pval,ppval,psv,phiarg,phix,ea,eax,psve,pvalmin,psl
!============
      if(iexvc(ig,jg).gt.0 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
      if(itype.le.0) go to 100
      if(itype.gt.2) go to 100
      if(itype.eq.2) go to 10
      if(ifunc.eq.6) go to 110
!
      phiarg = (psilim-psv)*delpsi
!
      if(ifunc.eq.2) go to 5
      phiarg = abs(1._R8- ((psv-psimin)*delpsi))
!
      n1 = 2
      n2 = 1
      pval = p0*phiarg**alphap*(1._R8+acoef(110)*phiarg*alphap           &  
     &                            /(alphap+1))
      ppval =-p0*alphap*phiarg**(alphap-1._R8)                           &  
     &         *(1._R8+acoef(110)*phiarg)*delpsi
!
      return
!
!.....ifunc=2  tokamak profiles (ORNL)
    5 continue
      phix = 1._R8-phiarg
      ea = exp(-alphap)
      eax = exp(-alphap*phix)
      if(alphap.eq.0) go to 41
      pval = p0*(ea*(1._R8+(1._R8-phix)*alphap)-eax)/(ea-1)/alphap       &  
     &     + acoef(881)*r0*smallt*acoef(882)                             &  
     &      + acoef(892)*p0*phiarg**2**(1._R8-phiarg)
      ppval = p0*(eax-ea)/(ea-1._R8)*delpsi                              &  
     &     - acoef(892)*p0*(2._R8*phiarg - 3._R8*phiarg**2)*delpsi
      go to 50
   41 pval = .5_R8*p0*(1._R8-phix)**2                                    &  
     &     + acoef(881)*r0*smallt*acoef(882)                             &  
     &      + acoef(892)*p0*phiarg**2*(1._R8-phiarg)
      ppval = -p0*delpsi*(1._R8-phix)                                    &  
     &     - acoef(892)*p0*(2._R8*phiarg - 3._R8*phiarg**2)*delpsi
   50 continue
      return
   10 continue
!
!.....get p from interpolating polynomial for surface averaged transport
      psve = max(psv,xsv(2))
      lval = 2
      do 20 lv=3,npsit
      lval = lv
      if(xsv(lval).gt.psve) go to 30
   20 continue
      go to 40
   30 continue
      lm = lval-1
      pval=asv(lm,2)+psve*(bsv(lm,2)+psve*(csv(lm,2)+psve*dsv(lm,2)))
      pvalmin = acoef(881)*acoef(882)*r0*smallt
      if(pval .lt. pvalmin) pval = pvalmin
      ppval = bsv(lm,2)+psve*(2._R8*csv(lm,2)+3._R8*dsv(lm,2)*psve)
      return
  110 continue
      psl =  (psv-psimin)*delpsi
!
!
!.....get p from interpolating polynomial for ifunc=6
      lval = 2
      do 120 lv=3,npsit6+1
      lval = lv
      if(xs6(lval).gt.psl) go to 130
  120 continue
      go to 40
  130 continue
      lm = lval-1
      pval=as6(lm,2)+psl*(bs6(lm,2)+psl*(cs6(lm,2)+psl*ds6(lm,2)))
      pvalmin = acoef(881)*acoef(882)*r0*smallt
      if(pval .lt. pvalmin) pval = pvalmin
      ppval = bs6(lm,2)+psl*(2._R8*cs6(lm,2)+3._R8*ds6(lm,2)*psl)
      return
!
   40 continue
!
!.....vacuum region
      pval = acoef(881)*r0*smallt*acoef(882)
      ppval = 0._R8
      return
  100 continue
!.....error exit
      ineg=53
      pval = 0._R8
      ppval = 0._R8
      write(nout,1000) psv,psimin,psilim
 1000 format(" error in peval,psv,psimin,psilim=",1p3e12.4)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
