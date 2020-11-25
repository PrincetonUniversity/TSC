      subroutine peval(psv,itype,pval,ppval,ig,jg)
!
      USE CLINAM
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.....for itype=1  evaluate from analytic function
!.......  itype=2 evaluate from mu profile
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,ig,jg,n1,n2,lval,lv,lm
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pval,ppval,psv,phiarg,phix,ea,eax,pconst,xparam,psve
      REAL*8 pvalmin,psl
!============
      if(itype.lt.1 .or. itype.gt.2) go to 100
!
      if(itype.eq.2) go to 10
!
!.....check if reading data from MDSPLUS Archive using trxpl
      if(acoef(4993).gt.0) go to 10
!
      if(ifunc.eq.6) go to 110
!
      phiarg = (psilim-psv)*delpsi
!
      if(ifunc.eq.2) go to 5
      if(iexvc(ig,jg).gt.0 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
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
!
!.....added 12/30/99 for FRC calculations
      if(acoef(894).gt.0) go to 51
      if(iexvc(ig,jg).gt.0 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
      phix = 1._R8-phiarg
      ea = exp(-alphap)
      eax = exp(-alphap*phix)
      if(alphap.eq.0) go to 41
      pval = p0*(ea*(1._R8+(1._R8-phix)*alphap)-eax)/(ea-1)/alphap       &  
     &     + fracn0*r0*smallt*acoef(882)                             &
     &      + acoef(892)*p0*phiarg**2**(1._R8-phiarg)
      ppval = p0*(eax-ea)/(ea-1._R8)*delpsi                              &  
     &     - acoef(892)*p0*(2._R8*phiarg - 3._R8*phiarg**2)*delpsi
      go to 50
   41 pval = .5_R8*p0*(1._R8-phix)**2                                    &  
     &     + fracn0*r0*smallt*acoef(882)                             &
     &      + acoef(892)*p0*phiarg**2*(1._R8-phiarg)
      ppval = -p0*delpsi*(1._R8-phix)                                    &  
     &     - acoef(892)*p0*(2._R8*phiarg - 3._R8*phiarg**2)*delpsi
   50 continue
      return
   51 continue
!
!.....special for FRC calculations
      pconst = p0*(acoef(894)*(tanh(acoef(891)+acoef(892))+1._R8)        &  
     &       - (tanh(acoef(892))+1._R8))/(2._R8*(1._R8-acoef(894)))
      xparam = (psv-psilim)/(psimin-psilim)
      pval = .5_R8*p0*(tanh(acoef(891)*xparam + acoef(892)) + 1._R8) +   &  
     & pconst
      ppval =.5_R8*p0*(   (1._R8/cosh(acoef(891)*xparam + acoef(892) ) )  &  
     & **2                                                               &  
     &      *acoef(891)/(psimin-psilim))
      return
   10 continue
      if(iexvc(ig,jg).gt.0 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
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
      pvalmin = fracn0*acoef(882)*r0*smallt
      if(pval .lt. pvalmin) pval = pvalmin
      ppval = bsv(lm,2)+psve*(2._R8*csv(lm,2)+3._R8*dsv(lm,2)*psve)
      return
  110 continue
      if(iexvc(ig,jg).gt.0 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
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
      pvalmin = fracn0*acoef(882)*r0*smallt
      if(pval .lt. pvalmin) pval = pvalmin
      ppval = bs6(lm,2)+psl*(2._R8*cs6(lm,2)+3._R8*ds6(lm,2)*psl)
      return
!
   40 continue
!
!.....vacuum region
      pval = fracn0*r0*smallt*acoef(882)
      ppval = 0._R8
      return
  100 continue
!.....error exit
      ineg=52
      pval = 0._R8
      ppval = 0._R8
      write(nout,1000) psv,psimin,psilim
 1000 format(" error in peval,psv,psimin,psilim=",1p3e12.4)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
