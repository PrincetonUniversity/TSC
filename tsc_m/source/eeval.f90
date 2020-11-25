      subroutine eeval(psv,itype,eval,epval,ig,jg)
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
      INTEGER itype,ig,jg,n3,n1,n2,lval,lv,lm
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 eval,epval,psv,phiarg,phix,ea,eax,psve,psl,evalmin
!============
      if(iexvc(ig,jg).gt.0 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
!
!.....check if reading data from MDSPLUS Archive using trxpl
      if(acoef(4992).gt.0) go to 10
      if(itype.le.0) go to 100
      if(itype.gt.2) go to 100
      if(itype.eq.2) go to 10
      if(ifunc.eq.6) go to 110
!
      phiarg = (psilim-psv)*delpsi
!
      if(ifunc.eq.2) go to 5
      n3 = 1._R8
      phiarg = (1._R8- ((psv-psimin)*delpsi)**n3)
      n1 = 2
      n2 = 1
      eval = e0*phiarg**alphae*(1._R8+acoef(110)*phiarg*alphae           &  
     &                            /(alphae+1))                           &  
     &     + fracn0*r0*smallt                                        &
     &     + acoef(892)*e0*phiarg**n1*(1._R8-phiarg)**n2
      epval =(-e0*alphae*phiarg**(alphae-1._R8)                          &  
     &         *(1._R8+acoef(110)*phiarg)     - acoef(892)*e0            &  
     &  *(n1-(n1+n2)*phiarg)*(1._R8-phiarg)**(n2-1)*phiarg**(n1-1))      &  
     &   *n3*((psv-psimin)*delpsi)**(n3-1)*delpsi
!
      return
!
!.....ifunc=2  tokamak profiles (ORNL)
    5 continue
      phix = 1._R8-phiarg
      ea = exp(-alphae)
      eax = exp(-alphae*phix)
      if(alphae.eq.0) go to 41
      eval = e0*(ea*(1._R8+(1._R8-phix)*alphae)-eax)/(ea-1)/alphae       &  
     &      + fracn0*r0*smallt                                       &
     &      + acoef(892)*e0*phiarg*(1._R8-phiarg)
      epval = e0*(eax-ea)/(ea-1._R8)*delpsi                              &  
     &     - acoef(892)*e0*(2._R8- 3._R8*phiarg)*phiarg**delpsi
      go to 50
   41 eval = .5_R8*e0*(1._R8-phix)**2                                    &  
     &      + fracn0*r0*smallt
      epval = -e0*delpsi*(1._R8-phix)
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
      eval=asv(lm,3)+psve*(bsv(lm,3)+psve*(csv(lm,3)+psve*dsv(lm,3)))
      if(eval .lt. fracn0*r0*smallt) eval = fracn0*r0*smallt
      epval = bsv(lm,3)+psve*(2._R8*csv(lm,3)+3._R8*dsv(lm,3)*psve)
      return
  110 continue
      psl =  (psv-psimin)*delpsi
!
!
!.....get e from interpolating polynomial for ifunc=6
      lval = 2
      do 120 lv=3,npsit6+1
      lval = lv
      if(xs6(lval).gt.psl) go to 130
  120 continue
      go to 40
  130 continue
      lm = lval-1
      eval=as6(lm,3)+psl*(bs6(lm,3)+psl*(cs6(lm,3)+psl*ds6(lm,3)))
      evalmin = fracn0*r0*smallt
      if(eval .lt. evalmin) eval = evalmin
      epval = bs6(lm,3)+psl*(2._R8*cs6(lm,3)+3._R8*ds6(lm,3)*psl)
      return
!
   40 continue
!
!.....vacuum region
      eval = fracn0*r0*smallt
      epval = 0._R8
      return
  100 continue
!.....error exit
      ineg=54
      eval = 0._R8
      epval = 0._R8
      write(nout,1000) psv,psimin,psilim
 1000 format(" error in eeval,psv,psimin,psilim=",1p3e12.4)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
