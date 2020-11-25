      subroutine tieval(psv,itype,tival,ig,jg,pval,eval,rval)
!
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.....for itype=1  evaluate from analytic function
!.......  itype=2 evaluate from mu profile
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,ig,jg,lval,lv,lm
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tival,pval,eval,rval,psv,psve
!============
      if(iexvc(ig,jg).gt.0 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
      if(itype.eq.2) go to 10
      if(itype.le.0) go to 100
      if(itype.gt.2) go to 100
!
      tival = (pval-eval)*udsh /                                         &  
     &      (rval*udsd/zgas)
!
      return
   10 continue
!
!.....interpolate for surface averaged transport
      psve = max(psv,xsv(2))
      lval = 2
      do 20 lv=3,npsit
      lval = lv
      if(xsv(lval).gt.psve) go to 30
   20 continue
      go to 40
   30 continue
      lm = lval-1
      tival = ((xsv(lv)-psve)*ti(lm) + (psve-xsv(lm))*ti(lv))            &  
     &      / (xsv(lv)-xsv(lm))
      if(tival.lt. tevv) tival = tevv
      if(psve.lt.xsv(2)) tival = ti(2)
      return
!
   40 continue
!
!.....vacuum region
      tival = tevv
      return
  100 continue
!.....error exit
      ineg=51
      tival = 0._R8
      write(nout,1000) psv,psimin,psilim
 1000 format(" error in tieval,psv,psimin,psilim=",1p3e12.4)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
