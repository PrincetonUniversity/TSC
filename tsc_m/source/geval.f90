      subroutine  geval (psv, itype, gval, gpval, ggp, ig,jg)
!.....number 8.30
!.......... 8.30 geval
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,ig,jg,icallg,isfa,ip,iz,lval,lv,lm
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 gval,gpval,ggp,psv,fkgrat,fgrat,ggpt,ff,temp,gvalt
      REAL*8 phiarg,eb,ebx,pf,ppf,fsq,denom,facp,facz,vpl,etal
      REAL*8 gxmjal,xmjal,pval,ppval,term1,term2,gpvalt
!============
        data  icallg, fkgrat, fgrat / 0, 0.5_R8, 0.6_R8/
!
!.....for itype=1  evaluate from analytic function
!.......  itype=2 evaluate from q profile
      if(itype.le.0) go to 100
      if(itype.gt.2) go to 100
      if(iexv(ig,jg).eq.1 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
      if(itype.eq.2) go to 10
!
      go to(1,2,3,4,5,1,5),ifunc
!
!.....ifunc=1  Tokamak Profiles (Princeton)
    1 continue
      ggpt = gp1*ff(3,psv,ig,jg) + gp2*ff(4,psv,ig,jg)
      ggp = ggpt
      temp = gzero**2+2._R8*(gp1*ff(1,psv,ig,jg)+gp2*ff(2,psv,ig,jg))
      if(temp.lt.0.0_R8) go to 100
      gvalt = sqrt(temp)
      gval = gvalt
      if(gvalt.eq.0) go to 100
      gpval = ggpt/gvalt
      return
!
!.....ifunc=2  tokamak profiles (ORNL)
    2 continue
      phiarg = (psv-psimin)*delpsi
      eb = exp(-alphag)
      ebx = exp(-alphag*phiarg)
      if(alphag.eq.0) go to 41
      pf = p0*(eb*(1._R8+(1._R8-phiarg)*alphag)-ebx)/(eb-1._R8)/alphag
      ppf = p0*(ebx-eb)/(eb-1)*delpsi
      go to 50
   41 pf = .5_R8*p0*(1._R8-phiarg)**2
      ppf = -p0*delpsi*(1._R8-phiarg)
   50 continue
      fsq = gzero**2 + 2._R8*(1._R8/betaj-1)*xplas**2*pf
!
        if (fsq.lt.0)   then
                        if (fkgrat.lt.0.001_R8)   go to 100
                        fsq = (fkgrat*gzero)**2
                        ggp = (fsq-gzero**2) * ppf/pf
        else
        ggp = (1._R8/betaj-1._R8) * xplas**2 * ppf
                        endif
        gval = sqrt(fsq)
        gpval= ggp / gval
        if (gval/gzero .lt. fgrat)   go to 200
        return
!
!
!.....ifunc=3  rfp profiles (LANL)
    3 continue
      if(psv.ge.psilim) go to 40
      phiarg = (psilim-psv)*delpsi
      if(delg.le.0._R8) go to 60
      gvalt = gzero + gprfp*phiarg**alphag
      gpval = -delpsi*gprfp*alphag*phiarg**(alphag-1._R8)
      gval = gvalt
      ggpt = gvalt*gpval
      ggp = ggpt
      return
 60   continue
      gvalt = gprfp*(1._R8+(delg-1._R8)*phiarg**alphag)
      gpval = -delpsi*alphag*gprfp*(delg-1._R8)*phiarg**(alphag-1._R8)
      gval = gvalt
      ggpt = gvalt*gpval
      ggp = ggpt
      return
!
!.....ifunc=4  spheromak profiles (Princeton)
    4 continue
      ggpt = gp1*ff(3,psv,ig,jg)
      ggp = ggpt
      temp = gzero**2 + 2._R8*gp1*ff(1,psv,ig,jg)
      if(temp.lt.0) go to 100
      gvalt = sqrt(temp)
      gval = gvalt
      if(gvalt.eq.0) go to 100
      gpval = ggpt/gvalt
      return
!
!.....ifunc=5   Ohmic Profiles
    5 continue
      do 96 isfa = 2,npsit
      ip = isfa
      if(psv.lt.xsv2(isfa)) go to 97
   96 continue
   97 continue
!
      iz = ip-1
      denom = xsv2(ip) - xsv2(iz)
      facp = (psv-xsv2(iz))/denom
      facz = 1._R8- facp
      vpl = tpi*(facz*qprof2(iz)*vp2(iz) + facp*qprof2(ip)*vp2(ip))
      etal =     facz*etpara(iz)        + facp*etpara(ip)
      gxmjal =   facz*gxmja2(iz)         + facp*gxmja2(ip)
      xmjal  =   facz*xmja2(iz)          + facp*xmja2(ip)
      denom = (1._R8+ gxmjal/(xmjal*gs(ig,jg)**2))
      call peval(psv,1,pval,ppval,ig,jg)
      term1 =  -ppval*vpl/(xmjal*denom)
      if(ifunc .eq. 5) then
      term2 =  -gp1/(tpi*etal*denom)
      endif
      if(ifunc .eq. 7) then
      term2 = -gp1*rjdb(ig,jg)*vpl/(xmjal*denom)
      endif
!
      ggp = term1 + term2 - rjcd(ig,jg)*vpl/(xmjal*denom)
   10 continue
!
!.....get g from interpolating polynomial for surface averaged transport
      lval = 2
      do 20 lv=3,npsit
      lval = lv
      if(xsv(lval).gt.psv) go to 30
   20 continue
   30 continue
      lm = lval-1
      gvalt=asv(lm,1)+psv*(bsv(lm,1)+psv*(csv(lm,1)+psv*dsv(lm,1)))
      gval = gvalt
      gpvalt = bsv(lm,1)+psv*(2._R8*csv(lm,1)+3._R8*dsv(lm,1)*psv)
      gpval = gpvalt
      if(itype.eq.1 .and. ifunc.eq.5) return
      if(itype.eq.1 .and. ifunc.eq.7) return
      ggp = gvalt*gpvalt
      return
!
   40 continue
!
!.....vacuum region
      ggp = 0._R8
      gval = gzero
      gpval = 0._R8
      return
  100 continue
!.....error exit
      ineg=13
      ggp = 0
      gval = gzero
      gpval = 0._R8
      write(nout,1000) psv,psimin,psilim
 1000 format(" error in ggp,psv,psimin,psilim=",1p3e12.4)
!
  200   icallg = icallg + 1
        if (mod(icallg,20).eq.0)  write (nout, 220)
  220  format(/'   delpsi   psimin      psv   phiarg   alphag       p0   &  
     &     pf      ppf    itype       ig       jg      fsq     gval      &  
     & ggp')
        write (nout, 240)  delpsi, psimin, psv,phiarg, alphag, p0, pf,   &  
     &                     ppf, itype, ig, jg, fsq, gval, ggp
  240   format (1p8e9.2, 3i9, 1p3e9.2)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
