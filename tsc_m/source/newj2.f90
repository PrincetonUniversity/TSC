      subroutine newj2(iter)
!****
!**** using flux from array old, this routine calculates the new
!**** right hand side of the grad-shafranov equation and loads it
!**** into the interior of the array psi
!**** right hand side=-tpi**2*(x*x*pp(old)  +
!****    + gp1*ff(3,old) + gp2*ff(4,old) )
!**** gp1 is calculated to fix the on axis q
!**** gp2 is calculated to fix the total toroidal plasma current
!**** qzero is the input q on axis
!****
      USE CLINAM
      USE SAPROP
      USE SVDCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iter,idbug,j,i,kk,k,isfa,ip,iz,ii,ig,l,ie,n,np,iem10
      INTEGER iem20,ngr,iabs,ngrvwl,ngrvcl,nwirem
      INTEGER INDEX
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rfac1,rfac2,rf,rfm,gpval,gppval,finterp,ajbs,ajcd,ajfw
      REAL*8 ajlh,ajec,sump,sum1,sum2,sum3,fac,psimid,xmid,ajmid,pval
      REAL*8 ppval,ff,denom,facp,facz,vpl,etal,gxmjal,xmjal,c1,c2
      REAL*8 c3,c4,c5,c6,c7,c0,aquad,bquad,cquad,rad,temp,root1
      REAL*8 root2,gval,aplfb,f1,psinn,pinterp,psinp,term,ellip
      REAL*8 delta,x1,x2,z1,z2,fluxdif,psi1,psi2,psi3,psi4,psisn1
      REAL*8 rfg,rfgm,anum,adenom,gcurtot,factor
!============
      data idbug/0/
      data rfac1/1.0_R8/
      data rfac2/0.0_R8/
      if(lrswtch.eq.0) go to 109
      psimin = -100._R8
      psilim = -99._R8
      delpsi = 1._R8
      go to 29
  109 continue
      rf = acoef(50)
      rfm = 1._R8-rf
!
      do 10 j=1,nzp
      do 10 i=1,nxp
      if(ifrst(7).eq.0) go to 11
      psi(i,j) = rf*psi(i,j) +rfm*psio(i,j)
   11 psio(i,j) = psi(i,j)
   10 continue
      ifrst(7) = 1
!
      call limpsi
      call magaxis
      call delpdef
      if(igone.ge.1) go to 29
!
      if(ifunc.ne.5.and.ifunc.ne.7) go to 15
      call flxvol(3)
      if(ineg.ne.0) return
      call defcoef
      call sprop
      call trcdef
!     call auxheat
      call defjdb
      call curdrive
      do 14 i=iminn,imaxx
      do 14 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 14
      call geval(psi(i,j),2,gs(i,j),gpval,gppval,i,j)
      rjcd(i,j) = 0._R8
      rjdb(i,j)=0._R8
!
      if(psi(i,j) .gt. psilim) go to 14
      do 640 kk=2,npsit
      k = kk
  640 if(psi(i,j) .lt. xsv2(kk)) go to 630
      k = npsit
!
  630 finterp = (psi(i,j)-xsv2(k-1))/(xsv2(k)-xsv2(k-1))
!
      ajbs = (ajavbs(k-1)+finterp*(ajavbs(k)-ajavbs(k-1)))*gs(i,j)
      ajcd = (ajavcd(k-1)+finterp*(ajavcd(k)-ajavcd(k-1)))*gs(i,j)
      ajfw = (ajavfw(k-1)+finterp*(ajavfw(k)-ajavfw(k-1)))*gs(i,j)
      ajlh = (ajavlh(k-1)+finterp*(ajavlh(k)-ajavlh(k-1)))*gs(i,j)
      ajec = (ajavec(k-1)+finterp*(ajavec(k)-ajavec(k-1)))*gs(i,j)
      rjcd(i,j) = ( ajcd + ajfw + ajbs + ajlh + ajec)
!
      rjdb(i,j)=(ajaveq(k-1)+finterp*(ajaveq(k)-ajaveq(k-1)))*gs(i,j)
   14 continue
   15 continue
!
      if(isvd.ge.1) call xptcalc
      if(psilim.le.psimin) ineg=35
      if(ineg.ne.0) return
!
      psilimo = psilim
!****
      sump=0.0_R8
      sum1=0.0_R8
      sum2=0.0_R8
      sum3 = 0._R8
      do 102 i=iminn,imaxx
      do 101 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 101
      fac=1.0_R8
      if(isym.eq.1 .and. j.ne.2) fac = 2._R8
      psimid = psi(i,j)
      if(psimid.ge.psilim) go to 101
      xmid = xary(i)
      ajmid = xmid*dxdz
      call peval(psimid,1,pval,ppval,i,j)
      sump = sump - ppval*fac*ajmid
!
      go to(91,92,93,94,95,91,95),ifunc
!
!.....ifunc=1 Tokamak Profiles (Princeton)
   91 continue
      sum1 = sum1 - (ff(3,psimid,i,j)/xmid**2)*fac*ajmid
      sum2 = sum2 - (ff(4,psimid,i,j)/xmid**2)*fac*ajmid
      go to 101
!
!.....ifunc=2  Tokamak Profiles (ORNL)
   92 sum1 = sum1 - (ff(5,psimid,i,j)/xmid**2)*fac*ajmid
      go to 101
!
!.....ifunc=3  RFP Profiles (LANL)
 93   continue
      if(delg.le.0._R8) go to 600
      sum1 = sum1 + (ff(7,psimid,i,j)/xmid**2)*fac*ajmid
      sum2 = sum2 + (ff(8,psimid,i,j)/xmid**2)*fac*ajmid
      go to 101
 600  continue
      sum1 = sum1 - (ff(6,psimid,i,j)/xmid**2)*fac*ajmid
      go to 101
!
!.....ifunc=4  Spheromak Profiles (Princeton)
   94 sum1 = sum1 - (ff(3,psimid,i,j)/xmid**2)*fac*ajmid
      go to 101
!
!.....ifunc=5   Ohmic Profiles
   95 continue
      do 96 isfa=2,npsit
      ip = isfa
      if(psimid.lt.xsv2(isfa)) go to 97
   96 continue
   97 continue
      iz = ip - 1
      denom = xsv2(ip) - xsv2(iz)
      facp = (psimid-xsv2(iz))/denom
      facz = 1._R8- facp
      vpl = tpi*(facz*qprof2(iz)*vp2(iz)+facp*qprof2(ip)*vp2(ip))
      etal =     facz*etpara(iz)        +facp*etpara(ip)
      gxmjal =   facz*gxmja2(iz)        +facp*gxmja2(ip)
      xmjal =    facz*xmja2(iz)         +facp*xmja2(ip)
      denom = (1._R8+ gxmjal/(xmjal*gs(i,j)**2))
      sum1 = sum1 + ppval*fac*ajmid*(-1._R8+ vpl/(xmjal*xmid**2*denom))
      sum3 = sum3 + fac*ajmid*rjcd(i,j)*vpl/(xmjal*xmid**2*denom)
      if(ifunc .eq. 5) then
      sum2 = sum2 + fac*ajmid/(tpi*etal*xmid**2*denom)
      endif
      if(ifunc .eq. 7) then
      sum2 = sum2 + fac*ajmid*rjdb(i,j)*vpl/(xmjal*xmid**2*denom)
      endif
!     if(iter.eq.1.and. i.eq.10 .and. j.eq.91) then
!     write(nterm,8811)i,j,iz,ip,npsit,xsv2(39),xsv2(40),xsv2(41),psimid
!8811 format(5i5,1p4e15.7)
!     endif
!
  101 continue
  102 continue
!
      go to(81,82,83,84,85,81,85),ifunc
!
!.....ifunc=1  Tokamak Profiles (Princeton)
   81 continue
!**** calculate gp1 and gp2 to fix total current and qzero
!**** define some constants
      if(abs(sum2).lt.1.E-12_R8)go to 1001
      call peval(psimin,1,pval,ppval,imag,jmag)
      c1=(tcuro-sump)/sum2
      c2=sum1/sum2
      c3=2.0_R8*ff(1,psimin,imag,jmag)
      c4=2.0_R8*ff(2,psimin,imag,jmag)
      c5=xmag*xmag*ppval
      c6=ff(3,psimin,imag,jmag)
      c7=ff(4,psimin,imag,jmag)
      c0=pi*xmag*qzero /tpi
!**** solve the quadratic equation for gp1
!**** define the coefficients
      aquad=c0*c0*(c6-c7*c2)**2
      bquad=2.0_R8*c0*c0*(c5+c7*c1)*(c6-c7*c2)+c4*c2-c3
      cquad=c0*c0*(c5+c7*c1)**2-c4*c1-gzero**2
!**** check for imaginary roots
      rad=bquad*bquad-4.0_R8*aquad*cquad
      if(rad.lt.0.0_R8)go to 1002
      temp=sqrt(rad)
      if(abs(aquad).lt.1.E-12_R8)go to 1003
      root1=0.5_R8*(temp-bquad)/aquad
      root2=-0.5_R8*(temp+bquad)/aquad
!**** set gp1=root2 only
      gp1=rfac1*root1+rfac2*root2
      gp2=c1-c2*gp1
      go to 1005
!
!**** errors
 1001 continue
      write(nout,3001)sump,sum1,sum2
 3001 format(" newj2 error sum2 too small=",1p3e14.6)
      ineg=10
      return
 1002 continue
      write(nout,3002)
 3002 format(" newj2 error imaginary roots")
      ineg=10
      return
 1003 continue
      write(nout,3003)
 3003 format(" newj2 error aquad too small")
      ineg=10
      return
 1005 continue
!**** check q0
!
      call geval(psimin,1,gval,gpval,gppval,imag,jmag)
      q0=- 2._R8*gval /(xmag*(c5+gp1*c6+gp2*c7))
      go to 29
!
!.....ifunc=2 Tokamak Profiles (ORNL)
   82 continue
      gp1 = tcuro/(sump+sum1)
      p0 = p0*gp1
      e0 = e0*gp1
      do ii=1,ntpts
      ppres(ii) = ppres(ii)*gp1
      enddo
      gp2 = 0._R8
      go to 29
!
!.....ifunc=3 RFP Profiles (LANL)
   83 continue
      if(delg.le.0._R8) go to 606
      if(abs(sum2).lt. 1.E-12_R8) go to 1001
      c1=(tcuro-sump)/sum2
      c2=sum1/sum2
      rad=c2*c2+4._R8*c1
      if(rad.lt. 0.0_R8)go to 1002
      gprfp=-0.5_R8*c2+0.5_R8*sqrt(rad)
      gp1=gprfp
      go to 29
 606  continue
      if(sum1.eq.0) go to 1001
      sum2 = -(tcuro-sump)/sum1
      if(sum2.lt.0) go to 1001
      gprfp = -sqrt(sum2)
      gp1 = gprfp
      gzero = gprfp
      go to 29
!
!.....ifunc=4 Spheromak Profiles (Princeton)
   84 if(sum1.eq.0) go to 1001
      gp1 = (tcuro-sump)/sum1
      gp2 = 0
      go to 29
!
!.....ifunc=5   Ohmic profiles
   85 continue
      gp1 = (tcuro-sum3-sum1)/sum2
!
!.....DIAG
!     write(nterm,8810) iter,tcuro,sum1,sum2,sum3
!8810 format(" iter,tcuro,sum1,sum2,sum3",i5,1p4e15.7)
   29 continue
!
!.....initialize current to zero
      do 30 i=2,nx
      do 30 j=2,nz
      ajphi(i,j) = acoef(850)*usdv/(etav*tpi*xary(i))
   30 vd(i,j) = xary(i)*dxdz*ajphi(i,j)
      aplfb = 0._R8
      do 804 ig=1,pngroup
  804 gcurfb(ig) = 0._R8
      if(idata.eq.6) call fedtsc

 8154 format(1p10e12.4)
!
!.....add  feedback current to plasma and coil currents
      if(numfb.le.0) go to 812
      if(acoef(901).gt.0 .and. iter.gt.int(acoef(908))) go to 812
      f1 = 1._R8
!
      do 811 l=1,numfb
!
!     if(tpro(istart) .lt. tfbon(l).or.
!    1   tpro(istart) .gt. tfbof(l) ) go to 811
      if(tfbon(l).ge.0 .and. tpro(istart).lt.tfbon(l)) go to 811
      if(tfbons(l).lt.0 .and. 0 .lt. int(-tfbons(l)) ) go to 811
      if(tfbof(l).ge.0 .and. tpro(istart).gt.tfbof(l)) go to 811
      if(tfbofs(l).lt.0 .and. 0 .gt. int(-tfbofs(l)) ) go to 811
!
      ie = ipext(l)
      if(ie.gt.1000) go to 1816
      if(ie.gt.20)  go to 1820
      if(ie.gt.10) go to 1910
      go to(806,806,806,811,1810,1811,1812,1813,1814,1815),ie
  806 continue
      n = 2*nfeedv(istart,l)-1
      np = n+1
      psinn = pinterp(xobs(n ),zobs(n ),iobs(n ),jobs(n ))
      psinp = pinterp(xobs(np),zobs(np),iobs(np),jobs(np))
      fac = f1
      if(fbfac1(l).gt.0) fac = fac*pcur(istart)/pcur(ntpts)
      term = fbfac(l)*(psinn-psinp-fbcon(l)*fac)
      if(nrfb(l).lt.0) go to 811
      if(nrfb(l).eq.0) go to 820
!..rxw/23/04/87
      go to 1803
 1810 continue
      term = fbfac(l)*(xmag-xmagw - fbcon(l))
      go to 1803
 1811 continue
      term = fbfac(l)*(zmag-zmagw - fbcon(l))
      go to 1803
 1812 term = fbfac(l)*(eps1c-eps10)
      go to 1803
 1813 term = fbfac(l)*(eps2c-eps20)
      go to 1803
 1814 term = fbfac(l)*(eps3c-eps30)
      go to 1803
 1815 term = fbfac(l)*(eps4c-eps40)
      go to 1803
 1910 continue
      call shape(ellip,delta,x1,x2,z1,z2)
      iem10 = ie - 10
      go to(1911,1912,1913,1914),iem10
 1911 fluxdif = fbfac(l)*(.5_R8*(x2+x1) - rzerw)
      go to 1803
 1912 fluxdif = fbfac(l)*(.5_R8*(x2-z1) - azerw)
      go to 1803
 1913 fluxdif = fbfac(l)*(ellip - ezerw)
      go to 1803
 1914 fluxdif = fbfac(l)*(delta - dzerw)
      go to 1803
 1820 continue
      call fluxmod(psi1,psi2,psi3,psi4,psisn1,2.40_R8,0.40_R8)
      iem20 = ie - 20
      go to(1821,1822,1823,1824),iem20
 1821 term = fbfac(l)*psi1
      go to 1803
 1822 term = fbfac(l)*psi2
      go to 1803
 1823 term = fbfac(l)*psi3
      go to 1803
 1824 term = fbfac(l)*psi4
      go to 1803
 1816 continue
      index = ie-1000 + ncoil-nwire
      term = ccoil(index)*udsi
 1803 continue
      ig = nrfb(l)
!....changes suggested by Weiner
      if(ipext(l).lt.5 .or. ipext(l).gt.6) then
      gcurfb(ig) = gcurfb(ig)+term*usdi
      else
      gcurfb(ig) = gcurfb(ig)+term*apl*1.E-6_R8*usdi
      endif
      go to 811
  820 aplfb = aplfb + term*usdi
  811 continue
  812 continue
      do 813 ig=1,pngroup
      rfg    = acoef(51)
      rfgm   = 1._R8- rfg
      gcurfb(ig) = rfgm*gcurfbo(ig)+gcurfb(ig)*rfg
      gcurfbo(ig) = gcurfb(ig)
  813 continue
      if(acoef(901) .gt. 0._R8.and. iter .gt. acoef(907)) then
      call fixshape(iter)
      endif
!
      tcuro = (pcur(istart)*usdi+aplfb)
      fac = sqrt(nx*(nz-1)/((2-isym)*4._R8*(alx-ccon)*alz))
      anum = (vloopv(istart)-acoef(850))*ndiv
      adenom = etav*fac*acoef(11)*udsv*tpi
      if(adenom .ne. 0) tcuro = tcuro + anum/adenom
!
!

      do 830 ii=1,nwire
      ngr = iabs(igroupw(ii))
      cwire0(ii) = 0._R8
!
!.....check for initialization using lrswtch
      if(lrswtch.eq.0) go to 607
      nn = ncoil-nwire+ii
      if(ngr.ne.lrswtch) go to 607
      if(zcoil(nn).gt.0) cwire0(ii) = acoef(12)*usdi/rswires(ii)
      if(zcoil(nn).lt.0) cwire0(ii) =-acoef(12)*usdi/rswires(ii)
  607 continue
      cwire0(ii) = cwire0(ii) + cwics(ii)*usdi
      cwire0(ii) = cwire0(ii)                                            &  
     &           +  (gcur(istart,ngr)*usdi+gcurfb(ngr))*aturnsw(ii)
!
      ngrvwl = ngrvw1(ii)
      if(ngrvwl.le.0) go to 829
      gcurtot = gcurfb(ngrvwl) + gcur(istart,ngrvwl)*usdi
      cwire0(ii) = cwire0(ii)+atnvw1(istart,ii)*gcurtot
!
  829 continue
      ngrvwl = ngrvw2(ii)
      if(ngrvwl.le.0) go to 839
      gcurtot = gcurfb(ngrvwl) + gcur(istart,ngrvwl)*usdi
      cwire0(ii) = cwire0(ii)+atnvw2(istart,ii)*gcurtot
!
  839 continue
      ngrvwl = ngrvw3(ii)
      if(ngrvwl.le.0) go to 849
      gcurtot = gcurfb(ngrvwl) + gcur(istart,ngrvwl)*usdi
      cwire0(ii) = cwire0(ii)+atnvw3(istart,ii)*gcurtot
!
  849 continue
      ngrvwl = ngrvw4(ii)
      if(ngrvwl.le.0) go to 859
      gcurtot = gcurfb(ngrvwl) + gcur(istart,ngrvwl)*usdi
      cwire0(ii) = cwire0(ii)+atnvw4(istart,ii)*gcurtot
     
!
  859 continue
!
!.....add initial current due to loop voltage
      if(iseries(ngr) .ne. 0) go to 858
      if(ilwire(ii).eq.1) cwire0(ii)=cwire0(ii)                          &  
     &                    +(vloopp + acoef(850))*usdv/rswire(ii)
  858 continue
!..rxw/05/10/87
      factor=1._R8
      if( (acoef(290).eq.2) .and. (save1(ii).gt.9._R8)                   &  
     &                      .and. (vgain(ii).ne.0)                       &  
     &                      .and. (idata.ne.7)                           &  
     &                      .and. (save2(ii).le.times) )                 &  
     &  factor=vgain(ii)/( vgain(ii) + rswires(ii) )
!===      ccoil(ncoil-nwire+ii) = cwire0(ii)
      ccoil(ncoil-nwire+ii) = factor * cwire0(ii)
      ccoils(ncoil-nwire+ii)= ccoil(ncoil-nwire+ii)*udsi

!...rxw/end
  830 continue
!
!.....recompute external coils
      if(ncoil.eq.nwire) go to 890
      do 896 ii=1,ncoil-nwire
      ngr = iabs(igroupc(ii))
      ccoil(ii) = aturnsc(ii)*(gcur(istart,ngr)*usdi + gcurfb(ngr))      &  
     &          + ccics(ii)*usdi
!
      ngrvcl = ngrvc1(ii)
      if(ngrvcl.le.0) go to 929
      gcurtot = gcurfb(ngrvcl) + gcur(istart,ngrvcl)*usdi
      ccoil(ii) = ccoil(ii)+atnvc1(istart,ii)*gcurtot
!
  929 continue
      ngrvcl = ngrvc2(ii)
      if(ngrvcl.le.0) go to 939
      gcurtot = gcurfb(ngrvcl) + gcur(istart,ngrvcl)*usdi
      ccoil(ii) = ccoil(ii)+atnvc2(istart,ii)*gcurtot
!
  939 continue
      ngrvcl = ngrvc3(ii)
      if(ngrvcl.le.0) go to 949
      gcurtot = gcurfb(ngrvcl) + gcur(istart,ngrvcl)*usdi
      ccoil(ii) = ccoil(ii)+atnvc3(istart,ii)*gcurtot
!
  949 continue
      ngrvcl = ngrvc4(ii)
      if(ngrvcl.le.0) go to 959
      gcurtot = gcurfb(ngrvcl) + gcur(istart,ngrvcl)*usdi
      ccoil(ii) = ccoil(ii)+atnvc4(istart,ii)*gcurtot
!
  959 continue
      ccoils(ii) = ccoil(ii)*udsi
  896 continue
  890 continue
!
!.....add currents due to wires in grid,  do not include plasma current
      if(nwire.eq.0) go to 151
      nwirem = nwire
!cj dir$ ivdep
      do 150 ii=1,nwirem
      n = ncoil-nwire+ii
      i = iwire(ii)
      j = jwire(ii)
      vd(i,j) = ccoil(n)*xary(i)
      ajphi(i,j) = ccoil(n)/dxdz
!
  150 continue
  151 continue
      if(lrswtch.ne.0.or.igone.ge.1) return
      if(ifunc.eq.5.or.ifunc.eq.7) gp2 = 0._R8
!****
!**** load rhs of equation
      do 500 i=iminn,imaxx
      do 400 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 400
      psimid = psi(i,j)
      xmid = xary(i)
      if(psimid.gt.psilim) go to 400
      if(psimid.lt.psimin) go to 400
      call geval(psimid,1,gval,gpval,gppval,i,j)
      call peval(psimid,1,pval,ppval,i,j)
      ajphi(i,j) = -(xmid**2*ppval+gppval)/xmid
      vd(i,j) =  xmid*dxdz*ajphi(i,j)
      if(ifunc.ne.5.and.ifunc.ne.7) go to 400
      fac = 1._R8
      if(isym.eq.1 .and. j.ne.2) fac = 2._R8
      gp2 = gp2 - (xmid**2*ppval + gppval)/xmid * dxdz  * fac*udsi
!
  400 continue
  500 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
