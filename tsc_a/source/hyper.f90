      subroutine hyper
!
!
!.....add hyperresistivity contributions to g and psi
!.....by solving implicit  equations
!
      USE CLINAM
      USE SCR1
      USE SCR5
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j,nloopm,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 df,epshyp,alam,alamf,xa,xc,x5,psiaa,psitw,psicc,psidd
      REAL*8 gpx,gpz,xbsq,v28t,denomx,alpha,denom,ajp2ij,ggx,ggz
      REAL*8 gsij,t1,t2,gzx5sq,x6,x4
!===========
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: etaa,etab,alama,            &  
     &                                       alamb,alamc,fac3p,fac4p,    &  
     &                                       fac41p,fac42p,fac43p,       &  
     &                                       fac44p,alamd,termg,psios,   &  
     &                                       pso2s,psave,gos,go2s,       &  
     &                                       gsave,trm1i,trm2i,trm3i
!     REAL*8, ALLOCATABLE, DIMENSION(:)   :: t6y,t4y,t8y,v4tv,v5tv,v6tv
      INTEGER :: istat

      if(.NOT.ALLOCATED(etaa)) ALLOCATE(etaa(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(etab)) ALLOCATE(etab(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(alama)) ALLOCATE(alama(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(alamb)) ALLOCATE(alamb(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(alamc)) ALLOCATE(alamc(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(alamd)) ALLOCATE(alamd(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(termg)) ALLOCATE(termg(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(psios)) ALLOCATE(psios(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(pso2s)) ALLOCATE(pso2s(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(psave)) ALLOCATE(psave(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(gos)) ALLOCATE(gos(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(go2s)) ALLOCATE(go2s(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(gsave)) ALLOCATE(gsave(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(trm1i)) ALLOCATE(trm1i(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(trm2i)) ALLOCATE(trm2i(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(trm3i)) ALLOCATE(trm3i(pnx,pnz), STAT=istat)
!     ALLOCATE (
!    &          etaa(pnx,pnz),etab(pnx,pnz),
!    &          alama(pnx,pnz),alamb(pnx,pnz),
!    &          alamc(pnx,pnz),alamd(pnx,pnz),
!    &          fac3p (pnx,pnz),fac4p (pnx,pnz),
!    &          fac41p(pnx,pnz),fac42p(pnx,pnz),fac43p(pnx,pnz),
!    &          fac44p(pnx,pnz),a
!    &          termg(pnx,pnz),
!    &          psios(pnx,pnz),pso2s(pnx,pnz),psave(pnx,pnz),
!    &          gos(pnx,pnz),go2s(pnx,pnz),gsave(pnx,pnz),
!    &          trm1i(pnx,pnz),trm2i(pnx,pnz),trm3i(pnx,pnz) ,
!    &          t6y(pnz),t4y(pnz),t8y(pnz),
!    &          v4tv(pnx),v5tv(pnx),v6tv(pnx),
!    &          STAT=istat )
      if (istat .ne. 0) stop "Allocation Error : hyper"
      

      do 50 i=2,nxp
      do 50 j=2,nzp
      etaa(i,j) = .5_R8*(etay(i,j)+etay(i,j-1))*dxisq/xary(i)
      etab(i,j) = .5_R8*(etay(i,j)+etay(i-1,j))*dzisq
   50 continue
!
!
      if(dt+dtold .le. 0) then
!     DEALLOCATE (etaa,etab,alama,alamb,alamc,alamd,termg,psios,
!    &            pso2s,psave,gos,go2s,gsave,trm1i,trm2i,trm3i
!    &            ,t6y,t4y,t8y,v4tv,v5tv,v6tv
!    &            )
      return
      end if

      df = acoef(67)
      sf = acoef(68)
       epshyp = acoef(59)
       nloopm = acoef(60)
      alam = alamf(psimin)
!
      do 417 i=iminn-1,imaxx
      xa = .5_R8*(xary(i)+xary(i+1))
      xc = .5_R8*(xary(i)+xary(i-1))
      x5 = xary(i)
      do 417 j=jminn-1,jmaxx
      psiaa=.5_R8*(psi(i+1,j)+psi(i,j))
      psitw=.5_R8*(psi(i,j+1)+psi(i,j))
      psicc=.5_R8*(psi(i-1,j)+psi(i,j))
      psidd=.5_R8*(psi(i,j-1)+psi(i,j))
      alama(i,j) = xa*alamf(psiaa)/(x5*deex**2)
      alamb(i,j) =    alamf(psitw)/    deez**2
      alamc(i,j) = xc*alamf(psicc)/(x5*deex**2)
      alamd(i,j) =    alamf(psidd)/    deez**2
!
  417 continue
!
!
!
!
      do 415 i=iminn,imaxx
      do 415 j=jminn,jmaxx
      gpx = (psi(i+1,j)-psi(i-1,j))/(2._R8*deex)
      gpz = (psi(i,j+1)-psi(i,j-1))/(2._R8*deez)
      xbsq = gpx**2 + gpz**2 + gs(i,j)**2
      bmagy(i,j) = xbsq/xary(i)**2
!
  415 continue
!
!
!
      do 608 i=iminn,imaxx
      v28t = dzisq
      v4tv(i) = xary(i)*dxisq/xarh(i)
      v6tv(i) = xary(i)*dxisq/xarh(i+1)
      v5tv(i) = -(v4tv(i)+v6tv(i)+2._R8*v28t)
      do 608 j=jminn,jmaxx
      psave(i,j) = psi(i,j)
      psios(i,j) = psi(i,j)
      gsave(i,j) = g(i,j)
      gos(i,j) = g(i,j)
      denomx =  (-(dt+dtold)*(alama(i,j)+alamb(i,j)                      &  
     &+.5_R8*(alam)*(1._R8/deez**2+1._R8/deex**2)+alamc(i,j)+alamd(i,j))  &  
     & *v5tv(i))
      trm1i(i,j) = 0.0_R8
      trm2i(i,j) = 0.0_R8
      trm3i(i,j) = 1.0_R8
      if(denomx.eq.0) go to 608
      alpha = sf*bmagy(i,j)/denomx
      denom = (1._R8+df/nx+alpha)
      trm1i(i,j) = 2._R8/denom
      trm2i(i,j) = (1._R8-df/nx)/denom
      trm3i(i,j) = alpha/denom
  608 continue
!
      nloop = 0
  609 continue
      do 610 i=iminn-1,imaxx+1
      do 610 j=jminn-1,jmaxx+1
      pso2s(i,j) = psios(i,j)
      psios(i,j) = psi(i,j)
      go2s(i,j) = gos(i,j)
      gos(i,j) = g(i,j)
  610 continue
!
  611 do 615 i=iminn,imaxx
      do 615 j=jminn,jmaxx
      ajp2ij = dzisq*(psi(i,j-1)+psi(i,j+1))                             &  
     &          + v4tv(i)*psi(i-1,j) + v6tv(i)*psi(i+1,j)                &  
     &          + v5tv(i)*psi(i,j)
      ajp2(i,j) = ajp2ij
      gpx = (psi(i+1,j)-psi(i-1,j))/(2._R8*deex)
      gpz = (psi(i,j+1)-psi(i,j-1))/(2._R8*deez)
      ggx = (xsqoj(i+1)*(g(i+1,j)+g(i+1,j+1))                            &  
     &      -xsqoj(i)*(g(i,j)+g(i,j+1)))/(2._R8*deex)
      ggz = (xsqoj(i)*(g(i,j+1)-g(i,j))                                  &  
     &      +xsqoj(i+1)*(g(i+1,j+1)-g(i+1,j)))/(2._R8*deez)
      gsij = .25_R8*(xsqoj(i)*(g(i,j+1)+g(i,j))                          &  
     &          + xsqoj(i+1)*(g(i+1,j+1)+g(i+1,j)))
      xbsq = gpx**2 + gpz**2 + gsij**2
      ajp3(i,j) = (ajp2ij*gsij                                           &  
     &           -(gpx*ggx+gpz*ggz))/xbsq
!
  615 continue
      if(isym.ne.1) go to 622
      do 621 i=iminn,imaxx
  621 ajp3(i,1) = ajp3(i,3)
  622 continue
!
      delpmx = 0.0_R8
      delgmx = 0._R8
      do 616 i=iminn,imaxx
      do 616 j=jminn,jmaxx
!
      t1 = (alama(i,j)*(ajp3(i+1,j)-ajp3(i,j))                           &  
     &     -alamc(i,j)*(ajp3(i,j)-ajp3(i-1,j)))
      t2 = (alamb(i,j)*(ajp3(i,j+1)-ajp3(i,j))                           &  
     &     -alamd(i,j)*(ajp3(i,j)-ajp3(i,j-1)))
      ajp4(i,j) = (t1+t2)/bmagy(i,j)
  616 continue
!
      do 617 i=iminn,imaxx
      do 107 j = 1,nzp
      termg(i,j) = 0.0_R8
  107 continue
      t8y(2) = 0._R8
      t8y(jmaxx+1) = 0._R8
      do 105 j=jminn+1,jmaxx+1
  105 t4y(j) = t6y(j)
      do 101 j=3-isym,jmaxx
!
      t8y(j)=(psios(i,j+1)+psios(i-1,j+1)-psios(i,j-1)-psios(i-1,j-1))   &  
     &  *.25_R8*(ajp4(i,j)+ajp4(i-1,j))/(xary(i)+xary(i-1))*(deex/deez)
  101 continue
  103 continue
      do 102 j=jminn+1,jmaxx+1
      t6y(j)=(deez/deex)*.125_R8*(psios(i+1,j)+psios(i+1,j-1)-psios(i-1,  &  
     & j)                                                                &  
     &     -psios(i-1,j-1))*(ajp4(i,j)+ajp4(i,j-1))/xary(i)
  102 continue
      do 108 j = jminn+1,jmaxx+1
      termg(i,j) = t8y(j) - t8y(j-1) + t6y(j) - t4y(j)
  108 continue
  104 continue
      if(alam .le. 0) go to 620
      do 618 j=jminn,jmaxx
!
      psi(i,j) = psios(i,j)*trm1i(i,j) - trm2i(i,j)*pso2s(i,j)           &  
     &         + trm3i(i,j)*(-(dt+dtold)*ajp4(i,j)*gs(i,j) + psave(i,j))    
      delpmx = max(delpmx,abs(psi(i,j)-psios(i,j))                       &  
     &           /((psilim-psimin)))
  618 continue
      if(i.lt.3) go to 623
      denom = abs(g(imag,jmag))
      do 619 j=jminn+1,jmaxx
!
      g(i,j) = gos(i,j)*trm1i(i,j) - trm2i(i,j)*go2s(i,j)                &  
     &         + trm3i(i,j)*( (dt+dtold)*termg(i,j) + gsave(i,j))
      delgmx = max(delgmx,abs(g(i,j)-gos(i,j))                           &  
     &           /denom)
  619 continue
      if(isym.eq.1) g(i,2) = g(i,3)
  623 continue
      if(isym.eq.1) psi(i,1) = psi(i,3)
  617 continue
      nloop = nloop+1
      if(nloop.gt.nloopm) go to 999
      if(delpmx .gt. epshyp .or. delgmx .gt. epshyp ) go to 609
  620 continue
!
!.....special arrays for lings plots
!.....December 1987   KML
      if(irfp.ne.1)go to 601
      do 589 i=2,nxp
      x5 = xarh(i)
      gzx5sq = gzero/xarh(i)**2
      x6 = xarh(i+1)
      x4 = xarh(i-1)
      do 588 j=2,nzp
      pling1(i,j) = etay(i,j)*udsr
      pling4(i,j) = ajp3(i,j)
      pling2(i,j) = 0._R8
      pling3(i,j) = 0._R8
      pling5(i,j) = ptsave(i,j)
      if(ajp2(i,j) .le. 0) go to 577
      pling2(i,j) = -gs(i,j)*ajp4(i,j)/ajp2(i,j)*udsr
      pling3(i,j)=(etay(i,j)*ajp2(i,j)-gs(i,j)*ajp4(i,j))                &  
     &            /(ajp2(i,j)*etay(i,j))
  577 continue
      gling1(i,j) = gtsave(i,j) - gzx5sq*uo(i,j)
      gling2(i,j) = etaa(i,j)*x6*go(i+1,j) + etaa(i-1,j)*x4*go(i-1,j)    &  
     &            + etab(i,j)*go(i,j+1) + etab(i,j-1)*go(i,j-1)          &  
     &     -((etaa(i,j)+etaa(i-1,j))*x5+etab(i,j)+etab(i,j-1))*go(i,j)
      gling3(i,j) = termg(i,j)
!
  588 continue
      if(isym.eq.0) go to 589
      gling1(i,2) = gling1(i,3)
      gling2(i,2) = gling2(i,3)
      gling3(i,2) = gling3(i,3)
  589 continue
!
      do 599 ii=1,nwire
      i = iwire(ii)
      j = jwire(ii)
      pling1(i,j) = etav*udsr
  599 pling3(i,j) = 0._R8
  601 continue
!
      do 510 i=iminn,imaxx
      do 510 j=jminn,jmaxx
      ptsave(i,j) = ptsave(i,j) + (psi(i,j)-psave(i,j))/(dt+dtold)
      psi(i,j) = psave(i,j)
      gtsave(i,j) = gtsave(i,j) + (g(i,j)-gsave(i,j))/(dt+dtold)
      g(i,j) = gsave(i,j)
  510 continue
!
!     DEALLOCATE (etaa,etab,alama,alamb,alamc,alamd,termg,psios,
!    &            pso2s,psave,gos,go2s,gsave,trm1i,trm2i,trm3i
!    &            ,t6y,t4y,t8y,v4tv,v5tv,v6tv
!    &            )
      return
  999 continue
      ineg=23
      write(nterm,6666) nloop, delpmx, delgmx, epshyp
 6666 format(" nloop, delpmx, delgmx, epshyp =", i5,1p3e12.4)
!     DEALLOCATE (etaa,etab,alama,alamb,alamc,alamd,termg,psios,
!    &            pso2s,psave,gos,go2s,gsave,trm1i,trm2i,trm3i
!    &            ,t6y,t4y,t8y,v4tv,v5tv,v6tv
!    &            )

      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
