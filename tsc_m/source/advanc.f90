      subroutine advanc
!.....* * * number 2.20 * * *                                          c
!
!......calculates time derivative terms for cell centered quantities
!
!.....u = j*delsq(omega)  . . . div of momentum density
!.....w = j*w  . . . toroidal momentum density
!.....r = j*n  . . . mass
!.....g = (j/x**2)*g  . . . toroidal flux
!.....q = j*p**(3/5)  . . . entropy
!
      USE CLINAM
      USE SCADVAN
      USE SCR1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dxdzi,dxodz,rfac,x4,x5,x6,xa,aj4,aj5,aj6
      REAL*8 aja,aj5i,ajai,tw1,tw2,tg1,taw1,tag1,gjsd,gisa,gjph
      REAL*8 omegj,gzjph,gzmid,xmidsq,wjph,rjph,qjph,ajmid,gmid
      REAL*8 wiph,giph,riph,qiph,gziph,omegi,xc,ajc,gisc,omeg1
      REAL*8 omeg2,omegcent,omegdif,denom1,fac5,fac6,denom,fac9
      REAL*8 fac8,fac,terma,termb,termc,termd,fac1,fac2,wt
!============
!     dimension aboune(pnz),abounf(pnz)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aboune
      REAL*8, ALLOCATABLE, DIMENSION(:) :: abounf
!============      
      IF(.not.ALLOCATED(aboune)) ALLOCATE( aboune(pnz), STAT=istat)
      IF(.not.ALLOCATED(abounf)) ALLOCATE( abounf(pnz), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : advanc  ' 
!============      
!
!
!
!...initialize boundary flags : a = right  boundary
!....*** 0 =    boundary ***  : b = top    boundary
!....*** 1 = no boundary ***  : c = left   boundary
!                             : d = bottom boundary
!
      i    = icol
      do 90 j=3,nzp
      abouna(j) = dxa(i,j)/dxisq
      abounb(j) = dzb(i,j)/dzisq
      abounc(j) = dxc(i,j)/dxisq
      abound(j) = dzd(i,j)/dzisq
      abounf(j) = face(i,j)
   90 continue
      do 92 j=3,nzp
   92 aboune(j) = abouna(j)*abounb(j)*abounc(j)*abound(j)
!
!...set grid point being computed
!
!
!...useful differential quantities
!
      dxdzi = 1._R8/dxdz
      dxodz = deex/deez
      rfac  = (deez/deex)**2
!
!...coordinate functions
!
      x4 = xarh(i-1)
      x5 = xarh(i)
      x6 = xarh(i+1)
      xa = xary(i)
      x4 = abs(x4)
      aj4 = x4*dxdz
      aj5 = x5*dxdz
      aj6 = x6*dxdz
      aja = xa*dxdz
      aj5i  = 1._R8/aj5
      ajai  = 1._R8/aja
      tw1   = aj5i*x5**2
      tw2   = aj5i
      tg1   = aj5i/x5**2
      taw1   = ajai*xa**2
      tag1   = ajai/xa**2
      gjsd = deex/(deez*x5)
      gisa = deez/(deex*xa)
!
!.......................................................................
!
!...define  left fluxes
!
!.......................................................................
!
      do 8 j=3,nzp
      wf2pa(j) = wf1pa(j)
      gf2pa(j) = gf1pa(j)
      rf2pa(j) = rf1pa(j)
      qf2pa(j) = qf1pa(j)
      uf2pa(j) = uf1pa(j)
    8 continue
!
!.....define initial top fluxes
      gjph  = .5_R8*(giz(2)+giz(3))
      omegj = omeg(i,3)-omeg(i,2)
      gzjph = gzero/xsqoj(i)
      wf1ta(2) = -.5_R8*(giz(2)+giz(3))*(piz(2)-pim(2))*tw1
      gf1ta(2) =(-.5_R8*(wiz(2)+wiz(3))*(piz(2)-pim(2))*tg1              &  
     &+(-(gjph-gzjph)*omegj*dzisq+gjph*(aiz(2)-aim(2))*aj5i))
!
      gzmid = gzero
      xmidsq = x5**2
      do 11 j=3,nzp
!
!
!.......................................................................
!
!...calculate the top fluxes
!
!.......................................................................
!
      wjph = .5_R8*(wiz(j)+wiz(j+1))
      gjph = .5_R8*(giz(j)+giz(j+1))
      rjph = .5_R8*(riz(j)+riz(j+1))
      qjph = .5_R8*(qiz(j)+qiz(j+1))
      gzjph = gzero/xsqoj(i)
!
      omegj = omeg(i,j+1)-omeg(i,j)
!
      wf1ta(j)= - gjph*(piz(j)-pim(j))*tw1
      gf1ta(j)=(- wjph*(piz(j)-pim(j))*tg1                               &  
     &+(-(gjph-gzjph)*omegj*dzisq+gjph*(aiz(j)-aim(j))*aj5i))
      rf1ta(j)=   rjph*(-omegj*dzisq+(aiz(j)-aim(j))*aj5i)
      qf1ta(j)=   qjph*(-omegj*dzisq+(aiz(j)-aim(j))*aj5i)
      ajmid = .5_R8*(ajphi(i,j)+ajphi(i-1,j))*x5
      gmid = .5_R8*(giz(j+1)+giz(j))*xsqoj(i)
      uf1ta(j) = (ajmid*.25_R8*(piz(j+1)+pim(j+1)-piz(j-1)-pim(j-1))     &  
     &         + (gmid-gzmid)*(giz(j+1)-giz(j))*xsqoj(i)                 &  
     &         + xmidsq*(pr(i,j+1)-pr(i,j)))*gjsd
   11 continue
      do 10 j=3,nzp
!
!.......................................................................
!
!.....calculate the right fluxes
!
!.......................................................................
!
      wiph  = .5_R8*(wip(j)+wiz(j))
      giph  = .5_R8*(gip(j)+giz(j))
      riph  = .5_R8*(rip(j)+riz(j))
      qiph  = .5_R8*(qip(j)+qiz(j))
      gziph = (gzero/(.5_R8*(xsqoj(i)+xsqoj(i+1))))
!
      omegi = omeg(i+1,j)-omeg(i,j)
      wf1pa(j)=   giph*(piz(j)-piz(j-1))*taw1
      gf1pa(j)= ( wiph*(piz(j)-piz(j-1))*tag1                            &  
     &-((giph-gziph)*omegi*dxisq+giph*(aiz(j)-aiz(j-1))*ajai))
      rf1pa(j)= - riph*(omegi*dxisq+(aiz(j)-aiz(j-1))*ajai)
      qf1pa(j)= - qiph*(omegi*dxisq+(aiz(j)-aiz(j-1))*ajai)
      ajmid = .5_R8*(ajphi(i,j)+ajphi(i,j-1))*xa
      gmid = .5_R8*(gip(j)*xsqoj(i+1)+giz(j)*xsqoj(i))
      xmidsq = xa**2
      uf1pa(j) = (ajmid*.25_R8*(pip(j)+pip(j-1)-pim(j)-pim(j-1))         &  
     &         + (gmid-gzmid)*(gip(j)*xsqoj(i+1)-giz(j)*xsqoj(i))        &  
     &         + xmidsq*(pr(i+1,j)-pr(i,j)))*gisa
!
   10 continue
      if(i.eq.2) return
      xc = xary(i-1)
      ajc = xc*dxdz
      gisc = deez/(deex*xc)
      do 20 j=3,nzp
      omeg1 = (omeg(i+1,j)-omeg(i,j))*abouna(j)
      omeg2 = (omeg(i,j)-omeg(i-1,j))*abounc(j)
      omegcent = omeg1+omeg2
      if(omegcent.gt.0) then
        omegdif = omeg2
                        else
        omegdif = omeg1
                          endif
!
!.......................................................................
!
!.....compute utsave and gtsave arrays for use in advan23
!
!.......................................................................
!
      gtsave(i,j) = gf1pa(j)*abouna(j)-gf2pa(j)*abounc(j)                &  
     &            + gf1ta(j)*abounb(j)-gf1ta(j-1)*abound(j)              &  
     &  -dxisq*x5*omegdif*dxdz*(1._R8/xa**2-1._R8/xc**2)*gzero
!
      utsave(i,j) = -((uf1pa(j)*abouna(j) - uf2pa(j  )*abounc(j))        &  
     &            +   (uf1ta(j)*abounb(j) - uf1ta(j-1)*abound(j)))
!
   20 continue
!.......................................................................
!
!...compute energy transport terms for isurf.eq.0
!
!.......................................................................
      if(isurf.ne.0.or.(ipres.eq.1.and.(idens.eq.1 .or. idens.eq.2))) go to 63
      call ener2d
   63 continue
      do 64 j=3,nzp
      if(isurf.eq.0 .and. idens.eq.0) go to 65
      fac7a(j) = 1._R8+eqrate
      rt(j) = 0._R8
      term9(j) = rsurf(i,j)*eqrate
      term3a(j) = 0._R8
      term4(j) = 0._R8
   65 continue
      if(isurf.eq.0 .and. ipres.eq.0) go to 64
!
!.....set terms to force pressure to equal value in qsurf array
      fac7(j) = 1._R8+eqrate
      qt2(j) = 0._R8
      qt4(j) = 0._R8
      qt(j) = 0._R8
      term8(j) = qsurf(i,j)*eqrate
      term5(j) = 0._R8
      term6(j) = 0._R8
   64 continue
      do 70 j=3,nzp
!.......................................................................
!
!.....advance r,q,w  (leap frog method)
!
!.......................................................................
      denom1= (fac7a(j) - .5_R8*(dt+dtold)*term4(j))
      fac5 = (1._R8+.5_R8*(dt+dtold)*term4(j))/denom1
      fac6 = (dt+dtold)/denom1
      term9(j) = term9(j)/denom1
!
      denom = (fac7(j) - .5_R8*(dt+dtold)*term6(j))
      fac9 = (1._R8+0.5_R8*(dt+dtold)*term6(j))/denom
      fac8 = (dt+dtold)/denom
      term8(j) = term8(j)/denom
!
      r(i,j) = fac5*ro(i,j) + fac6*(rt(j)+term3a(j)) + term9(j)
!
      q(i,j) = fac9*qo(i,j) + fac8*(qt(j)+term5(j)) + term8(j)
!
      ro(i,j) = riz(j)
      qo(i,j) = qiz(j)
   70 continue
!
      fac = amux*x5**2*acoef(62)
      terma = fac*aja*dxisq/xa**2
      termb = fac*aj5*dzisq/x5**2
      termc = fac*ajc*dxisq/xc**2
      termd = fac*aj5*dzisq/x5**2
      denom = 1._R8+.5_R8*(dt+dtold)*((terma+termb+termc+termd)/aj5      &  
     &                          + acoef(90))
      do 80 j=3,nzp
      fac1 = (2._R8-denom)*abounf(j)/denom
      fac2 = (dt+dtold)*abounf(j)/denom
!
      wt = (wf1pa(j)-wf2pa(j)+wf1ta(j)-wf1ta(j-1))
!
      w(i,j) = fac1*wo(i,j) + fac2*(wt+                                  &  
     &       wiz(j-1)/aj5*termd                                          &  
     &      +wim(j  )/aj4*termc                                          &  
     &      +wip(j  )/aj6*terma                                          &  
     &      +wiz(j+1)/aj5*termb )
!
      wo(i,j) = wiz(j)
   80 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
