      subroutine ener2d
!
!.....calculates energy terms when there is no surface averaging
!.....and ipres=1
!
      USE CLINAM
      USE SCR1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dxdzi,dxodz,rfac,x4,x5,x6,xa,xc,aj4,aj5,aj6,aja,ajc
      REAL*8 aj4i,aj5i,aj6i,ajai,ajci,tw1,tw2,tg1,tg2,tr1,tq1,taw1
      REAL*8 taw2,tag1,tag2,tar1,taq1,xsqa,xsqb,xsqc,xsqd,gjsd,gisa
      REAL*8 gisc,facq,gsi1,gsi2,gsj1,gsj2,ejpx1,ejpx2,ejpz1,ejpz2
      REAL*8 trm1,trm2,trm3,trm4,trm5,trm6,ce,pa,pb,pc,pd,ra,rb,rc
      REAL*8 rd,bsqa,bsqb,bsqc,bsqd,ta,tb,tc,td,gisqa,gisqc,gjsqb
      REAL*8 gjsqd,c1a,c1c,c4b,c4d,dparm,coa,cob,coc,cod,pjas,pias
      REAL*8 pjcs,pics,pjbs,pibs,pjds,pids
!============
      i = icol
!
!
      dxdzi = 1._R8/dxdz
      dxodz = deex/deez
      rfac  = (deez/deex)**2
!
!...coordinate functions
!
      x4 = xary(i-1) - 0.5_R8*deex
      x5 = xarh(i)
      x6 = xarh(i+1)
      xa    =       xary(i)
      xc = xary(i-1)
      x4 = abs(x4)
      aj4   = x4*dxdz
      aj5   = x5*dxdz
      aj6   = x6*dxdz
      aja   = xa*dxdz
      ajc = xc*dxdz
      aj4i  = 1._R8/aj4
      aj5i  = 1._R8/aj5
      aj6i  = 1._R8/aj6
      ajai  = 1._R8/aja
      ajci  = 1._R8/ajc
      tw1   = aj5i*x5**2
      tw2   = aj5i
      tg1   = aj5i/x5**2
      tg2   = aj5i
      tr1   = aj5i
      tq1   = aj5i
      taw1   = ajai*xa**2
      taw2   = ajai
      tag1   = ajai/xa**2
      tag2   = ajai
      tar1   = ajai
      taq1   = ajai
      xsqa = xa**2
      xsqb = x5**2
      xsqc = xc**2
      xsqd = x5**2
      gjsd = deex/(deez*x5)
      gisa = deez/(deex*xa)
      gisc = deez/(deex*xc)
      do 30 j=3,nzp
!
!...ohmic heating contribution
!
      facq = (2._R8/5._R8)*pr(i,j)**(-2._R8/5._R8)
!
      gsi1=(gip(j)*x6-giz(j)*x5)/deex*abouna(j)
      gsi2=(giz(j)*x5-gim(j)*x4)/deex*abounc(j)
      gsj1=(giz(j+1)   -      giz(j)        )*x5/deez*abounb(j)
      gsj2=(giz(j)   -      giz(j-1)        )*x5/deez*abound(j)
       ejpx1=2._R8/(1._R8/etay(i  ,j-1)+1._R8/etay(i  ,j  ))
       ejpx2=2._R8/(1._R8/etay(i-1,j-1)+1._R8/etay(i-1,j  ))
       ejpz1=2._R8/(1._R8/etay(i  ,j  )+1._R8/etay(i-1,j  ))
       ejpz2=2._R8/(1._R8/etay(i  ,j-1)+1._R8/etay(i-1,j-1))
      qt2(j) = facq*aj5*.25_R8*                                          &  
     &  (etay(i  ,j)*ajphi(i  ,j)**2+etay(i  ,j-1)*ajphi(i  ,j-1)**2     &  
     &  +etay(i-1,j)*ajphi(i-1,j)**2+etay(i-1,j-1)*ajphi(i-1,j-1)**2)    &  
     &  + facq*0.5_R8*(x5**2/aj5)*(gsj1**2*ejpx1                         &  
     &  +                       gsj2**2*ejpx2                            &  
     &  +                       gsi1**2*ejpz1                            &  
     &  +                       gsi2**2*ejpz2)
!
!...viscous heating contribution  (9/13/92)
      trm1 = amux*(uiz(j)*aj5i)**2
      trm2 = acoef(10)*amux*(.25_R8/x5**2)*(                             &  
     &  (biz(j)*xa**2/aja)**2                                            &  
     & +(biz(j-1)*xa**2/aja)**2                                          &  
     & +(bim(j)*xc**2/ajc)**2                                            &  
     & +(bim(j-1)*xc**2/ajc)**2 )
      trm3 = acoef(62)*amux/x5**2*                                       &  
     & (((wip(j)/aj6  - wim(j)/aj4) / (2._R8*deex))**2                   &  
     & +((wiz(j+1) - wiz(j-1))/(2*deez*aj5))**2)
      trm4 = acoef(92)/x5**2*                                            &  
     &   .5_R8*(((abip(j) - aiz(j))/deex)**2                             &  
     &     + ((abip(j-1)-aiz(j-1))/deex)**2                              &  
     &     + ((abip(j) - abip(j-1))/deez)**2                             &  
     &     + ((aiz(j)  - aiz(j-1))/deez)**2)
      trm5 = acoef(91)*((.5_R8*(omeg(i+1,j) - omeg(i-1,j))/deex)**2      &  
     &                 + (.5_R8*(omeg(i,j+1) - omeg(i,j-1))/deez)**2)
      trm6 = acoef(90)/x5**2*(w(i,j)/aj5)**2
!
!...note: cross terms in drag are missing
      qt2(j) = qt2(j) + facq*aj5*(trm1+trm2+trm3+trm4+trm5+trm6)
!
   30 continue
      do 40 j=3,nzp
      fac7(j) = 1._R8
      fac7a(j) = 1._R8
      term8(j) = 0._R8
   40 continue
!
!
      do 50 j=3,nzp
!...compute heat conduction contribution
!
      ce = acoef(890)*(30._R8)*(0.4_R8)*q(i,j)/(pr(i,j)*aj5)
!
      pa = .5_R8*(pr(i,j)+pr(i+1,j))
      pb = .5_R8*(pr(i,j)+pr(i,j+1))
      pc = .5_R8*(pr(i,j)+pr(i-1,j))
      pd = .5_R8*(pr(i,j)+pr(i,j-1))
!
      ra = .5_R8*(roj(i,j)+roj(i+1,j))
      rb = .5_R8*(roj(i,j)+roj(i,j+1))
      rc = .5_R8*(roj(i,j)+roj(i-1,j))
      rd = .5_R8*(roj(i,j)+roj(i,j-1))
!
!
      bsqa =    .5_R8*(bmagy(i  ,j  ) + bmagy(i+1,j  ))                  &  
     &       +( .5_R8*(   gs(i  ,j  ) +    gs(i  ,j-1)) )**2/xsqa
      bsqb =    .5_R8*(bmagy(i  ,j  ) + bmagy(i  ,j+1))                  &  
     &       +( .5_R8*(   gs(i  ,j  ) +    gs(i-1,j  )) )**2/xsqb
      bsqc =    .5_R8*(bmagy(i  ,j  ) + bmagy(i-1,j  ))                  &  
     &       +( .5_R8*(   gs(i-1,j  ) +    gs(i-1,j-1)) )**2/xsqc
      bsqd =    .5_R8*(bmagy(i  ,j  ) + bmagy(i  ,j-1))                  &  
     &       +( .5_R8*(   gs(i  ,j-1) +    gs(i-1,j-1)) )**2/xsqd
!
      ta = pa/(bsqa+acoef(891)*pa)
      tb = pb/(bsqb+acoef(891)*pb)
      tc = pc/(bsqc+acoef(891)*pc)
      td = pd/(bsqd+acoef(891)*pd)
!
!
      gisqa = .5_R8*(etay(i,j)+etay(i,j-1))*dxisq*aja
      gisqc = .5_R8*(etay(i-1,j)+etay(i-1,j-1))*dxisq*ajc
      gjsqb = .5_R8*(etay(i,j)+etay(i-1,j))*dzisq*aj5
      gjsqd = .5_R8*(etay(i-1,j-1)+etay(i,j-1))*dzisq*aj5
!
      c1a = abouna(j)*gisqa*ta*ra
      c1c = abounc(j)*gisqc*tc*rc
      c4b = abounb(j)*gjsqb*tb*rb
      c4d = abound(j)*gjsqd*td*rd
!
      qt4(j) = ce*                                                       &  
     &     (pro(i+1,j)/roj(i+1,j)*(c1a)                                  &  
     &     +pro(i-1,j)/roj(i-1,j)*(c1c)                                  &  
     &     +pro(i,j+1)/roj(i,j+1)*(c4b)                                  &  
     &     +pro(i,j-1)/roj(i,j-1)*(c4d)                                  &  
     &     +pro(i,j)/roj(i,j)*(-c1a-c1c-c4b-c4d))
!
      qt(j) = qf1pa(j)*abouna(j)-qf2pa(j)*abounc(j)                      &  
     &       + qf1ta(j)*abounb(j)-qf1ta(j-1)*abound(j)+qt2(j)+qt4(j)
!
      rt(j) = rf1pa(j)*abouna(j)-rf2pa(j)*abounc(j)                      &  
     &   + rf1ta(j)*abounb(j)-rf1ta(j-1)*abound(j)
   50 continue
!
!.....compute parallel diffusion
!
      do 60 j=3,nzp
      dpara(j) = dpar*.5_R8*(etay(i,j)+etay(i,j-1))*abouna(j)
      dparb(j) = dpar*.5_R8*(etay(i,j)+etay(i-1,j))*abounb(j)
      dparc(j) = dpar*.5_R8*(etay(i-1,j)+etay(i-1,j-1))*abounc(j)
      dpard(j) = dpar*.5_R8*(etay(i-1,j-1)+etay(i,j-1))*abound(j)
   60 continue
      dparm = 2._R8*etav/ndiv
      do 61 j=3,nzp
      if(dpara(j).gt.dparm) dpara(j) = dparm
      if(dparb(j).gt.dparm) dparb(j) = dparm
      if(dparc(j).gt.dparm) dparc(j) = dparm
      if(dpard(j).gt.dparm) dpard(j) = dparm
   61 continue
!
!
      do 62 j=3,nzp
      coa=dpara(j)/(aja*.5_R8*(bmagy(i,j)+bmagy(i+1,j)))*(piz(j)-piz(j-  &  
     & 1))
      cob=dparb(j)/(aj5*.5_R8*(bmagy(i,j)+bmagy(i,j+1)))*(piz(j)-pim(j))     
       
      coc=dparc(j)/(ajc*.5_R8*(bmagy(i,j)+bmagy(i-1,j)))*(pim(j)-pim(j-  &  
     & 1))
      cod=dpard(j)/(aj5*.5_R8*(bmagy(i,j)+bmagy(i,j-1)))                 &  
     &    *(piz(j-1)-pim(j-1))
!
      pjas = coa*(piz(j)-piz(j-1))
      pias = coa*(pip(j)+pip(j-1)-pim(j)-pim(j-1))*.25_R8
      pjcs = coc*(pim(j)-pim(j-1))
      pics = coc*(piz(j)+piz(j-1)-psi(i-2,j)-psi(i-2,j-1))*.25_R8
      pjbs = cob*(pim(j+1)+piz(j+1)-pim(j-1)-piz(j-1))*.25_R8
      pibs = cob*(piz(j)-pim(j))
      pjds = cod*(pim(j)+piz(j)-psi(i-1,j-2)-psi(i,j-2))*.25_R8
      pids = cod*(piz(j-1)-pim(j-1))
!
      term3a(j)= riz(j-1)*(pids+.25_R8*pias-.25_R8*pics)/aj5             &  
     &      + rim(j)*(pjcs+.25_R8*pjbs-.25_R8*pjds)/aj4                  &  
     &      + rip(j)*(pjas-.25_R8*pjbs+.25_R8*pjds)/aj6                  &  
     &      + riz(j+1)*(pibs-.25_R8*pias+.25_R8*pics)/aj5                &  
     &      + rim(j-1)*(-.25_R8*pics-.25_R8*pjds)/aj4                    &  
     &      + rip(j-1)*( .25_R8*pias+.25_R8*pjds)/aj6                    &  
     &      + rim(j+1)*( .25_R8*pics+.25_R8*pjbs)/aj4                    &  
     &      + rip(j+1)*(-.25_R8*pias-.25_R8*pjbs)/aj6
      term4(j)= -(pjas+pjcs+pibs+pids)/aj5
      term5(j)= qiz(j-1)*(pids+.25_R8*pias-.25_R8*pics)/aj5              &  
     &      + qim(j)*(pjcs+.25_R8*pjbs-.25_R8*pjds)/aj4                  &  
     &      + qip(j)*(pjas-.25_R8*pjbs+.25_R8*pjds)/aj6                  &  
     &      + qiz(j+1)*(pibs-.25_R8*pias+.25_R8*pics)/aj5                &  
     &      + qim(j-1)*(-.25_R8*pics-.25_R8*pjds)/aj4                    &  
     &      + qip(j-1)*( .25_R8*pias+.25_R8*pjds)/aj6                    &  
     &      + qim(j+1)*( .25_R8*pics+.25_R8*pjbs)/aj4                    &  
     &      + qip(j+1)*(-.25_R8*pias-.25_R8*pjbs)/aj6
      term6(j) = -(pjas+pjcs+pibs+pids)/aj5
   62 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
