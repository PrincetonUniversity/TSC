!#include "library_names.h"
      subroutine advsfa(i_newton)
!
!.....advances surface averaged adiabatic quantities in time 1 time step
!
!
!.......................................................................
!.....notation:
!
!.amat(m,n,j)*y(n,j+1) - bmat(m,n,j)*y(n,j) + cmat(m,n,j)*y(n,j-1)+dvec(m,j)=0
!
!.....where j is flux surface lables and  y(1,j) = adn(j)  ... differential number density
!                                         y(2,j) = adp(j)  . . total entropy density
!                                         y(3,j) = ade(j)  . . electron entropy density
!                                         y(4,j) = adi(j)  . . rotational transform
!                                         y(5,j) = ado(j)  . . toroidal angular momentum
!                                         y(6,j) = d2pi(j) . . 2nd deriv. of ion  press
!                                         y(7,j) = d2pe(j) . . 2nd deriv. of elec press
!.......................................................................
!
!
      USE CLINAM
      USE RUNAWAY
      USE SAPROP
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iadvsfae,npsp,j,m,n,ier,noeqn,noeqnp,l,jj,i_newton,npsihm
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rdpsi2,rdpsq,psi1,psi2,phirpgc,fac1,fac2,fvac,adnv
      REAL*8 adpv,adev,adiv,fac4,fac3,fac5,aiph,aimh,rsph,rsmh,psph,psmh
      REAL*8 esph,esmh,pmh,pph,emh,eph,rmh,rph,pval,pvalm,eps
      REAL*8 facj2,facc1,t1,t2,t3,t4,t5,facj,ppval,rval
      REAL*8 eval, epval
      REAL*8 rpval,fac,fluxn, factor
      REAL*8 sum, AREAL
      REAL*8 piped, peped, piedge, peedge
      REAL*8 piref(ppsi),peref(ppsi)
!============
      data iadvsfae /0/
!============
      INTERFACE
      subroutine f04aae(a1,n1,a2,n2,n3,n4,a3,n5,a4,n6)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n1,n2,n3,n4,n5,n6
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 a1(n1,*),a2(n2,*),a3(n5,*),a4(*)
      INTEGER ipiv(n1)
      INTEGER i, j
      REAL*8 d
      END subroutine f04aae

      END INTERFACE
!
!.....define old time variables
!
!
!......define energy source due to auxialliary heating
!
!     call auxheat
      if(i_newton .eq. 1) call auxheat
!......define edge density source
      if((idens.eq.0 .or. idens.eq.3) .and. i_newton .eq. 1) call edgesource
!
      npsp = npsit+2
      if(npsp.gt.ppsi) npsp = ppsi
      do 5 j=1,npsp
      if(i_newton .eq. 1) then
      adpo(j) = adp(j)
      adeo(j) = ade(j)
      adno(j) = adn(j)
      adio(j) = adi(j)
      adoo(j) = ado(j)
      endif
      do 5 m=1,pfour
      dvec(m,j) = 0._R8
      do 5 n=1,pfour
      amat(m,n,j) = 0._R8
      bmat(m,n,j) = 0._R8
      cmat(m,n,j) = 0._R8
    5 continue
!
      piref = 0
      peref = 0
      npsihm = int(pwidthc*AREAL(npsit))
!
!     define reference ion and electron profiles if pedestal is present
!.....temporary limit pedestal temperature to 5000ev
      if(npsihm .gt. 0) then
        piped = min((adp(npsihm)-ade(npsihm))/vpg(npsihm),tped*usdh*ane(npsihm))
        peped = min(ade(npsihm)/vpg(npsihm), tped*usdh*ane(npsihm))
        piedge =  (acoef(882)-1)*fracn0*r0*smallt
        peedge =   fracn0*r0*smallt
        do j=1,npsihm
          piref(j) = piped
          peref(j) = peped
        enddo
        do j=npsihm+1,npsit
          factor = AREAL((j-npsihm)**2)/AREAL((npsit - npsihm)**2)
          piref(j) = piped + factor*(piedge - piped)
          peref(j) = peped + factor*(peedge - peped)
        enddo 
        do j=npsit+1,npsp
          piref(j) = piedge
          peref(j) = peedge
        enddo
      endif
     
      rdpsi2 = .5_R8/dpsi
      rdpsq = 1._R8/dpsi**2
      psi1 = (.5_R8*dpsi)
      psi2 = (1.5_R8*dpsi)
!
!
!.....implicit parameter
      phirpgc = 1._R8-phirpg
!
!.....implicit solution for rs,es,ps
!
      fac1 = dtsf*phirpg*rdpsi2
      fac2 = dtsf*phirpg*rdpsq
!
      do 230 j=2,npsp
      fvac = 0._R8
!
!.....check if zone is in vacuum
      if(j.le.npsit) go to 101
      fvac = 1.E6_R8
      adnv = fracn0*r0*vp(j)
      adpv = acoef(882)*fracn0*r0*smallt*vpg(j)
      adev = fracn0*r0*smallt*vpg(j)
      adiv = 1._R8/qprof2(npsit)
  101 continue
      fac4 = (2._R8/3._R8)*vpg(j)/vp(j)
      fac3 = fac2*fac4
      fac5 = fac4*dtsf*acoef(3101)
!
      aiph = .5_R8*(adi(j)+adi(j+1))
      aimh = .5_R8*(adi(j)+adi(j-1))
      rsph = .5_R8*(adn(j)+adn(j+1))
      rsmh = .5_R8*(adn(j)+adn(j-1))
      psph = .5_R8*(adp(j)+adp(j+1))
      psmh = .5_R8*(adp(j)+adp(j-1))
      esph = .5_R8*(ade(j)+ade(j+1))
      esmh = .5_R8*(ade(j)+ade(j-1))
      pmh = .5_R8*(adp(j)/vpg(j)+adp(j-1)/vpg(j-1))
      pph = .5_R8*(adp(j)/vpg(j)+adp(j+1)/vpg(j+1))
      emh = .5_R8*(ade(j)/vpg(j)+ade(j-1)/vpg(j-1))
      eph = .5_R8*(ade(j)/vpg(j)+ade(j+1)/vpg(j+1))
      rmh = .5_R8*(adn(j)/vp (j)+adn(j-1)/vp(j-1))
      rph = .5_R8*(adn(j)/vp (j)+adn(j+1)/vp(j+1))
      pval = j*dpsi
      pvalm = (j-1._R8)*dpsi
 
!
!.......................................................................
!.......mass density
!.......................................................................
!
      amat(1,1,j) = fac1*( - cs0(j))                                     &  
     &            - fac2*cs1(j)*rsph /vp(j+1)
      cmat(1,1,j) = -fac1*( - cs0(j-1))                                  &  
     &            - fac2*cs1(j-1)*rsmh /vp(j-1)
      bmat(1,1,j) = -1._R8+fvac - fac1*(-cs0(j) + cs0(j-1))              &  
     &    - fac2*(rsph*cs1(j)+rsmh*cs1(j-1))/vp(j)
!
      amat(1,2,j) = -fac2*rsph*cs2(j)/vpg(j+1)
      cmat(1,2,j) = -fac2*rsmh*cs2(j-1)/vpg(j-1)
      bmat(1,2,j) = -fac2*(rsph*cs2(j) + rsmh*cs2(j-1))/vpg(j)
!
      amat(1,3,j) = -fac2*rsph*cs3(j)/vpg(j+1)
      cmat(1,3,j) = -fac2*rsmh*cs3(j-1)/vpg(j-1)
      bmat(1,3,j) = -fac2*(rsph*cs3(j) + rsmh*cs3(j-1))/vpg(j)
!
      amat(1,4,j) = -fac2*rsph*cs4(j)*gxmja(j+1)*xmja(j+1)
      cmat(1,4,j) = -fac2*rsmh*cs4(j-1)*gxmja(j-1)*xmja(j-1)
      bmat(1,4,j) = -fac2*(rsmh*cs4(j-1)+rsph*cs4(j))*gxmja(j)*xmja(j)
!
      dvec(1,j) = -adno(j) + fvac*adnv - dtsf*phirpgc*rdpsi*             &  
     &          (rsph*cs(j) - rsmh*cs(j-1))                              &  
     &     - dtsf*vp(j)*(srave(j) + sraveb(j) + simpe(j) + sraveedg(j)   &  
     &                 + sravejet(j) )
!
 
!.......................................................................
!.......energy density
!.......................................................................
!
!
!
      amat(2,2,j) = fac3*vp2(j)*(dsi2(j)+dse2(j))/vpg(j+1)
      cmat(2,2,j) = fac3*vp2(j-1)*(dsi2(j-1)+dse2(j-1))/vpg(j-1)
      bmat(2,2,j) = -1._R8+fvac + fac3*(vp2(j)*(dsi2(j)+dse2(j))         &  
     &           +vp2(j-1)*(dsi2(j-1)+dse2(j-1)) )/vpg(j)
!
      amat(2,1,j) = fac3*vp2(j)*(dsi1(j)+dse1(j))/vp(j+1)
      cmat(2,1,j) = fac3*vp2(j-1)*(dsi1(j-1)+dse1(j-1))/vp(j-1)
      bmat(2,1,j) = fac3/vp(j)*(vp2(j)*(dsi1(j)+dse1(j))                 &  
     &                         +vp2(j-1)*(dsi1(j-1)+dse1(j-1)) )
!
      amat(2,3,j) = fac3*vp2(j)*(dsi3(j)+dse3(j))/vpg(j+1)    
      cmat(2,3,j) = fac3*vp2(j-1)*(dsi3(j-1)+dse3(j-1))/vpg(j-1)
      bmat(2,3,j) = fac3/vpg(j)*(vp2(j)*(dsi3(j)+dse3(j))                &  
     &            + vp2(j-1)*(dsi3(j-1)+dse3(j-1)) )
!
      amat(2,4,j) = fac3*vp2(j)*(dse4(j)+dsi4(j))                        &  
     &                            *gxmja(j+1)*xmja(j+1)
      cmat(2,4,j) = fac3*vp2(j-1)*(dse4(j-1)+dsi4(j-1))                  &  
     &                            *gxmja(j-1)*xmja(j-1)
      bmat(2,4,j) = fac3*(vp2(j)*(dse4(j)+dsi4(j))                       &  
     &            +vp2(j-1)*(dse4(j-1)+dsi4(j-1)))*gxmja(j)*xmja(j)
!
!
      dvec(2,j) = -adpo(j) + fvac*adpv                                   &  
     &  + dtsf*rdpsi*                                                    &  
     &    fac4*(vp2(j)  *( phirpgc*(dsi(j)  +dse(j)  )+ phirpg*(dsi0(j)  +dse0(j  )))   &  
     &         -vp2(j-1)*( phirpgc*(dsi(j-1)+dse(j-1))+ phirpg*(dsi0(j-1)+dse0(j-1))))  &  
     &  -dtsf*fac4*( (.5_R8*(as(j-1)+as(j))+vlooph(j))/tpi*              &  
     &    (gxmja2(j) - gxmja2(j-1))*rdpsi                                &  
     & + vp(j)*(savei(j)+savee(j)+(savia(j) + savea(j))*ialpha           &  
     & + savefw(j) + savebm(j) + savibm(j) + savelh(j) + savilh(j)       &  
     & + savifw(j)                                                       &
     & - sradion(j)*usdp/usdt - savebre(j) - savecyc(j) - saveimp(j) ))
!
!....add terms due to hyper-conductivity
      amat(2,6,j) = fac5*vp2(j)*d2s(j)
      cmat(2,6,j) = fac5*vp2(j-1)*d2s(j-1)
      bmat(2,6,j) = fac5*(vp2(j)*d2s(j)+vp2(j-1)*d2s(j-1))
      amat(2,7,j) = fac5*vp2(j)*d2s(j)
      cmat(2,7,j) = fac5*vp2(j-1)*d2s(j-1)
      bmat(2,7,j) = fac5*(vp2(j)*d2s(j)+vp2(j-1)*d2s(j-1))
!
!
!
!.......................................................................
!......electron energy density
!.......................................................................
!
      facj2 = -.5_R8*rdpsi*dtsf*fac4*phirpg*(.5_R8*(as(j-1)+as(j))+      &  
     & vlooph(j))
      facc1 = fac1*.5_R8*(cs(j-1)+cs(j))*vpg(j)*(2._R8/3._R8)
!
      amat(3,3,j) = fac3*vp2(j)*dse3(j)/vpg(j+1) + facc1/vpg(j+1) 
      cmat(3,3,j) = fac3*vp2(j-1)*dse3(j-1)/vpg(j-1)-facc1/vpg(j-1)
      bmat(3,3,j) = -1._R8+fvac + fac3*(vp2(j)*dse3(j)                   &  
     &  +vp2(j-1)*dse3(j-1))/vpg(j)-(1._R8+avez(j))*equila(j)*dtsf
!
      amat(3,1,j) = fac3*vp2(j)*dse1(j)/vp(j+1)
      cmat(3,1,j) = fac3*vp2(j-1)*dse1(j-1)/vp(j-1)
      bmat(3,1,j) = fac3/vp(j)*                                          &  
     &   (vp2(j)*dse1(j) + vp2(j-1)*dse1(j-1))
!
      amat(3,2,j) = fac3*vp2(j)*dse2(j)/vpg(j+1)                         &  
     &  - facc1/vpg(j+1)
      cmat(3,2,j) = fac3*vp2(j-1)*dse2(j-1)/vpg(j-1)                     &  
     &  + facc1/vpg(j-1)
      bmat(3,2,j) = fac3/vpg(j)*                                         &  
     &  (vp2(j)*dse2(j) + vp2(j-1)*dse2(j-1)) + avez(j)*equila(j)*dtsf
!
      amat(3,4,j) = fac3*vp2(j)*dse4(j)*gxmja(j+1)*xmja(j+1)
      cmat(3,4,j) = fac3*vp2(j-1)*dse4(j-1)*gxmja(j-1)*xmja(j-1)
      bmat(3,4,j) = fac3*(vp2(j)*dse4(j)                                 &  
     &                   +vp2(j-1)*dse4(j-1))*gxmja(j)*xmja(j)
!
      dvec(3,j) = -adeo(j) + fvac*adev                                   &  
     &      + dtsf*rdpsi*                                                &  
     &        fac4*(vp2(j)*  (phirpgc*dse(j) + phirpg*dse0(j))           &
     &            - vp2(j-1)*(phirpgc*dse(j-1)+phirpg*dse0(j-1)) )       &  
     &      - dtsf*fac4*((.5_R8*(as(j-1)+as(j)) + vlooph(j))/tpi*        &  
     &           (gxmja2(j) - gxmja2(j-1))*rdpsi                         &  
     &   +phirpgc*(cs(j-1)+cs(j))*(pph-pmh-eph+emh)                      &  
     &  *vp(j)*rdpsi2 + vp(j)*(savee(j) + savea(j)*ialpha                &  
     & + savefw(j) + savebm(j) + savelh(j)                               &  
     &  -sradion(j)*usdp/usdt - savebre(j) - savecyc(j) - saveimp(j)))
!
!....add terms due to hyper-conductivity
      amat(3,7,j) = fac5*vp2(j)*d2s(j)
      cmat(3,7,j) = fac5*vp2(j-1)*d2s(j-1)
      bmat(3,7,j) = fac5*(vp2(j)*d2s(j)+vp2(j-1)*d2s(j-1))
!
!.......................................................................
!.......poloidal flux density . . .
!.......................................................................
!
      amat(4,1,j) = -fac2*as1(j)/vp(j+1)
      cmat(4,1,j) = -fac2*as1(j-1)/vp(j-1)
      bmat(4,1,j) = -fac2*(as1(j-1)+as1(j))/vp(j)
!
      amat(4,2,j) = -fac2*as2(j)/vpg(j+1)
      cmat(4,2,j) = -fac2*as2(j-1)/vpg(j-1)
      bmat(4,2,j) = -fac2*(as2(j-1)+as2(j))/vpg(j)
!
      amat(4,3,j) = -fac2*as3(j)/vpg(j+1)
      cmat(4,3,j) = -fac2*as3(j-1)/vpg(j-1)
      bmat(4,3,j) = -fac2*(as3(j-1)+as3(j))/vpg(j)
!
!
      t1 = rdpsi*etpara(j)*etafac(j)*tpi**2/(xmja2(j)*aiph)**2
      t2 = rdpsi*etpara(j-1)*etafac(j-1)*tpi**2/(xmja2(j-1)*aimh)**2
      t3 = gxmja(j+1)*xmja(j+1)
      t4 = gxmja(j-1)*xmja(j-1)
      t5 = gxmja(j)*xmja(j)
      facj = -phirpg*dtsf*rdpsi
!
      amat(4,4,j) = facj*t1*t3
      cmat(4,4,j) = facj*t2*t4
      bmat(4,4,j) = -1._R8+fvac + facj*t5*(t1+t2)
      dvec(4,j)=-adio(j)-phirpgc*dtsf*rdpsi*(as(j)-as(j-1))+fvac*adiv
      if(j.gt.2) go to 231
      amat(4,4,j) = facj*t1*t3
      cmat(4,4,j) = 0._R8
      bmat(4,4,j) = -1._R8+fvac + facj*t5*(t1+2._R8*t2)
  231 continue
!
!........................................................................
!.......rotation density
!.......................................................................
      bmat(5,5,j) = -1._R8
!
!
!.........................................................................
!.......second derivative of ion pressure - reference ion pressure
!.........................................................................
!
!
      amat(6,2,j) =   vp2(j)*d2s(j)/vpg(j+1)
      amat(6,3,j) = - vp2(j)*d2s(j)/vpg(j+1)
      bmat(6,6,j) = vp(j)*dpsi**2
      bmat(6,2,j) =  (vp2(j)*d2s(j)+vp2(j-1)*d2s(j-1))/vpg(j)
      bmat(6,3,j) = -(vp2(j)*d2s(j)+vp2(j-1)*d2s(j-1))/vpg(j)
      cmat(6,2,j) =   vp2(j-1)*d2s(j-1)/vpg(j-1)
      cmat(6,3,j) = - vp2(j-1)*d2s(j-1)/vpg(j-1)
      dvec(6,j)   = - ( vp2(j)  *d2s(j)  *piref(j+1) - (vp2(j)*d2s(j)+vp2(j-1)*d2s(j-1))*piref(j)   &
                      + vp2(j-1)*d2s(j-1)*piref(j-1))
!
!.........................................................................
!.......second derivative of electron pressure - reference electron pressure
!.........................................................................
!
!
      amat(7,3,j) =  vp2(j)*d2s(j)/vpg(j+1)
      bmat(7,7,j) =  vp(j)*dpsi**2
      bmat(7,3,j) = (vp2(j)*d2s(j)+vp2(j-1)*d2s(j-1))/vpg(j)
      cmat(7,3,j) =  vp2(j-1)*d2s(j-1)/vpg(j-1)
      dvec(7,j)   = - ( vp2(j)  *d2s(j)  *peref(j+1) - (vp2(j)*d2s(j)+vp2(j-1)*d2s(j-1))*peref(j)   &
                      + vp2(j-1)*d2s(j-1)*peref(j-1))
  230 continue
!
      ier = 1._R8
!
      noeqn = 7
!.......................................................................
!.....compute emat and fvec matrices
!.......................................................................
      noeqnp = noeqn + 1
      do 200 j=2,npsp
      do 190 m=1,noeqn
      do 190 n=1,noeqn
!
!.....store a in tmp4 temporarily
!
      tmp4(m,n) = amat(m,n,j)
      sum = 0._R8
      if(j.eq.2) go to 175
      do 170 l=1,noeqn
  170 sum = sum + cmat(m,l,j)*emat(l,n,j-1)
  175 continue
  190 tmp1(m,n) = bmat(m,n,j) - sum
!
      do 160 m=1,noeqn
      sum = 0._R8
      if(j.eq.2) go to 155
      do 150 n=1,noeqn
  150 sum = sum + cmat(m,n,j)*fvec(n,j-1)
!
!.....store d+c*f temporarily in tmp4
!
  155 continue
      tmp4(m,noeqnp) = dvec(m,j) + sum
  160 continue
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!.....call leqt1f(tmp1,noeqnp,noeqn,pfour,tmp4,0,tmp3,ier)
      ier    = 1
      call f04aae(tmp1,pfour,tmp4,pfour,noeqn,noeqnp,                    &  
     &            tmp5,pfour,tmp3,ier)
      if(ier.ne.0) then
       write(nout,5555) kcycle,j
       write(nterm,5555) kcycle,j
 5555  format("error in advsfa  ,kcycle,j=",i8,i6)
       do n=1,pfour
       do m=1,pfour
        write(nout,5556) m,n,amat(m,n,j),bmat(m,n,j),cmat(m,n,j)
       enddo
       enddo
       write(nout,5557) j,cs(j),cs(j-1),dsi(j),dsi(j-1),                 &  
     &                    dse(j),dse(j-1)
 5556  format("m,n= (",i1,",",i1, ")   a,b,c(m,n)",1p3e12.4)
 5557  format(i3," cs(j)...dse(j-1)",1p6e12.4)
       write(nout,5558) dsi(j),dsi(j-1),dse(j),dse(j-1),vp(j),vp(j-1),   &  
     &                  vp2(j),vp2(j-1)
 5558  format("dsi(j)...vp2(j-1)",1p8e10.2)
       if(j.le.npsit) iadvsfae = iadvsfae + 1
      endif
      if(iadvsfae.gt.100) ineg=15
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
      do 165 m=1,noeqn
      do 163 n=1,noeqn
      emat(m,n,j) = tmp5(m,n)
  163 continue
      fvec(m,j) = tmp5(m,noeqnp)
  165 continue
  200 continue
!.......................................................................
!
!.....define boundary values
!.......................................................................
      adnv = fracn0*r0*vp(npsp)
      adpv = acoef(882)*fracn0*r0*smallt*vpg(npsp)
      adev = fracn0*r0*smallt*vpg(npsp)
      adiv = 1._R8/qprof2(npsit)
      pvec(1,npsp) = adnv
      pvec(2,npsp) = adpv
      pvec(3,npsp) = adev
      pvec(4,npsp) = adiv
      pvec(5,npsp) = 0.
      pvec(6,npsp) = 0.
      pvec(7,npsp) = 0.
!
!.......................................................................
!.....now       use emat and fvec to compute pvec
!.......................................................................
      do 400 jj=2,npsp-1
      j = npsp + 1 - jj
      do 390 m=1,noeqn
      sum = 0._R8
      do 380 n=1,noeqn
  380 sum = sum + emat(m,n,j)*pvec(n,j+1)
  390 pvec(m,j) = sum + fvec(m,j)
      adn(j) = pvec(1,j)
      adp(j) = pvec(2,j)
      ade(j) = pvec(3,j)
      adi(j) = pvec(4,j)
      ado(j) = pvec(5,j)
      d2pi(j)= pvec(6,j)
      d2pe(j)= pvec(7,j)
  400 continue
!
!.....adjust variables to ensure positivity
!
      do 600 j=2,npsit
!
!.....override q calculation with that from surface average
      adi(j) = 2._R8/(qprof2(j)+qprof2(j-1))
!
!.....prescribe pressure for ipres .ne. 0
!
      if(ipres.le.0) go to 590
      call peval(xsv(j),1, pval,ppval,imag,jmag)
      adp(j) = pval*vpg(j)
      call eeval(xsv(j),1, eval,epval,imag,jmag)
      ade(j) = eval*vpg(j)
  590 continue
!
!.....prescribe density for idens .eq. 1 or 2
!
      if(idens.le.0 .or. idens.ge.3) go to 599
      call reval(xsv(j),idens,0,rval,rpval,imag,jmag)
      adn(j) = rval*vp(j)
      if(ipres.gt.0) go to 599
!
!.....reduce zone energy if density is decreasing
      if(adno(j).le.0 .or. kcycle.le.0) go to 599
      fac = adn(j)/adno(j)
      if(fac.gt.1) goto 599
      adp(j) = adp(j)*fac
      ade(j) = ade(j)*fac
  599 continue
      adnv = fracn0*r0*vp(j)
      if(adn(j).le.adnv) adn(j)=adnv
!
      adpv = acoef(882)*smallt*adn(j)/vp(j)*vpg(j)
      adev = smallt*adn(j)/vp(j)*vpg(j)
      if(adp(j).le.adpv) adp(j)=adpv
      if(ade(j).le.adev) ade(j)=adev
      if(ade(j).ge.adp(j)) ade(j) = adp(j)
  600 continue
!
!
!...SPECIAL DIAGNOSTIC PRINTOUT
!     write(nterm,9600) (j,adp(j),ade(j),adn(j),j=1,5)
!9600 format(" j,adp,ade,adn",i3,1p3e12.4)
!
      adn(1) = adn(2)
      adp(1) = adp(2)
      ade(1) = ade(2)
      adi(1) = adi(2)
!
!     calculate effective thermal conductivity due to hyper term
      chiesec = 1.e-8
      chiisec = 1.e-8
      if(acoef(3101) .gt. 0) then
        do j=2,npsit
          chiesec(j) = acoef(3101)*(d2pe(j+1)-d2pe(j))                        &
                 /(phirpg*rdpsq*(ade(j+1)/vpg(j+1) - ade(j)/vpg(j))) / udst
          chiisec(j) = acoef(3101)*(d2pi(j+1)-d2pi(j))                        & 
                 /(phirpg*rdpsq*((adp(j+1)-ade(j+1))/vpg(j+1)                 &
                                -(adp(j)  -ade(j)  )/vpg(j+1)))      /udst
        enddo
        chiesec(1) = chiesec(2)
        chiisec(1) = chiisec(2)
      endif
!
!.....flux of particles crossing last boundary
      npsp = npsit+1
      fluxn = .5_R8*(adno(npsp)+adno(npsit))/vp2(npsit)                  &  
     &      *(phirpgc*cs(npsit) + phirpg*(cs0(npsit)                     &  
     &  + cs4(npsit)*(adi(npsp)*gxmja(npsp)*xmja(npsp)                   &  
     &             - adi(npsit)*gxmja(npsit)*xmja(npsit))*rdpsi          &  
     &  + cs2(npsit)*(adp(npsp)/vpg(npsp)-adp(npsit)/vpg(npsit))*rdpsi   &  
     &  + cs1(npsit)*(adn(npsp)/vp (npsp)-adn(npsit)/vp (npsit))*rdpsi   &  
     &  + cs3(npsit)*(ade(npsp)/vpg(npsp)-ade(npsit)/vpg(npsit))*rdpsi))     
!
!.....heat flux crossing last boundary in watts
      j = npsit
      dsi(j-1) =                                                         &  
     &     dsi0(j-1)                                                     &  
     &   + dsi4(j-1)*(adi(j)*gxmja(j)*xmja(j)                            &  
     &              - adi(j-1)*gxmja(j-1)*xmja(j-1))*rdpsi               &  
     &   + dsi2(j-1)*(adp(j)/vpg(j)-adp(j-1)/vpg(j-1))*rdpsi             &  
     &   + dsi1(j-1)*(adn(j)/vp(j)-adn(j-1)/vp(j-1))*rdpsi               &  
     &   + dsi3(j-1)*(ade(j)/vpg(j)-ade(j-1)/vpg(j-1))*rdpsi
!
      dse(j-1) =                                                         &  
     &     dse0(j-1)                                                     &  
     &   + dse4(j-1)*(adi(j)*gxmja(j)*xmja(j)                            &  
     &              - adi(j-1)*gxmja(j-1)*xmja(j-1))*rdpsi               &  
     &   + dse2(j-1)*(adp(j)/vpg(j)-adp(j-1)/vpg(j-1))*rdpsi             &  
     &   + dse1(j-1)*(adn(j)/vp(j)-adn(j-1)/vp(j-1))*rdpsi               &  
     &   + dse3(j-1)*(ade(j)/vpg(j)-ade(j-1)/vpg(j-1))*rdpsi
      hfluxp = vp2(npsit-1)*(dsi(npsit-1)+dse(npsit-1))*udsp/udst
!
      hfluxp = 0.5_R8*(hfluxp + hfluxpo)
      hfluxpo = hfluxp
!>>>>>>>debug
!     write(73,1072)
!1072 format("   j      dsi         dsi0        dse         dse0")
!1073 format(i5,1p4e12.4)
!     do j=2,npsit
!     dsi(j-1) =                                                         &
!    &     dsi0(j-1)                                                     &
!    &   + dsi4(j-1)*(adi(j)*gxmja(j)*xmja(j)                            &
!    &              - adi(j-1)*gxmja(j-1)*xmja(j-1))*rdpsi               &
!    &   + dsi2(j-1)*(adp(j)/vpg(j)-adp(j-1)/vpg(j-1))*rdpsi             &
!    &   + dsi1(j-1)*(adn(j)/vp(j)-adn(j-1)/vp(j-1))*rdpsi               &
!    &   + dsi3(j-1)*(ade(j)/vpg(j)-ade(j-1)/vpg(j-1))*rdpsi

!     dse(j-1) =                                                         &
!    &     dse0(j-1)                                                     &
!    &   + dse4(j-1)*(adi(j)*gxmja(j)*xmja(j)                            &
!    &              - adi(j-1)*gxmja(j-1)*xmja(j-1))*rdpsi               &
!    &   + dse2(j-1)*(adp(j)/vpg(j)-adp(j-1)/vpg(j-1))*rdpsi             &
!    &   + dse1(j-1)*(adn(j)/vp(j)-adn(j-1)/vp(j-1))*rdpsi               &
!    &   + dse3(j-1)*(ade(j)/vpg(j)-ade(j-1)/vpg(j-1))*rdpsi

!     write(73,1073) j, dsi(j),dsi0(j),dse(j),dse0(j)
!     enddo
!     stop
!>>>>>>debug




!
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok
