      subroutine cdbm_chi(model,dpdr,dptdr,wexb,chi0,fs1,fee,chicdbm)
      
      USE CLINAM
      USE SAPROP

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      
      REAL*8, PARAMETER :: aee = 1.602176487e-19

      INTEGER j,model,jmin_cdbm,jlim_cdbm,npsihm,j1,j2

      REAL*8 delta2,rhodens,va,eps0,wpe2,xchi,dchi
      REAL*8 fej,fsj,calf,ckcdbm,sa_min
      REAL*8 delta,ellip,x1,x2,z1,z2
      REAL*8 curv,AREAL,fs_jmin,csval,pval1,ppval1,eval1,peval1
      REAL*8 aexb,gexb,rmaj,shearl,alphaj
      REAL*8, DIMENSION(npsit) :: dqdr,fs,fe,fk,wexb,dEdr,dpdr,dptdr
      REAL*8, DIMENSION(npsit) :: fs1,fee,chi0s
      REAL*8, DIMENSION(npsit) :: alpha,shear,alphat
      REAL*8, DIMENSION(npsit) :: chicdbm,chi0,chicdbms
!      
      eps0 = udsp/(2.998e8_R8)**2
      npsihm = int(pwidthc * AREAL(npsit))
!     dqdr(1) = 0.0_R8
      do j=1,npsit
        if(j .le. 2) then
          dqdr(j) = (-qprof2(j+2)+4.*qprof2(j+1)-3.*qprof2(j))          &
     &              / (rminora(j+2)-rminora(j))
        elseif (j .gt. npsit-2) then
          dqdr(j) = (3.*qprof2(j)-4.*qprof2(j-1)+qprof2(j-2))           &
     &              / (rminora(j)-rminora(j-2))
        else
          dqdr(j)=(-qprof2(j+2)+8.*qprof2(j+1)-8.*qprof2(j-1)+qprof2(j-2))&
     &                /(3.*(rminora(j+2)-rminora(j-2)))
        end if
      end do
!
      do j=1,npsit
         fk(j) = 1.0_R8
         fe(j) = 1.0_R8
         fs(j) = 1.0_R8
         chi0(j)=1.0_R8
      end do
!


      j1 = int(acoef(3204) * AREAL(npsit))
      j2 = int(acoef(3203) * AREAL(npsit))
      ckcdbm = acoef(3201)
      if (ckcdbm .eq. 0_R8) ckcdbm=12._R8
      if(zeff .le. 1.) zeff = 1.5
      calf = 1.0
      gexb=acoef(3202)
      do j=1,npsit
        rmaj = xmag+rminora(j)
!       rhodens = 2.0*1.67e-27_R8*(ane(j)/zeff)
!       va = sqrt((gzero/rmaj)**2*udsp/rhodens)
!       if (j .EQ. 1) then
!       write(*,*)'------------'
!       write(*,*) 'rho1, SI= ',rhodens
!       write(*,*) 'ni1=',ane(j)/zeff
!       write(*,*) 'vA1= ',va      
!       end if
!        va = sqrt((gzero/rmaj)**2*udsp/rhodens)
        call peval(xsv2(j),2,pval1,ppval1,imag,jmag)
        call eeval(xsv2(j),2,eval1,peval1,imag,jmag)
        rhodens = 2.5*1.67e-27_R8*(pval1-eval1)/(ti(j)*usdh)
        va = sqrt((gzero/rmaj)**2*udsp/rhodens)
!       if (j .EQ. 1) then
!       write(*,*) 'pvali = ',(pval1-eval1)*udsp
!       write(*,*) 'rho2,SI = ',rhodens
!       write(*,*) 'ni2=',udsd*udsp*(pval1-eval1)/ti(j),(pval1-eval1)/(ti(j)*usdh)
!       write(*,*) 'vA2= ',va      
!       end if
!        write(*,*)'----------------'
!        write(*,*)'rhodens',rhodens,'va',va
        wpe2 = ane(j)*(1.6e-19_R8)**2/(9.1e-31_R8*eps0)
        delta2 = (2.998e8_R8)**2/wpe2
        alpha(j) = -2.0*usdp*qprof2(j)**2*rmaj**3*dpdr(j)/gzero**2
        alphat(j) = -2.0*usdp*qprof2(j)**2*rmaj**3*dptdr(j)/gzero**2
        curv = -(rminora(j)/rmaj)*(1.0-1.0/qprof2(j)**2)
        shear(j) = dqdr(j)*rminora(j)/qprof2(j)
        shearl = sqrt(shear(j)**2 + 0.01)
        chi0(j) = sqrt(abs(alpha(j)))**3*delta2*va/(qprof2(j)*rmaj)
!       write(*,*) 'j',j,'chi0',chi0(j)
!       write(*,*)'--rhodens',rhodens,'va',va
!       write(*,*) '--chi0',chi0(j)
        if (j .ge. j2) then
           alphaj = -2.0*usdp*qprof2(j)**2*rmaj**3*dpdr(j2)/gzero**2
           chi0(j)=sqrt(abs(alphaj))**3*delta2*va/(qprof2(j)*rmaj)
        end if
        if (model .ne. 0) then
!          fsj = trcofs(shear(j),calf*alphat(j),curv,rmaj,rminora(j))
           fsj = trcofs(shearl,calf*alphat(j),curv,rmaj,rminora(j))
           if(fsj .le. 0.) fsj = 1._R8
!          if (shear(j)-alpha(j) .gt. 0.5) fsj = 1.0_R8
           if (j .gt. j1) fsj = fs(j1)
           fs(j) = fsj
        end if
        aexb = (qprof2(j)*rmaj/(shearl*va))**2
!       if (acoef(3202) .eq. 0.0_R8) gexb = 6.*shearl/alphat(j)
        if (acoef(3202) .eq. 0.0_R8) gexb = 6.*shearl**3/alphat(j)
        if (model .eq. 2 .or. model .eq. 3)                             &
     &                fe(j) = 1.0_R8/(1.0_R8+gexb*aexb*wexb(j)**2) 
        if (model .eq. 4 .or. model .eq. 5)                             &
     &                fe(j) = gexb*fexb(wexb(j),shear(j),alpha(j))
        if (model .eq. 1 .or. model .eq. 3 .or. model .eq. 5) then
           call shapex(ellip,delta,x1,x2,z1,z2,xsv2(j))
           fk(j) = (2.*sqrt(ellip)/(1.+ellip**2))**1.5
        end if         
        if (j .gt. j1) fe(j) = 1.0_R8 
      end do
!
!      jmin_cdbm = 1
!      sa_min = 0._R8
!      do j=1,npsihm-20
!         if (shear(j)-alpha(j) .le. sa_min .and. j .gt. 10) then
!           sa_min = shear(j)-alpha(j)
!           jmin_cdbm = j
!         end if
!      end do
      chi0s(1)=chi0(1)
      fs1(1)=fs(1)
      fee(1)=fe(1)
      chi0s(2)=(chi0(1)+chi0(2)+chi0(3))/3._R8
      fs1(2)=(fs(1)+fs(2)+fs(3))/3.
      fee(2)=(fe(1)+fe(2)+fe(3))/3.
      chi0s(npsihm-1)=(chi0(npsihm-2)+chi0(npsihm-1)+chi0(npsihm))/3._R8
      fs1(npsihm-1)=(fs(npsihm-2)+fs(npsihm-1)+fs(npsihm))/3.
      fee(npsihm-1)=(fe(npsihm-2)+fe(npsihm-1)+fe(npsihm))/3.
      chi0s(npsihm)=chi0(npsihm)
      fs1(npsihm)=fs(npsihm)
      fee(npsihm)=fe(npsihm)
      do j=3,npsihm-2
         chi0s(j)=(chi0(j-2)+chi0(j-1)+chi0(j)+chi0(j+1)+chi0(j+2))/5.
         fs1(j)=(fs(j-2)+fs(j-1)+fs(j)+fs(j+1)+fs(j+2))/5.
         fee(j)=(fe(j-2)+fe(j-1)+fe(j)+fe(j+1)+fe(j+2))/5.
      end do
      do j=1,npsihm
         chicdbms(j) = ckcdbm*fs1(j)*fk(j)*fee(j)*chi0s(j)
      end do      
      chicdbm(1)=chicdbms(1)
      chicdbm(2)=(chicdbms(1)+chicdbms(2)+chicdbms(3))/3._R8      
      chicdbm(npsihm-1)=(chicdbms(npsihm-2)+chicdbms(npsihm-1)          &
     &                 + chicdbms(npsihm))/3._R8      
      chicdbm(npsihm)=chicdbms(npsihm)      
      do j=3,npsihm-2
         chicdbm(j)=(chicdbms(j-2)+chicdbms(j-1)+chicdbms(j)            &
     &              + chicdbms(j+1)+chicdbms(j+2))/5.0_R8
         if (j .gt. j1) chicdbm(j) = chicdbm(j1)
      end do

      do j=npsihm,npsit
         fs1(j)=fs1(npsihm-2)
         fee(j)=fee(npsihm-2)  
         chicdbm(j)=chicdbm(npsihm-2)
      end do
!      
!      return
      
 CONTAINS

      REAL FUNCTION trcofs(shear,alpha,curv,rmakk,rmikk)

     IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 alpha,curv,shear,sa,fs1,fs2,rmakk,rmikk


      if (alpha .ge. 0.0) then
         sa = shear - alpha
         if (curv .gt. 0.) then
            fs2 = sqrt(rmakk*curv/rmikk)**3/(shear**2)
         else
            fs2 = 0.
         endif
      else
         sa = alpha - shear
         if (curv .lt. 0.) then
            fs2 = sqrt(-rmakk*curv/rmikk)**3/(shear**2)
         else
            fs2 = 0.
         endif
      endif

      if (sa .ge. 0.0_R8) then
         fs1 = (1.+9.*sqrt(2.0)*sa**2.5)                     &
     &        /(sqrt(2.0)*(1.0-2.0*sa+3.0*sa**2+2.0*sa**3))
      else
         fs1 = 1.0_R8/sqrt(2.0*(1.0-2.0*sa)*(1.0-2.0*sa+3.0*sa**2))
      endif
!     
!     trcofs = fs1
      trcofs = max(fs1,fs2)

!     return
      END FUNCTION trcofs

      REAL FUNCTION cexb(shear,alpha)

!	ExB shearing effect - weak shear approximation

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 alpha,shear,sa
      REAL*8 fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9
      
      sa = shear -alpha
      
      fac2 = (1.0_R8-2.0_R8*sa)
      fac3 = 6.0_R8*sa**2/(1.0_R8-2.0_R8*sa)
      fac1 = 8.0_R8*alpha*fac2/(2.0_R8 + fac3)
!
      fac4 = 3.0_R8/8.0_R8*(1.0_R8-sa**2)
      fac5 = (1.0_R8+6.0_R8*sa**2)/2.0_R8/alpha
!      
      fac6 = (1.0_R8-5.0_R8*sa**2)/4.0_R8
      fac7 = 2.0_R8*sa**2/alpha
      fac8 = (fac6+fac7)**2 + 3.0_R8/8.0_R8      
      fac9 = fac4+fac5 + 0.5_R8*fac8
      cexb = fac1 * fac9
      
      END FUNCTION cexb

      REAL FUNCTION fexb(wexb,shear,alpha)

!	ExB shearing effect - strong shear approximation
!	WARNING: this part has not been properly debugged for lack of test cases   
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 alpha,wexb,shear,alpha1,alpha2,beta,gamma
      

      if(abs(alpha) .lt. 1.e-3_R8) then
         alpha1 = 1.0e-3_R8
      else
         alpha1 = abs(alpha)
      endif

      beta = 0.5_R8*alpha1**(-0.602_R8)                                        &
     & *(13.018_R8-22.28915_R8*shear+17.018_R8*shear**2)                          &
     & /(1._R8-0.277584_R8*shear+1.42913_R8*shear**2)

      alpha2 = (-10.0_R8/3.0_R8)*alpha+(16.0_R8/3.0_R8)

      if(shear .lt. 0.) then
      gamma = (1.0_R8/(1.1_R8*sqrt(1.-shear-2.*shear**2-3.0_R8*shear**3)))+0.75_R8
      else
      gamma = ((1.0_R8-0.5_R8*shear)                                            &
     & /(1.1_R8-2.0_R8*shear+alpha2*shear**2+4.0_R8*shear**3))+0.75_R8
      endif

      fexb = exp(-beta*wexb**gamma)

!     return
      END FUNCTION fexb

 END SUBROUTINE cdbm_chi


