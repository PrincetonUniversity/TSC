      subroutine cdbm_dErdr(dpidr,vrotNB,vtori,vpoli,Erad,wExB)
      
      USE CLINAM
      USE SAPROP

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      
      REAL*8, PARAMETER :: aee = 1.602176487e-19

      INTEGER i,j,isw

      REAL*8, DIMENSION(npsit) :: dpidr,dErdr,bpol,vrotNB,Erad,wExB
      REAL*8, DIMENSION(npsit) :: Erads,Eraf,Erfacs,dErdrs,vtori,vpoli
      REAL*8 delta2,rhodens,va,eps0,wpe2,factor
      REAL*8 dum_br,dum_bz,rint,rmaj,pval1,ppval1,eval1,peval1
      REAL*8, DIMENSION(nx) :: bpol_dum
      REAL*8 gradsq,dpsidx,dpsidz,gsval,psval,psixz,psixx,psizz
!  
!      if(kk .eq. 1) rmajora(1) = xmag      
      if(zeff .le. 1.) zeff = 1.5
!
!      do ii=1,nx
!         call grap(3,zmag,xary(ii),gradsq,dpsidx,dpsidz,gsval,psval,    &
!     &              psixz,psixx,psizz,isw)
!         dum_br = dpsidz/xary(ii)
!         dum_bz = -dpsidx/xary(ii)
!         bpol_dum(ii) = sqrt(dum_br**2+dum_bz**2)
!      end do
!      do kk=1,npsit
!         do ii=1,nx
!            if (xsv(kk).gt.psi(ii,jmag).and.xsv(kk).le.psi(ii+1,jmag)) then
!               rint=(xsv(kk)-psi(ii,jmag))/(psi(ii+1,jmag)-psi(ii,jmag))
!               dum_br=bpol_dum(ii)+(bpol_dum(ii+1)-bpol_dum(ii))*rint
!            end if
!         end do
!         bpol(kk)=dum_br
!      end do
!
!     rmajora(1) = xmag
      do j=1,npsit         
        rmaj = xmag+rminora(j)
        if(j .le. 2) then
           bpol(j) = -(-xsv2(j+2)+4.*xsv2(j+1)-3.*xsv2(j))              &
     &                 / ((rminora(j+2)-rminora(j))*rmaj)
        elseif (j .ge. 2 .and. j .le. npsit-2) then
           bpol(j) = -(xsv2(j+2)+8.*xsv2(j+1)-8.*xsv2(j-1)+xsv2(j-2))   &
     &                 / (3.*(rminora(j+2)-rminora(j-2))*rmaj)
        elseif (j .gt. npsit-2) then
           bpol(j) = -(3.*xsv2(j)-4.*xsv2(j-1)+xsv2(j-2))               &
     &                 / ((rminora(j)-rminora(j-2))*rmaj)
        end if
      end do
!
      do j=2,npsit
        call peval(xsv2(j),2,pval1,ppval1,imag,jmag)
        call eeval(xsv2(j),2,eval1,peval1,imag,jmag)
         Erads(j) = dpidr(j)*ti(j)*usdh/(pval1-eval1)/(aee*avez(j))     &
     &            + bpol(j)*(vtori(j)+vrotNB(j))                        &
     &            - vpoli(j)*gzero/rmajora(j)
!        Erfacs(j) = Erads(j)/(rmajora(j)*bpol(j))
         Erfacs(j) = Erads(j)/(rminora(j)*bpol(j))
      end do
!     Erads(1) = 0.0_R8
!     Erfacs(1)= 0.0_R8
!     smoothen Erad profile    
      Erad(2)=(Erads(1)+Erads(2)+Erads(3)) / 3._R8
      Eraf(2)=(Erfacs(1)+Erfacs(2)+Erfacs(3)) / 3._R8
      Erad(1)=Erad(2)
      Eraf(1)=Eraf(2)
      do j=3,npsit-2
         Erad(j) = (Erads(j-2)+Erads(j-1)+Erads(j)                      &
     &           + Erads(j+1)+Erads(j+2))/5.
         Eraf(j) = (Erfacs(j-2)+Erfacs(j-1)+Erfacs(j)                  &
     &           + Erfacs(j+1)+Erfacs(j+2))/5.
      end do
      Erad(npsit-1)=(Erads(npsit-2)+Erads(npsit-1)+Erads(npsit))/3._R8
      Erad(npsit)=Erads(npsit)
      Eraf(npsit-1)=(Erfacs(npsit-2)+Erfacs(npsit-1)+Erfacs(npsit))/3._R8
      Eraf(npsit)=Erfacs(npsit)
!
      do j=1,npsit
        if(j .le. 2) then
           dErdrs(j) = (-Eraf(j+2)+4.0*Eraf(j+1)-3.0*Eraf(j))          &
     &              /(rminora(j+2)-rminora(j))
        elseif (j .gt. 2 .and. j .le. npsit-2) then
           dErdrs(j) = (-Eraf(j+2)+8.*Eraf(j+1)-8.*Eraf(j-1)+Eraf(j-2))  &
     &                   /(3.*(rminora(j+2)-rminora(j-2)))
        elseif (j .gt. npsit-2) then
           dErdrs(j)=(3.*Eraf(j)-4.*Eraf(j-1)+Eraf(j-2))                &
     &               /(rminora(j)-rminora(j-2))
        endif
      end do
!     dErdrs(1) = 0.0_R8
!      
      dErdr(1)=dErdrs(1)
      dErdr(2)=(dErdrs(1)+dErdrs(2)+dErdrs(3))/3.
      do j=3,npsit-2
         dErdr(j)=(dErdrs(j-2)+dErdrs(j-1)+dErdrs(j)+dErdrs(j+1)        &
     &                  +dErdrs(j+2))/5.
      end do
!
      wExB(1) = 0.0_R8
      do j = 2,npsit
         wExB(j) = dErdr(j)*bpol(j)*rmajora(j)**2/gzero
      enddo


      return
      end

