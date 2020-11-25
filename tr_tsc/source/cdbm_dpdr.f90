      subroutine cdbm_dpdr(pvalt,dpdr,dpidr,dptdr)
      
      USE CLINAM
      USE SAPROP

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      
      REAL*8, PARAMETER :: aee = 1.602176487e-19

      INTEGER i,j,model,npsihm,n_int

      REAL*8 pval1,eval1,ppval1,peval1
      REAL*8 pvala,pvalb,pvalc
      REAL*8 pvalia,pvalib,pvalic
      REAL*8, DIMENSION(npsit) :: dpdr,dpidr,dptdr,pval,pvali,pvalt
      REAL*8, DIMENSION(npsit) :: dpdrs,dpidrs,dptdrs,pvals,pvalis,pvalts 
! 
      do j=1,npsit
         call peval(xsv2(j),2,pval1,ppval1,imag,jmag)
         call eeval(xsv2(j),2,eval1,peval1,imag,jmag)
         pval(j) = pval1*udsp
         pvali(j) = (pval1-eval1)*udsp
         pvalt(j) = pvalt(j)*udsp
      end do
!

! apply smoothing function to pressure profiles to reduce oscillations in chi
         pvals(1) = pval(1)
         pvalis(1) = pvali(1)
         pvalts(1) = pvalt(1)
         pvals(2)=(pval(1)+pval(2)+pval(3))/3._R8
         pvalis(2)=(pvali(1)+pvali(2)+pvali(3))/3._R8
         pvalts(2)=(pvalt(1)+pvalt(2)+pvalt(3))/3._R8
      do j=3,npsit-2
         pvals(j)=(pval(j-2)+pval(j-1)+pval(j)+pval(j+1)+pval(j+2))/5._R8
         pvalis(j)=(pvali(j-2)+pvali(j-1)+pvali(j)                      &
     &             + pvali(j+1)+pvali(j+2))/5._R8
         pvalts(j)=(pvalt(j-2)+pvalt(j-1)+pvalt(j)                      &
     &             + pvalt(j+1)+pvalt(j+2))/5._R8
      end do
!
! calculate thermal, ion and total  pressure gradients
!     dpdrs(1) = 0.0
!     dpidrs(1)=0.0
!     dptdrs(1)=0.0
      do j=1,npsit
         if(j .le. 2) then
            dpdrs(j) = (-pvals(j+2)+4.0*pvals(j+1)-3.0*pvals(j))        &
     &                  /(rminora(j+2)-rminora(j))
            dpidrs(j) = (-pvalis(j+2)+4.0*pvalis(j+1)-3.0*pvalis(j))    &
     &                  /(rminora(j+2)-rminora(j))
            dptdrs(j) = (-pvalts(j+2)+4.0*pvalts(j+1)-3.0*pvalts(j))    &
     &                  /(rminora(j+2)-rminora(j))
         elseif (j .gt. 2 .and. j .le. npsit-2) then
            dpdrs(j)=(-pvals(j+2)+8.*pvals(j+1)-8.*pvals(j-1)+pvals(j-2))&
     &                   / (3.0*(rminora(j+2)-rminora(j-2)))
            dpidrs(j)=(-pvalis(j+2)+8.*pvalis(j+1)-8.*pvalis(j-1)       &
     &                + pvalis(j-2))/(3.0*(rminora(j+2)-rminora(j-2)))
            dptdrs(j)=(-pvalts(j+2)+8.*pvalts(j+1)-8.*pvalts(j-1)       &
     &                + pvalts(j-2))/(3.0*(rminora(j+2)-rminora(j-2)))
         elseif (j .gt. npsit-2) then
            dpdrs(j) = (3.*pvals(j)-4.*pvals(j-1)+pvals(j-2))           &
     &                  / (rminora(j)-rminora(j-2))
            dpidrs(j) = (3.*pvalis(j)-4.*pvalis(j-1)+pvalis(j-2))       &
     &                  / (rminora(j)-rminora(j-2))
            dptdrs(j) = (3.*pvalts(j)-4.*pvalts(j-1)+pvalts(j-2))       &
     &                  / (rminora(j)-rminora(j-2))
         end if
      end do
!
! apply smoothing function to pressure gradient
      dpdr(1) = dpdrs(1)
      dpidr(1) = dpidrs(1)
      dptdr(1) = dptdrs(1)
      dpdr(2) = (dpdrs(1)+dpdrs(2)+dpdrs(3)) /3.0
      dpidr(2) = (dpidrs(1)+dpidrs(2)+dpidrs(3)) /3.0
      dptdr(2) = (dptdrs(1)+dptdrs(2)+dptdrs(3)) /3.0
!     dpdr(1) = dpdr(2)
!     dpidr(1) = dpidr(2)
!     dptdr(1) = dptdr(2)
      do j=3,npsit-2
         dpdr(j) = (dpdrs(j-2)+dpdrs(j-1)+dpdrs(j)                      &
     &              + dpdrs(j+1)+dpdrs(j+2)) /5.0
         dpidr(j) = (dpidrs(j-2)+dpidrs(j-1)+dpidrs(j)                  &
     &               + dpidrs(j+1)+dpidrs(j+2)) / 5.0
         dptdr(j) = (dptdrs(j-2)+dptdrs(j-1)+dptdrs(j)                  &
     &               + dptdrs(j+1)+dptdrs(j+2)) / 5.0
      end do
      dpdr(npsit) = dpdrs(npsit)
      dpidr(npsit) = dpidrs(npsit)
      dptdr(npsit) = dptdrs(npsit)
      dpdr(npsit-1) = (dpdrs(npsit-2)+dpdrs(npsit-1)+dpdrs(npsit)) /3.0
      dpidr(npsit-1) = (dpidrs(npsit-2)+dpidrs(npsit-1)+dpidrs(npsit)) /3.0
      dptdr(npsit-1) = (dptdrs(npsit-2)+dptdrs(npsit-1)+dptdrs(npsit)) /3.0
!
      return
      end

