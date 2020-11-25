      function fexb(wexb,shear,alpha)


!	ExB shearing effect for te CDBM model

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 fexb,alpha,wexb,shear,alpha1,alpha2,beta,gamma
      

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

      return
      end
