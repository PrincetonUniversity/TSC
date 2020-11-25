      function cexb(shear,alpha)


!	ExB shearing effect for te CDBM model

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 cexb,alpha,shear,sa
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
      
      return
      end
