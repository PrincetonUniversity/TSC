      function trcofs(shear,alpha,curv,rmakk,rmikk)

     IMPLICIT NONE
      
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 trcofs,alpha,curv,shear,sa,fs1,fs2,rmakk,rmikk


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

      return
      end


