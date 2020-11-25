!     -----------------------------------------------------------------
      SUBROUTINE FePlus(fe, nuCollx, Dcollx, Dqlx, Vpar, lgFe, ivZero,   &  
     &     nv)
!     solve for fe for positive velocity, using the convention that
!     the new value is the old value plus the mean integrand
!     times the delta-v to the new value
      USE params
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nv, ivZero, iv
      REAL*8    ExpMax, RsltMin
!     REAL*8    exp
      REAL*8                                                             &  
     &        fe(nv), nuCollx(nv), lgFe(nv),                             &  
     &        Dcollx(nv), Dqlx(nv), Vpar(nv)
!                                       fe(ivZero) is assumed
!                                       prescribed elsewhere
!     DATA    ExpMax / 100. /           Alpha and Ted like 85 better than 100
      DATA    ExpMax /  85._R8/
      INTEGER DoCut6
      DATA    DoCut6 /0/
      lgFe(ivZero) = 0._R8
      RsltMin = exp ( - ExpMax )
!     RsltMin = dexp ( - ExpMax )
      do 10 iv = ivZero, nv - 1
         lgFe(iv + 1) = lgFe(iv) +                                       &  
     &   0.5_R8*(nuCollx(iv+1) / (Dqlx(iv+1) + Dcollx(iv+1)) +           &  
     &        nuCollx(iv  ) / (Dqlx(iv  ) + Dcollx(iv  )))*              &  
     &        (Vpar(iv + 1) - Vpar(iv))
10    continue
      do 20 iv = ivZero + 1, nv-2
        if( lgFe(iv) .lt. ExpMax) then
          fe(iv) = fe(ivZero) * exp( - lgFe(iv) )
        else
          fe(iv) = fe(ivZero) * RsltMin
        endif
20    continue
      fe(nv-1) = 0.00_R8
      fe(nv  ) = 0.00_R8
      if (DoCut6 .eq. 1) then
         do iv=1,nv
            if (abs(vpar(iv)) .gt. 0.6_R8) fe(iv) = 0.00_R8
         enddo
      endif
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
