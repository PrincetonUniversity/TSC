!
!     -----------------------------------------------------------------
!
      SUBROUTINE FastFrac(npsi)
      USE FeBins
      USE params
      USE PlPr
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER npsi, i,j, jwt
      REAL*8    ExpMax, RsltMin
      REAL*8    duNorm, duVth2, duFracN, duFracE, duMaxwN, duMaxwE
      REAL*8    dumFe, dumMx, dumV2, duDelV, exp1
      REAL*8    wt(2)
      DATA    wt(1)         , wt(2)          /                           &  
     &        0.666666666666_R8, 1.333333333333_R8/
      DATA    ExpMax / 100._R8/
!     FastFrac: Fast Fraction of N FstFracN and of Energy FstFracE
!     are computed and writtend in Subroutine FastFrac(npsi)
!     du        prefix meaning dummy
!     duNorm    Normalization constant for this flux surface
!               such that f_e(v) = duNorm * exp(-v^2/(2 v_T^2) for
!               Maxwellian
!     duVth2    v_T^2 normalizing the exponent for this flux surface
!     duDelV    The width of the velocity grid at the center
!     duFracN   The fraction of particles in the heated portion of the tail
!     duFracE   The fraction of energy in the heated portion of the tail
!     duMaxwN   "particles" in Maxwellian; (local units not relevant)
!     duMaxwE   "energy" in Maxwellian; (local units not relevant)
!     FeNorm    array dimensioned NPSIDIM containing NeAry/(2 PI)/Vtherm(ip)
!
!     Simpsons Rule is (f1 + 4f2 + f3)h/3 or
!                     (2f1 + 4f2 + 2f3 + 4f4 + 2f5 )h/3 less
!                   - ( f1                      f5 )h/3 which can be ignored
!     because the end points are exponentially small.  The weight (wt) is
!     from the expression  [ 2 + 2( j / 2 ) - j ]  = 1, {j=1,3,5...}
!                                                  = 2, {j=2,4,6...}
!
      duDelV  = abs ( Vpar(IvZero+1) - Vpar(IvZero) )
      RsltMin = exp ( -ExpMax )
      do 20 i = 1, npsi
        duFracN = 0._R8
        duFracE = 0._R8
        duMaxwN = 0._R8
        duMaxwE = 0._R8
        duNorm  = FeNorm(i)
        duVth2  = Vtherm(i)*Vtherm(i)
        if ( vtherm(i) .ge. duDelV ) then
          do 10 j = 1, nv
            jwt =  ( 2 + 2 * ( j / 2 ) - j )
            dumFe  =              fe(j,i,iITR) * wt(jwt) * dvsym(j)
            dumV2  = vpar(j)**2 * fe(j,i,iITR) * wt(jwt) * dvsym(j)
            duFracN = duFracN + dumFe
            duFracE = duFracE + dumV2
 
            exp1 =  vpar(j)*vpar(j)/(2._R8*duVth2)
            if(exp1 .lt. ExpMax) then
              exp1 = exp( -exp1 )
            else
              exp1 = RsltMin
            endif
            exp1 = exp1 * duNorm
!           exp1 = exp( -vpar(j)*vpar(j)/(2.*duVth2) )* duNorm
              dumMx  =              exp1 * wt(jwt) * dvsym(j)
              dumV2  = vpar(j)**2 * exp1 * wt(jwt) * dvsym(j)
              duMaxwN = duMaxwN + dumMx
              duMaxwE = duMaxwE + dumV2
 10       continue
 
          FstFracN(i) = ( duFracN - duMaxwN ) / duMaxwN
            FstFracE(i) = ( duFracE - duMaxwE ) / duMaxwE
        else
          FstFracN(i) = 0._R8
          FstFracE(i) = 0._R8
        endif
 20   continue
      if (PrFlg(FSTFRCWR) .ge. TRUE)  then
!       call LSCpause
        write(nLSCcomm,1000)
 1000       format (                                                     &  
     &'  Psi Index   Fraction Tail Particles  Fraction Tail Energy',/)
         write(nLSCcomm,1001)                                            &  
     &                     (  i, FstFracN(i), FstFracE(i)  , i=1,npsi)
 1001       format(t8, i4, t28, 1pe10.3, t50, 1pe10.3)
!       call LSCpause
      endif
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
