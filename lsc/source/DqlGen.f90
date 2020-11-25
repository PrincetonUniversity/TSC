!     -----------------------------------------------------------------
      SUBROUTINE DqlGen
!     deposit wave energy into dql matrix
      USE DqlBins
      USE FeBins
      USE params
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iv, ip, iNotSmoo, iYesSmoo
      INTEGER iry, izn
      DATA iNotSmoo, iYesSmoo / 1 , 2 /
!     deposit energy density on ray iray and zone izone
!     into appropriate
!     velocity bin.
!     See eq 22 and 28 of Valeo/Eder.
!     D_{ql} = \pi/2 (e/m)^2 Sum_j <E_z_j^2> \delta(\omega - k_par v_par)
!     where the j represents all rays and intersections with the flux surface
!     in question (index supressed here) and the <> mean flux surface avg.
!     The factor 2 is appropriate for E meaning an amplitude; it makes the
!     E^2 into an rms quantity.
!     The delta function connects the k-par dependence of E with the v
!     dependence of D.
!     \delta(w) = (2\pi h)^{-.5} exp(-(w/h)^2/2.) where h is a variable
!     width.  By writing in terms of velocity index and nsmoo
!
!     D(v_i) = \pi/2 (e/m)^2 Sum_j E_j^2 (v_par/(\omega \Delta v_par) *
!         (2\pi nsmoo)^{-.5} exp (- (i-j)^2/(2 nsmoo^2))
!
!     .                                 Clear the unsmoothed part of Dql
      call DqlClear
!     .                                 Fill up the unsmoothed part of Dql
      do 20 iry  = 1, nrays
 
         do 10 izn   = 1, nzones-1
            iv = ivind(izn  , iry )
!                                       If  ivind=0, the ray was stopped before
!                                       reaching this izone (izn).
!                                       If Power =0, no contribution to Dql
!                                       No calculation is appropriate.
            if (             iv .eq. 0   .or.                            &  
     &          Power(izn,iry)  .eq. 0.00_R8) go to 11
            ip = izind(izn  , iry )
!
!     We take the peak power as equal to the average power.
!
!           if ( iFudgDmp(iry) .eq. izn) then
!
!              Dql(iv, ip, iNotSmoo) = Dql(iv, ip, iNotSmoo) +
!     ^            Power(izn  , iry ) * rFudgDmp(iry)
!     ^           * Ezsq(izn  , iry ) * abs(vpar(iv)) / dvsym(iv)
!            else
!              Dql(iv, ip, iNotSmoo) = Dql(iv, ip, iNotSmoo) +
!     ^            Power(izn  , iry )
!     ^           * Ezsq(izn  , iry ) * abs(vpar(iv)) / dvsym(iv)
!            endif
 
              Dql(iv, ip, iNotSmoo) = Dql(iv, ip, iNotSmoo) +            &  
     &            Power(izn  , iry ) * rFudgDmp(izn,iry)                 &  
     &           * Ezsq(izn  , iry ) * abs(vpar(iv)) / dvsym(iv)
!
 10      continue
 11   continue
 20   continue
!
!     .                                 Multiply by normalization
      call svmult(NVELDIM * NPSIDIM, Dql, Dqlnorm, Dql)
!
!     .                                 Copy unsmoothed into space for smoothed
      do 30 iv=1,NVELDIM
        do 25 ip=1,NPSIDIM
          Dql(iv,ip,iYesSmoo) = Dql(iv,ip,iNotSmoo)
 25     continue
 30   continue
!
!     .                                 Smooth the new copy
      call Smooth(Dql(1,1,iYesSmoo),                                     &  
     &            nv, NVELDIM, npsi, nsmoo, qlsm)
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
